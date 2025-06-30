#!/usr/bin/env python3
######################################################################
# parallel_extract_kraken_reads.py
#
# This script is a modified version of extract_kraken_reads.py from
# KrakenTools by Jennifer Lu.
# It has been adapted to process multiple FASTA/FASTQ files in parallel
# against a single Kraken output/report file to improve performance.
#
# Original Copyright (C) 2019-2023 Jennifer Lu, jennifer.lu717@gmail.com
# This file is part of KrakenTools (GPLv3 License)
#
# Modifications by: Google Gemini (2025)
#
# This program extracts reads from multiple sequence files based on
# classifications in a Kraken output file.
######################################################################
import os, sys, argparse, multiprocessing
import gzip
from time import gmtime, strftime
from functools import partial
from Bio import SeqIO

# (The Tree class and other functions from the previous version
# remain unchanged here)

class Tree(object):
    'Tree node.'
    def __init__(self, taxid, level_num, level_id, children=None, parent=None):
        self.taxid = taxid
        self.level_num = level_num
        self.level_id = level_id
        self.children = []
        self.parent = parent
        if children is not None:
            for child in children:
                self.add_child(child)
    def add_child(self, node):
        assert isinstance(node,Tree)
        self.children.append(node)

def process_kraken_report(report_line):
    l_vals = report_line.strip().split('\t')
    if len(l_vals) < 5: return []
    try: int(l_vals[1])
    except ValueError: return []
    try:
        taxid = int(l_vals[-3])
        level_type = l_vals[-2]
        map_kuniq = {'species':'S', 'genus':'G','family':'F',
            'order':'O','class':'C','phylum':'P','superkingdom':'D',
            'kingdom':'K'}
        if level_type not in map_kuniq: level_type = '-'
        else: level_type = map_kuniq[level_type]
    except ValueError:
        taxid = int(l_vals[-2])
        level_type = l_vals[-3]
    spaces = 0
    for char in l_vals[-1]:
        if char == ' ': spaces += 1
        else: break
    level_num = int(spaces/2)
    return[taxid, level_num, level_type]

def process_kraken_output(kraken_line):
    l_vals = kraken_line.split('\t')
    if len(l_vals) < 5: return [-1, '']
    tax_id = l_vals[2]
    if "taxid" in tax_id:
        tax_id = tax_id.split("taxid ")[-1][:-1]
    read_id = l_vals[1]
    if (tax_id == 'A'): tax_id = 81077
    else: tax_id = int(tax_id)
    return [tax_id, read_id]

def extract_reads_worker(task_args, read_ids_to_save, filetype, fastq_out):
    """
    Worker function to extract reads from a single sequence file.
    """
    seq_file, output_file = task_args
    count_seqs = 0
    count_output = 0

    sys.stdout.write(f">> Worker started for: {os.path.basename(seq_file)}\n")
    sys.stdout.flush()

    # Determine if input is gzipped
    if seq_file.endswith('.gz'):
        s_file = gzip.open(seq_file, 'rt')
    else:
        s_file = open(seq_file, 'r')

    # Open output file
    o_file = open(output_file, 'w')

    output_format = "fastq" if fastq_out else "fasta"

    for record in SeqIO.parse(s_file, filetype):
        count_seqs += 1
        test_id = str(record.id)
        # Handle read IDs with suffixes like /1 or /2
        test_id_nosuffix = test_id[:-2] if ("/1" in test_id or "/2" in test_id) else test_id

        if test_id in read_ids_to_save or test_id_nosuffix in read_ids_to_save:
            count_output += 1
            SeqIO.write(record, o_file, output_format)

    s_file.close()
    o_file.close()

    sys.stdout.write(f">> Worker finished for: {os.path.basename(seq_file)}. Extracted {count_output} reads to {output_file}\n")
    sys.stdout.flush()
    return (output_file, count_output)

def main():
    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(description="Parallel extraction of reads from multiple FASTA/FASTQ files based on Kraken classification.")
    parser.add_argument('-k', dest='kraken_file', required=True, help='Kraken output file to parse')
    parser.add_argument('-s', '--seq-files', dest='seq_files', required=True, nargs='+', help='Space-delimited list of FASTA/FASTQ files to process.')
    parser.add_argument('-t', "--taxid", dest='taxid', required=True, nargs='+', help='Taxonomy ID(s) of reads to extract (space-delimited)')
    parser.add_argument('--output-suffix', dest='output_suffix', default=".extracted.fq", help='Suffix to append to input filenames for output files [default: .extracted.fq]')
    parser.add_argument('--out-dir', dest='out_dir', default='.', help='Destination folder to store output files [default: current directory]')
    parser.add_argument('-r','--report', dest='report_file', default="", help='Kraken report file. [required only if --include-parents/children is specified]')
    parser.add_argument('--include-parents', dest="parents", action='store_true', default=False, help='Include reads classified at parent levels of the specified taxids')
    parser.add_argument('--include-children', dest='children', action='store_true', default=False, help='Include reads classified more specifically than the specified taxids')
    parser.add_argument('--exclude', dest='exclude', action='store_true', default=False, help='Extract all reads NOT matching specified taxids')
    parser.add_argument('--fastq-output', dest='fastq_out', action='store_true', default=False, help='Print output in FASTQ format [requires input FASTQ, default: FASTA]')
    parser.add_argument('--max', dest='max_reads', default=100000000, type=int, help='Maximum number of reads to save [default: 100,000,000]')
    parser.add_argument('--threads', dest='threads', type=int, default=os.cpu_count(), help=f'Number of parallel processes to use [default: all available CPUs ({os.cpu_count()})]')

    args=parser.parse_args()

    sys.stdout.write(f"PROGRAM START TIME: {strftime('%m-%d-%Y %H:%M:%S', gmtime())}\n")

    # --- Create output directory if it doesn't exist ---
    os.makedirs(args.out_dir, exist_ok=True)
    sys.stdout.write(f">> Output files will be saved in: {os.path.abspath(args.out_dir)}\n")


    # --- Step 0: Identify all Taxonomy IDs to be included ---
    save_taxids = {int(tid) for tid in args.taxid}
    if args.parents or args.children:
        if not args.report_file:
            sys.stderr.write(">> ERROR: --report must be specified with --include-parents/children.\n")
            sys.exit(1)
        sys.stdout.write(f">> STEP 0: PARSING REPORT FILE {args.report_file}\n")
        # (The tree-building logic from the original script remains unchanged)
        base_nodes = {}
        with open(args.report_file, 'r') as r_file:
            prev_node = None
            for line in r_file:
                report_vals = process_kraken_report(line)
                if not report_vals: continue
                [taxid, level_num, level_id] = report_vals
                if taxid == 0: continue
                if taxid == 1:
                    level_id = 'R'
                    root_node = Tree(taxid, level_num, level_id)
                    prev_node = root_node
                    if taxid in save_taxids: base_nodes[taxid] = root_node
                    continue
                while level_num != (prev_node.level_num + 1):
                    prev_node = prev_node.parent
                curr_node = Tree(taxid, level_num, level_id, None, prev_node)
                prev_node.add_child(curr_node)
                prev_node = curr_node
                if taxid in save_taxids: base_nodes[taxid] = curr_node
        if args.parents:
            for tid in base_nodes:
                curr_node = base_nodes[tid]
                while curr_node.parent != None:
                    curr_node = curr_node.parent
                    save_taxids.add(curr_node.taxid)
        if args.children:
            for tid in base_nodes:
                curr_nodes = list(base_nodes[tid].children)
                while curr_nodes:
                    curr_n = curr_nodes.pop()
                    save_taxids.add(curr_n.taxid)
                    if curr_n.children:
                        curr_nodes.extend(curr_n.children)

    # --- Step 1: Parse Kraken File to get all Read IDs ---
    sys.stdout.write(f"\tFound {len(save_taxids)} taxonomy IDs to parse\n")
    sys.stdout.write(f">> STEP 1: PARSING KRAKEN FILE FOR READ IDS: {args.kraken_file}\n")

    read_ids_to_save = set()
    exclude_taxids = save_taxids if args.exclude else set()
    if args.exclude: save_taxids = set()

    with open(args.kraken_file, 'r') as k_file:
        for i, line in enumerate(k_file):
            if (i + 1) % 1000000 == 0:
                sys.stdout.write(f'\r\t{ (i + 1) / 1000000:.2f} million kraken lines processed...')
                sys.stdout.flush()
            [tax_id, read_id] = process_kraken_output(line)
            if tax_id == -1: continue

            if args.exclude:
                if tax_id not in exclude_taxids:
                    read_ids_to_save.add(read_id)
            elif tax_id in save_taxids:
                read_ids_to_save.add(read_id)

            if len(read_ids_to_save) >= args.max_reads:
                break
    sys.stdout.write(f'\r\tFinished processing kraken file. Found {len(read_ids_to_save)} total unique read IDs to extract.\n')

    if not read_ids_to_save:
        sys.stdout.write(">> No reads found for the specified taxid(s). Exiting.\n")
        sys.exit(0)

    # --- Step 2: Determine file type from first input file ---
    first_file = args.seq_files[0]
    if first_file.endswith('.gz'):
        with gzip.open(first_file, 'rt') as f:
            first_char = f.read(1)
    else:
        with open(first_file, 'r') as f:
            first_char = f.read(1)

    if first_char == ">": filetype = "fasta"
    elif first_char == "@": filetype = "fastq"
    else:
        sys.stderr.write("ERROR: Sequence files must be in FASTA or FASTQ format.\n")
        sys.exit(1)

    if filetype != 'fastq' and args.fastq_out:
        sys.stderr.write('ERROR: for FASTQ output, input files must be FASTQ\n')
        sys.exit(1)

    # --- Step 3: Create job list and process in parallel ---
    tasks = []
    for seq_file in args.seq_files:
        if not os.path.exists(seq_file):
            sys.stderr.write(f"WARNING: Input file not found, skipping: {seq_file}\n")
            continue
        base, _ = os.path.splitext(os.path.basename(seq_file))
        if seq_file.endswith('.gz'):
             base, _ = os.path.splitext(base) # handle .fastq.gz
        
        output_filename = f"{base}{args.output_suffix}"
        # *** Construct full output path using the destination directory ***
        output_path = os.path.join(args.out_dir, output_filename)
        tasks.append((seq_file, output_path))


    sys.stdout.write(f">> STEP 2: Spawning {args.threads} workers to process {len(tasks)} files in parallel...\n")

    # Use functools.partial to create a function with fixed arguments
    worker_func = partial(extract_reads_worker,
                          read_ids_to_save=read_ids_to_save,
                          filetype=filetype,
                          fastq_out=args.fastq_out)

    with multiprocessing.Pool(processes=args.threads) as pool:
        results = pool.map(worker_func, tasks)

    # --- Finalize ---
    sys.stdout.write("\n--- SUMMARY ---\n")
    total_reads = 0
    for output_file, count in results:
        total_reads += count
        sys.stdout.write(f"\tGenerated file: {output_file} ({count} reads)\n")

    sys.stdout.write(f"\nTotal reads extracted across all files: {total_reads}\n")
    sys.stdout.write(f"PROGRAM END TIME: {strftime('%m-%d-%Y %H:%M:%S', gmtime())}\n")
    sys.exit(0)

if __name__ == "__main__":
    main()
    
