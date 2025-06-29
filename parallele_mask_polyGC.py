#!/usr/bin/env python3
"""
Created on Wed Jun 25 14:52:22 2025

@author: carol
"""

import argparse
import multiprocessing
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import re
import sys

def process_single_file(input_file, output_file_path, threshold):
    """
    Processes a single FASTQ/FASTA file to mask poly-G/C sequences.
    """
    pattern = re.compile(f"(G{{{int(threshold)+1},}}|C{{{int(threshold)+1},}})")
    records_to_write = []
    try:
        for record in SeqIO.parse(input_file, "fastq"):
            seq_str = str(record.seq)
            masked_seq = pattern.sub(lambda m: 'N' * len(m.group()), seq_str)
            record.seq = Seq(masked_seq)
            records_to_write.append(record)
        SeqIO.write(records_to_write, output_file_path, "fastq")
    except ValueError:
        # If fastq parsing fails, try fasta
        try:
            for record in SeqIO.parse(input_file, "fasta"):
                seq_str = str(record.seq)
                masked_seq = pattern.sub(lambda m: 'N' * len(m.group()), seq_str)
                record.seq = Seq(masked_seq)
                records_to_write.append(record)
            SeqIO.write(records_to_write, output_file_path, "fasta")
        except Exception as e:
            print(f"Error processing {input_file}: {e}")
            return False
    return True

def mask_poly_gc_parallel(args):
    input_files = []
    for pattern in args.i:
        input_files.extend(glob.glob(pattern))

    if not input_files:
        print("No input files found matching the provided patterns.")
        return

    # Create a temporary directory for intermediate output files
    temp_dir = "temp_masked_output"
    os.makedirs(temp_dir, exist_ok=True)

    # Prepare arguments for parallel processing
    tasks = []
    for i, input_file in enumerate(input_files):
        # Determine output format based on input file extension for temporary files
        if input_file.lower().endswith((".fastq", ".fq")):
            output_format = "fastq"
        elif input_file.lower().endswith((".fasta", ".fa", ".fna")):
            output_format = "fasta"
        else:
            print(f"Skipping {input_file}: Unsupported file extension. Please use .fastq, .fq, .fasta, .fa, or .fna.")
            continue

        temp_output_file = os.path.join(temp_dir, f"masked_part_{i}.{output_format}")
        tasks.append((input_file, temp_output_file, args.t))

    if not tasks:
        print("No valid input files to process.")
        os.rmdir(temp_dir)
        return

    # Run processes in parallel
    with multiprocessing.Pool(processes=args.p) as pool:
        results = pool.starmap(process_single_file, tasks)

    if args.o:
        # Concatenate all temporary output files into a single final output file
        with open(args.o, "w") as outfile:
            for i, (input_file, temp_output_file, _) in enumerate(tasks):
                if results[i] and os.path.exists(temp_output_file):
                    with open(temp_output_file, "r") as infile:
                        outfile.write(infile.read())
                    os.remove(temp_output_file)  # Clean up temporary file
        print(f"Masking complete. All masked sequences written to {args.o}", file=sys.stderr)
    else:
        # Output to stdout
        for i, (input_file, temp_output_file, _) in enumerate(tasks):
            if results[i] and os.path.exists(temp_output_file):
                with open(temp_output_file, "r") as infile:
                    sys.stdout.write(infile.read())
                os.remove(temp_output_file) # Clean up temporary file
        print("Masking complete. All masked sequences written to stdout.", file=sys.stderr)

    os.rmdir(temp_dir) # Remove the temporary directory


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Mask poly-G/C sequences in FASTQ/FASTA files in parallel.")
    parser.add_argument('-i', '--input', dest='i', required=True, nargs='+',
                        help="One or more FASTQ/FASTA files (or patterns like *.fastq) to be masked. Example: 'file1.fastq file2.fastq' or '*.fastq'")
    parser.add_argument('-o', '--output', dest='o', required=False,
                        help="The name of the single output file.")
    parser.add_argument('-t', '--threshold', dest='t', required=True, type=int,
                        help="Threshold to determine polyGC, at least 10.")
    parser.add_argument('-p', '--processes', dest='p', type=int, default=multiprocessing.cpu_count(),
                        help=f"Number of parallel processes to use (default: {multiprocessing.cpu_count()}).")

    args = parser.parse_args()

    if args.t < 10:
        print("Error: The threshold (-t) must be at least 10.")
    else:
        mask_poly_gc_parallel(args)

