#!/usr/bin/env python3
"""
Created on Sun Jun 29 15:38:16 2025

@author: Google Gemini (2025)
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline # For running local blastn
from Bio.Blast import NCBIXML # For parsing XML output
import argparse
import sys
import os # For temporary files

def find_and_reorient_dnaA(assembly_file, dnaA_sequence_file, output_file):
    """
    Finds the dnaA sequence in a circular assembly using blastn (query vs subject),
    reorients the assembly to start at the dnaA sequence,
    and ensures dnaA is on the forward strand based on the best BLAST hit.

    Args:
        assembly_file (str): Path to the input assembly file (FASTA or GenBank).
                             This file will be used as the -subject for blastn.
                             (Assumes this file contains a single sequence record).
        dnaA_sequence_file (str): Path to the file containing the dnaA DNA sequence (FASTA).
                                  This file will be used as the -query for blastn.
        output_file (str): Path for the output modified assembly file.

    Returns:
        Bio.SeqRecord or None: The modified SeqRecord if successful, None otherwise.
    """
    # --- Input Validation and Reading ---
    try:
        # Read the assembly. Assuming a single record in the file.
        record = SeqIO.read(assembly_file, "fasta") # <<< IMPORTANT: Change "fasta" to "genbank" if your input assembly is GenBank
    except ValueError as e:
        print(f"Error reading assembly file '{assembly_file}': {e}. Make sure it contains a single record and the format is correct.", file=sys.stderr)
        return None
    except FileNotFoundError:
        print(f"Error: Assembly file '{assembly_file}' not found.", file=sys.stderr)
        return None

    try:
        # Read the dnaA sequence from its file. Assuming a single record in the file.
        dnaA_seq_record = SeqIO.read(dnaA_sequence_file, "fasta") # dnaA sequence is typically in FASTA
        dnaA_seq = dnaA_seq_record.seq
    except ValueError as e:
        print(f"Error reading dnaA sequence file '{dnaA_sequence_file}': {e}. Make sure it contains a single sequence.", file=sys.stderr)
        return None
    except FileNotFoundError:
        print(f"Error: dnaA sequence file '{dnaA_sequence_file}' not found.", file=sys.stderr)
        return None

    original_seq = record.seq
    genome_len = len(original_seq)

    # --- Run blastn with -subject ---
    query_fasta = dnaA_sequence_file # Use the provided dnaA file as query
    subject_fasta = assembly_file    # Use the provided assembly file as subject
    blast_output_xml = "blastn_results.xml" # Temporary output file

    print(f"Running blastn with '{query_fasta}' as query and '{subject_fasta}' as subject...")
    try:
        # Define BLAST command line. -outfmt 5 for XML output (easily parsable by Biopython)
        # -max_target_seqs 1 to get only the best match
        blastn_cline = NcbiblastnCommandline(
            query=query_fasta,
            subject=subject_fasta, # Use -subject instead of -db
            outfmt=5, # XML output
            out=blast_output_xml,
            max_target_seqs=1, # We only need the top hit for reorientation
            evalue=0.001 # Set a reasonable E-value cutoff
        )
        stdout, stderr = blastn_cline() # Execute the command
        
        if stderr:
            print(f"BLASTn warnings/errors: {stderr}", file=sys.stderr)
        print("BLASTn search completed.")

    except Exception as e:
        print(f"Error running blastn: {e}. Make sure 'blastn' is in your PATH.", file=sys.stderr)
        if os.path.exists(blast_output_xml):
            os.remove(blast_output_xml)
        return None

    # --- Parse BLAST Results ---
    best_hit_start = -1
    best_hit_strand = 0 # 1 for forward, -1 for reverse
    hsp = None # Initialize hsp to be accessible later if found

    try:
        with open(blast_output_xml, "r") as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            
            for blast_record in blast_records: # Expect only one record as one query
                if blast_record.alignments:
                    alignment = blast_record.alignments[0] # Top alignment
                    
                    if alignment.hsps:
                        hsp = alignment.hsps[0] # Top HSP
                        
                        # Determine the start position on the assembly (0-based)
                        # and INFER THE STRAND based on start/end coordinates
                        if hsp.sbjct_start < hsp.sbjct_end: # Alignment on forward strand of subject
                            best_hit_start = hsp.sbjct_start - 1 # Convert to 0-based
                            best_hit_strand = 1 # Set to forward strand
                            inferred_strand_label = "Plus" # For printing
                        else: # Alignment on reverse strand of subject
                            best_hit_start = hsp.sbjct_end - 1 # Convert to 0-based (hsp.sbjct_end is smaller for reverse strand)
                            best_hit_strand = -1 # Set to reverse strand
                            inferred_strand_label = "Minus" # For printing
                        
                        print(f"Best BLAST hit found:")
                        print(f"  Query ID: {blast_record.query_id}")
                        print(f"  Subject ID: {alignment.title}")
                        print(f"  E-value: {hsp.expect}")
                        print(f"  Identity: {hsp.identities}/{hsp.align_length} ({hsp.identities/hsp.align_length:.2%})")
                        print(f"  Subject start: {hsp.sbjct_start}, Subject end: {hsp.sbjct_end}")
                        # Use the inferred strand label for printing instead of hsp.sbjct_strand
                        print(f"  Subject strand (inferred): {inferred_strand_label}")
                        break # Processed the best hit, exit loop
                    else:
                        print("No HSPs found for the best alignment.", file=sys.stderr)
                else:
                    print("No alignments found for the dnaA sequence. dnaA might not be in the assembly or too dissimilar.", file=sys.stderr)
                break # Processed the first (and likely only) BLAST record
            
    except Exception as e:
        print(f"Error parsing BLAST results: {e}", file=sys.stderr)
        best_hit_start = -1 # Indicate failure
    
    finally:
        # Clean up temporary BLAST output file
        if os.path.exists(blast_output_xml):
            os.remove(blast_output_xml)
        print("Cleaned up temporary BLAST output file.")



    # --- Reorientation Logic ---
    modified_record = None

    if best_hit_start != -1:
        rotation_point = best_hit_start
        
        if best_hit_strand == 1: # dnaA found on forward strand of assembly
            print(f"Reorienting assembly to start at position {rotation_point} on forward strand.")
            reoriented_seq = original_seq[rotation_point:] + original_seq[:rotation_point]
            modified_record = SeqRecord(
                reoriented_seq,
                id=record.id + "_reoriented",
                name=record.name,
                description=record.description + " (reoriented to dnaA start)"
            )
            
            # --- Feature Reorientation for Forward Strand ---
            if hasattr(record, 'features') and record.features:
                print("Attempting to reorient features (forward strand).")
                new_features = []
                for feature in record.features:
                    try:
                        new_start_pos = (feature.location.start - rotation_point) % genome_len
                        new_end_pos = (feature.location.end - rotation_point) % genome_len
                        
                        new_location = feature.location.__class__(
                            new_start_pos, new_end_pos, strand=feature.location.strand,
                            circular=feature.location.circular if hasattr(feature.location, 'circular') else False
                        )
                        new_features.append(feature.__class__(
                            location=new_location, type=feature.type, qualifiers=feature.qualifiers
                        ))
                    except Exception as fe:
                        print(f"Warning: Could not reorient feature '{feature.type}' at {feature.location}. Error: {fe}", file=sys.stderr)
                modified_record.features = new_features

        elif best_hit_strand == -1: # dnaA found on reverse strand of assembly
            print(f"dnaA found on reverse strand. Reverse complementing and reorienting assembly to start at position {rotation_point}.")
            
            # First, reverse complement the entire assembly
            rc_assembly_seq = original_seq.reverse_complement()
            
            # For a reverse hit, hsp.sbjct_start is the larger coordinate (the "end" of the feature)
            # and hsp.sbjct_end is the smaller coordinate (the "start" of the feature)
            # relative to the original sequence.
            # When you reverse complement the *entire* genome, the base at `X` moves to `genome_len - 1 - X`.
            # So, the original `hsp.sbjct_start` (1-based, larger coordinate) moves to `genome_len - hsp.sbjct_start` (0-based)
            # This will be the new start of the dnaA on the *reverse-complemented* sequence.
            
            rotation_point_on_rc = genome_len - hsp.sbjct_start # Convert 1-based hsp.sbjct_start (larger coord) to 0-based on RC sequence
            
            reoriented_seq = rc_assembly_seq[rotation_point_on_rc:] + rc_assembly_seq[:rotation_point_on_rc]

            modified_record = SeqRecord(
                reoriented_seq,
                id=record.id + "_reoriented_rc",
                name=record.name,
                description=record.description + " (reoriented to dnaA start, strand reversed)"
            )
            
            # --- Feature Reorientation for Reverse Complemented + Rotated ---
            if hasattr(record, 'features') and record.features:
                print("Attempting to reorient and reverse-complement features.")
                new_features = []
                for feature in record.features:
                    try:
                        # 1. Reverse complement the feature's location relative to the original sequence length
                        rc_start_orig = genome_len - feature.location.end
                        rc_end_orig = genome_len - feature.location.start

                        # 2. Flip the strand
                        new_strand = -feature.location.strand if feature.location.strand else None

                        # 3. Apply the rotation to these new reverse-complemented coordinates
                        # The rotation_point_on_rc is the start of dnaA on the *newly reverse-complemented sequence*.
                        rotated_rc_start = (rc_start_orig - rotation_point_on_rc) % genome_len
                        rotated_rc_end = (rc_end_orig - rotation_point_on_rc) % genome_len
                        
                        new_location = feature.location.__class__(
                            rotated_rc_start, rotated_rc_end, strand=new_strand,
                            circular=feature.location.circular if hasattr(feature.location, 'circular') else False
                        )
                        new_features.append(feature.__class__(
                            location=new_location, type=feature.type, qualifiers=feature.qualifiers
                        ))
                    except Exception as fe:
                        print(f"Warning: Could not reorient and reverse-complement feature '{feature.type}' at {feature.location}. Error: {fe}", file=sys.stderr)
                modified_record.features = new_features

    else:
        print("No suitable BLAST hit found for dnaA sequence. Cannot reorient assembly.", file=sys.stderr)
        return None

    if modified_record:
        # Write the modified assembly to the output file
        # <<< IMPORTANT: Change "fasta" to "genbank" if you want to output GenBank with features
        SeqIO.write(modified_record, output_file, "fasta") 
        print(f"Modified assembly saved to: {output_file}")
    
    return modified_record

# --- Main execution block for command-line interface ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Reorient a circular DNA assembly based on a dnaA gene sequence using blastn (query vs subject).",
        formatter_class=argparse.RawTextHelpFormatter # For better formatting of help message
    )
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help="Path to the input assembly file (e.g., genome.fasta or genome.gb).\n"
             "Assumes a single, circular contig. Format can be FASTA or GenBank.\n"
             "This file will be used as the -subject for blastn."
    )
    parser.add_argument(
        '-g', '--gene',
        type=str,
        required=True,
        help="Path to a FASTA file containing the dnaA sequence to search for.\n"
             "Expected to be a single sequence record. This will be the -query for blastn."
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help="Path for the output modified assembly file (e.g., reoriented_genome.fasta or reoriented_genome.gb).\n"
             "Output format will be FASTA unless specified otherwise in the script."
    )
    parser.add_argument(
        '-m', '--max_mismatches',
        type=int,
        required=True,
        help="The maximum number of allowed mismatches."
    )

    args = parser.parse_args()

   
    
    # --- Execute the reorientation function with the parsed arguments ---
    modified_record = find_and_reorient_dnaA(
        args.input,
        args.gene,
        args.output)
    

if not modified_record:
    sys.exit(1) # Exit with an error code if the operation failed
else:
    print("Script execution completed successfully.")
    

   