#!/usr/bin/env python

import pandas as pd
import argparse

def extract_sample(entry):
    """Extract sample ID from tracking file entry."""
    if ":" in entry:
        return entry.split(":")[1].split(".")[0]  # Extract sample before "."
    return "0"  # Convert "-" to "0" for binary matrix

def process_tracking_file(input_file, output_prefix):
    """Process GFFCompare tracking file to generate binary presence table."""
    # Read tracking file
    df = pd.read_csv(input_file, sep='\t', header=None)

    # Define columns based on tracking file format
    df.columns = ['transcript_id', 'locus_id', 'ref_gene', 'class_code'] + [f"sample_{i}" for i in range(1, len(df.columns) - 3)]

    # Extract sample presence as binary values (1 if transcript detected, 0 otherwise)
    for col in df.columns[4:]:  # Sample columns start from index 4
        df[col] = df[col].fillna("-").astype(str).apply(lambda x: "1" if ":" in x else "0")

    # Convert to binary matrix
    df_binary = df[['transcript_id'] + df.columns[4:].tolist()]

    # Save output file
    output_file = f"{output_prefix}_transcript_presence.tsv"
    df_binary.to_csv(output_file, sep='\t', index=False)

    print(f"Processed file saved as: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GFFCompare tracking file to binary presence matrix.")
    parser.add_argument("input_file", help="Path to the GFFCompare tracking file")
    parser.add_argument("output_prefix", help="Prefix for the output file")
    args = parser.parse_args()

    process_tracking_file(args.input_file, args.output_prefix)