#!/usr/bin/env python
import argparse
import gzip
from Bio import SeqIO

def open_fastq(filename):
    """Open .fastq or .fastq.gz files correctly"""
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")  # text mode for FASTQ
    else:
        return open(filename, "r")

def filter_fastq(input_file, output_file, keep_ids):
    """Write reads to output if their ID is in keep_ids"""
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(open_fastq(input_file), "fastq"):
            if record.id.split()[0] in keep_ids:
                SeqIO.write(record, out_handle, "fastq")

def main():
    parser = argparse.ArgumentParser(description="Keep only reads with common IDs between R1 and R2, and apply to R1, R2, I1, and I2.")
    parser.add_argument("-1", "--r1", required=True, help="R1 FASTQ or FASTQ.GZ file")
    parser.add_argument("-2", "--r2", required=True, help="R2 FASTQ or FASTQ.GZ file")
    parser.add_argument("-I1", "--index1", required=True, help="Index1 FASTQ or FASTQ.GZ file")
    parser.add_argument("-I2", "--index2", required=True, help="Index2 FASTQ or FASTQ.GZ file")
    parser.add_argument("-o", "--output_prefix", required=True, help="Prefix for output FASTQ files")

    args = parser.parse_args()

    print("Reading R1 and R2 to find common read IDs...")
    r1_ids = {record.id.split()[0] for record in SeqIO.parse(open_fastq(args.r1), "fastq")}
    r2_ids = {record.id.split()[0] for record in SeqIO.parse(open_fastq(args.r2), "fastq")}
    common_ids = r1_ids & r2_ids

    print(f"Found {len(common_ids):,} common read IDs.")

    print("Filtering all four FASTQ files...")
    filter_fastq(args.r1, f"{args.output_prefix}_R1.fastq", common_ids)
    filter_fastq(args.r2, f"{args.output_prefix}_R2.fastq", common_ids)
    filter_fastq(args.index1, f"{args.output_prefix}_I1.fastq", common_ids)
    filter_fastq(args.index2, f"{args.output_prefix}_I2.fastq", common_ids)

    print("Filtering complete.")

if __name__ == "__main__":
    main()
