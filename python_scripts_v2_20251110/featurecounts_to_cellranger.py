import pandas as pd
import scipy.io
import gzip
import os
import sys

# Input and output
input_file = sys.argv[1]  # featureCounts txt file
output_dir = sys.argv[2]  # output folder for features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Read in the featureCounts file (skip first two lines)
df = pd.read_csv(input_file, sep='\t', comment='#', skiprows=1)

# First column is Gene ID or Gene Name
genes = df.iloc[:, 0].astype(str)

# Remaining columns are cell barcodes (header = CB)
barcodes = df.columns[6:]

# Count matrix (gene x cell)
counts = df.iloc[:, 6:].astype(int)

# Write matrix.mtx.gz (Matrix Market format)
scipy.io.mmwrite(os.path.join(output_dir, "matrix.mtx"), counts)

# Write features.tsv.gz
with gzip.open(os.path.join(output_dir, "features.tsv.gz"), 'wt') as f:
    for gene in genes:
        f.write(f"{gene}\t{gene}\tGene Expression\n")

# Write barcodes.tsv.gz
with gzip.open(os.path.join(output_dir, "barcodes.tsv.gz"), 'wt') as f:
    for bc in barcodes:
        f.write(f"{bc}\n")

# Compress matrix.mtx
with open(os.path.join(output_dir, "matrix.mtx"), 'rb') as f_in:
    with gzip.open(os.path.join(output_dir, "matrix.mtx.gz"), 'wb') as f_out:
        f_out.writelines(f_in)

# Clean up uncompressed matrix
os.remove(os.path.join(output_dir, "matrix.mtx"))

print("conversion complete. Files written to:", output_dir)
