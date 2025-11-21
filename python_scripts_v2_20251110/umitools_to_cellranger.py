import pandas as pd
import numpy as np
import os
import argparse
from scipy import sparse
from scipy.io import mmwrite

def main(input_file, output_dir):
    # Load umi_tools count data
    df = pd.read_csv(input_file, sep="\t")

    # Do NOT strip GRCh38_ or mm10_ prefixes
    genes = pd.Index(df["gene"].unique(), name="gene")
    barcodes = pd.Index(df["cell"].unique(), name="barcode")

    # Map to indices
    gene_to_row = {gene: i for i, gene in enumerate(genes)}
    barcode_to_col = {bc: i for i, bc in enumerate(barcodes)}

    # Build COO sparse matrix
    row = df["gene"].map(gene_to_row)
    col = df["cell"].map(barcode_to_col)
    data = df["count"]
    matrix = sparse.coo_matrix((data, (row, col)), shape=(len(genes), len(barcodes)))

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Write features.tsv.gz
    features_df = pd.DataFrame({
        "gene_id": genes,
        "gene_name": genes,
        "feature_type": ["Gene Expression"] * len(genes)
    })
    features_df.to_csv(os.path.join(output_dir, "features.tsv.gz"),
                       sep="\t", index=False, header=False, compression="gzip")

    # Write barcodes.tsv.gz
    barcodes.to_series().to_csv(os.path.join(output_dir, "barcodes.tsv.gz"),
                                index=False, header=False, compression="gzip")

    # Write matrix.mtx.gz
    matrix_file = os.path.join(output_dir, "matrix.mtx")
    mmwrite(matrix_file, matrix)
    os.system(f"gzip -f {matrix_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert umi_tools counts.tsv to Cell Ranger-style matrix.")
    parser.add_argument("-i", "--input", required=True, help="Path to umitools_counts.tsv")
    parser.add_argument("-o", "--output", required=True, help="Output directory for Cell Ranger-style files")
    args = parser.parse_args()

    main(args.input, args.output)
