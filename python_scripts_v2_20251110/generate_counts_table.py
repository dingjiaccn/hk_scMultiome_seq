import pysam                     
from collections import defaultdict
import os                          
from scipy import io               
from scipy import sparse                  
import gzip                        
import sys                      


bam_path = sys.argv[1]       # e.g., RNA_tagged_dedup.featureCounts.bam
output_dir = sys.argv[2]     # e.g., 10x_output

os.makedirs(output_dir, exist_ok=True)

# Read BAM and count CB-XT pairs
bamfile = pysam.AlignmentFile(bam_path, "rb")
cell_gene_counts = defaultdict(lambda: defaultdict(int))
all_genes = set()
all_cells = set()

for read in bamfile.fetch(until_eof=True):
    if read.has_tag("CB") and read.has_tag("XT"):
        cb = read.get_tag("CB")
        gene = read.get_tag("XT")
        cell_gene_counts[cb][gene] += 1
        all_genes.add(gene)
        all_cells.add(cb)

genes_sorted = sorted(all_genes)
cells_sorted = sorted(all_cells)
gene_to_idx = {gene: i for i, gene in enumerate(genes_sorted)}
cell_to_idx = {cell: i for i, cell in enumerate(cells_sorted)}

rows, cols, data = [], [], []
for cell, gene_dict in cell_gene_counts.items():
    for gene, count in gene_dict.items():
        rows.append(gene_to_idx[gene])
        cols.append(cell_to_idx[cell])
        data.append(count)

matrix = sparse.coo_matrix((data, (rows, cols)), shape=(len(genes_sorted), len(cells_sorted)))

# Save uncompressed matrix file
mtx_path = os.path.join(output_dir, "matrix.mtx")
io.mmwrite(mtx_path, matrix)  # Corrected from `counts` to `matrix`

# Compress it to .gz
with open(mtx_path, 'rb') as f_in:
    with gzip.open(mtx_path + ".gz", 'wb') as f_out:
        f_out.writelines(f_in)

# Remove uncompressed file
os.remove(mtx_path)

# Write barcodes.tsv.gz
with gzip.open(os.path.join(output_dir, "barcodes.tsv.gz"), "wt") as f:
    for bc in cells_sorted:
        f.write(f"{bc}\n")

# Write features.tsv.gz
with gzip.open(os.path.join(output_dir, "features.tsv.gz"), "wt") as f:
    for gene in genes_sorted:
        f.write(f"{gene}\t{gene}\tGene Expression\n")
