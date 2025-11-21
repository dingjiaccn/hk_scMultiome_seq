#!/usr/bin/env python3
import pysam
import argparse

def parse_barcode(read_name, assay):
    """
    Extract barcodes from the start of read name:
    e.g. CGTACTAGATTAATTCCAGACACAATTAGATCTACTAAATTGGTTT:NB551502:...
    """
    if assay.upper() == "RNA":
        barcode = read_name.split(":")[0]
        sample_idx_seq = barcode[0:8]
        gem_idx_seq = barcode[8:24]
        tn5_idx_seq = barcode[24:34]
    if assay.upper() == "ATAC":
        barcode = read_name.split(":")[0]
        sample_idx_seq = barcode[0:8]
        gem_idx_seq = barcode[8:24]
        tn5_idx_seq = barcode[24:32]
    return sample_idx_seq, gem_idx_seq, tn5_idx_seq

def get_umi(read_name):
    barcode = read_name.split(":")[0]
    return barcode[-10:]

def tag_bam(input_bam, output_bam, assay):
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)

    for read in bam_in.fetch(until_eof=True):
        if read.is_unmapped:
            continue  # skip unmapped reads

        try:
            sample_idx_seq, gem_idx_seq, tn5_idx_seq = parse_barcode(read.query_name, assay)
        except Exception as e:
            print(f"Failed to parse barcode from {read.query_name}")
            continue

        read.set_tag("SB", sample_idx_seq, value_type="Z")
        read.set_tag("GB", gem_idx_seq, value_type="Z")
        read.set_tag("WB", tn5_idx_seq, value_type="Z")

        tcb_seq = sample_idx_seq+gem_idx_seq+tn5_idx_seq
        read.set_tag("CB", tcb_seq, value_type="Z")

        if assay.upper() == "RNA":
            umi = get_umi(read.query_name)
            read.set_tag("UB", umi, value_type="Z")

        bam_out.write(read)

    bam_in.close()
    bam_out.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add SB (sample), GB (10x gel), WB (well), CB (SB+GB+WB), and UB (RNA) tags to BAM file.")
    parser.add_argument("-i", "--input-bam", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output-bam", required=True, help="Output BAM file with tags added")
    parser.add_argument("-a", "--assay", required=True, choices=["RNA", "ATAC"], help="Assay type: RNA or ATAC")
    args = parser.parse_args()

    tag_bam(args.input_bam, args.output_bam, args.assay)
