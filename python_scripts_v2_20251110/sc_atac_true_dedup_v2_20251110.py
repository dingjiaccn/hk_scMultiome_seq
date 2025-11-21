
#!/usr/bin/env python3
from __future__ import print_function
import pysam
import sys

inbam = sys.argv[1]
outbam = sys.argv[2]

readsin = pysam.AlignmentFile(inbam, "rb")
readsout = pysam.AlignmentFile(outbam, "wb", template=readsin)

last_reference = None
observed_reads = set()

for read in readsin:
    # skip unmapped or secondary/supplementary reads
    if read.is_unmapped or read.mate_is_unmapped or read.is_secondary or read.is_supplementary:
        continue

    reference_id = read.reference_id

    # reset per-chromosome set to save memory
    if last_reference != reference_id:
        print(f"Deduplicating {read.reference_name}...")
        last_reference = reference_id
        observed_reads = set()

    # --- per-cell key (CB tag required) ---
    if not read.has_tag("CB"):
        continue
    cb = read.get_tag("CB")

    # --- define fragment coordinates ---
    # tlen can be negative depending on read orientation
    if read.tlen < 0:
        fragstart = read.mpos - read.tlen   # rightmost read starts after leftmost
        fragend = read.mpos
    else:
        fragstart = read.pos
        fragend = read.pos + read.tlen

    # Build a unique fragment key per cell barcode
    read_id = (cb, fragstart, fragend, read.is_read1)

    if read_id not in observed_reads:
        observed_reads.add(read_id)
        readsout.write(read)

readsin.close()
readsout.close()
print("Per-cell deduplication completed.")
