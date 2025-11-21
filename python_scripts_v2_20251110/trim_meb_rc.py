#!/usr/bin/env python
import argparse
import subprocess
import sys
import gzip
import time
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os
import re
import itertools

def generate_mismatches(sequence_list, num_mismatches, allow_n=False, reverse_complement=False):
    """
    Generate a list of mimatched sequences to a given sequence. Must only contain ATGC.
    """
    mismatch_list = []
    letters = 'ACGT'

    if allow_n:
        letters += 'N'
    for sequence in sequence_list:
        if reverse_complement:
            newsequence = sequence.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
            newsequence = newsequence.upper()[::-1]
        else:
            newsequence = sequence.upper()
        
        mismatches = []
        for locs in itertools.combinations(range(len(newsequence)), num_mismatches):
            sequence_chars = [[char] for char in newsequence]
            for loc in locs:
                orig_char = newsequence[loc]
                sequence_chars[loc] = [l for l in letters if l != orig_char]
            for poss in itertools.product(*sequence_chars):
                mismatches.append(''.join(poss))

        mismatch_list.append(mismatches)
    return mismatch_list

# Parse arguments
parser = argparse.ArgumentParser(description='A program to append barcode to fastq for scidrop data.')
parser.add_argument('-1', '--input1', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R1.')
parser.add_argument('-o', '--output_prefix', required=True, help='Prefix to add to output files.')
parser.add_argument('-s', '--seq', required=True, help='Fixed starting sequence to trim (e.g., MeB adapter).')
args = parser.parse_args()

# Set up input and output files
output_file_2 = f'{args.output_prefix}_R2.fastq'
output_stat = f'{args.output_prefix}_R2_stat.txt'
output_r2 = open(output_file_2, 'w')

if2 = FastqGeneralIterator(args.input1)

start = time.time()

meb_seq = args.seq
meb_mismatch_list = generate_mismatches([meb_seq], 2, False, False)[0]
meb_mismatch_list.append(meb_seq[:len(meb_seq)])  # include the original sequence itself

total_reads = 0
valid_reads = 0

for (r2_name, r2_seq, r2_qual) in if2:
    total_reads += 1
    adaptor = r2_seq[0:len(meb_seq)]
    r2_seq_trimmed = r2_seq[len(meb_seq):]
    r2_qual_trimmed = r2_qual[len(meb_seq):]
    if adaptor in meb_mismatch_list:
        valid_reads += 1
        output_r2.write(f'@{r2_name}\n{r2_seq_trimmed}\n+\n{r2_qual_trimmed}\n')

output_r2.close()

# Write stats to file
with open(output_stat, 'w') as statfile:
    statfile.write(f'Total reads:\t{total_reads}\n')
    statfile.write(f'Valid reads:\t{valid_reads}\n')

print('Done trimming in %.2f minutes. R2_stat.txt written.' % ((time.time() - start) / 60.0))
