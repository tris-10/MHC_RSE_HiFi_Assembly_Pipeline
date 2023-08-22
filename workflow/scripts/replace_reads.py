#!/usr/bin/env python

"""
Write all updated reads followed by all original reads that were not updated to a new file.  Files can be
compressed/uncompressed.
"""

from Bio import SeqIO
from collections import defaultdict
import sys
import argparse
import paf_io
import fileinput


def main():
    parser = argparse.ArgumentParser(description='__doc__')
    parser.add_argument('original_reads', help='path to haplotype-specific, edited and coverage-capped fasta file')
    parser.add_argument('updated_reads', help='Path to cleaned edge read fasta')
    parser.add_argument('output_reads', help='Path to merged fasta file output')
    args = parser.parse_args()

    print('>Starting read replacement')

    updated_read_dict = defaultdict(list)
    orig_count = removed_count = updated_count = 0

    with fileinput.hook_compressed(args.updated_reads, 'rt') as ip_fasta:
        for read in SeqIO.parse(ip_fasta, 'fasta'):
            updated_read_dict[paf_io.get_original_name(read.id)].append(read)
            updated_count += 1

    with fileinput.hook_compressed(args.output_reads, 'wt') as op_fasta, fileinput.hook_compressed(args.original_reads, 'rt') as ip_fasta:
        for read in SeqIO.parse(ip_fasta, 'fasta'):
            orig_count += 1
            orig_name = paf_io.get_original_name(read.id)
            if orig_name not in updated_read_dict:
                SeqIO.write(read, op_fasta, 'fasta')
            elif len(updated_read_dict[orig_name]) > 0:
                SeqIO.write(updated_read_dict[orig_name], op_fasta, 'fasta')
                updated_read_dict[orig_name] = []
                removed_count += 1
            else:
                removed_count += 1
        for reads in updated_read_dict.values():
            if len(reads) > 0:
                SeqIO.write(reads, op_fasta, 'fasta')


    print('Wrote {0} updated reads and {1} original reads for a total of {2}. {3} original '
          'reads were replaced'.format(updated_count, orig_count, updated_count + orig_count, removed_count))
    print('>Finished read replacement')

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit(1)
