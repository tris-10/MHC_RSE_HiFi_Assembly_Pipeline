#!/usr/bin/env python

"""
Write out unedited PacBio reads that were the source of a list of haplotype-specific edited reads.
"""

from Bio import SeqIO
import sys
import argparse
import paf_io


def main():

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('unedited_fasta_in', help='Path to unedited PacBio reads in fasta format')
    parser.add_argument('hap_fasta_in', help='Path to haplotype-specific, edited PacBio reads in fasta format')
    parser.add_argument('hap_fasta_out', help='Path to haplotype-specific, unedited PacBio reads output in fasta format')
    args = parser.parse_args()

    read_name_set = set()

    with open(args.hap_fasta_in, 'r') as ip_fasta:
        for read in SeqIO.parse(ip_fasta, 'fasta'):
            read_name_set.add(paf_io.get_original_name(read.id))

    with open(args.unedited_fasta_in, 'r') as ip_fasta, open(args.hap_fasta_out, 'w') as op_fasta:
        for read in SeqIO.parse(ip_fasta, 'fasta'):
            if paf_io.get_original_name(read.id) in read_name_set:
                read.description = ''
                SeqIO.write(read, op_fasta, "fasta")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit(1)

