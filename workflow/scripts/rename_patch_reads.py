#!/usr/bin/env python

from Bio import SeqIO
import sys
import argparse
import re

def main():
    parser = argparse.ArgumentParser(description="Renames final read")
    parser.add_argument("split_patch_in")
    parser.add_argument("split_patch_out")
    args = parser.parse_args()

    with open(args.split_patch_in, 'r') as ipf, open(args.split_patch_out, 'w') as opf:
        for record in SeqIO.parse(ipf, 'fasta'):
            seq_name = re.sub("ccs\d+","ccs", record.id)
            record.id = seq_name
            record.description = ""
            SeqIO.write(record, opf, 'fasta')


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit()
