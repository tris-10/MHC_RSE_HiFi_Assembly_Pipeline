#!/usr/bin/env python

from Bio import SeqIO
import re
import os
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="This script reames split reads so they can be used is a scondary split")
    parser.add_argument("split_patch_in")
    parser.add_argument("split_support_in")
    parser.add_argument("split_patch_out")
    parser.add_argument("split_support_out")
    args = parser.parse_args()

    name_count = {}
    id_dict = {}
    with open(args.split_support_in, 'r') as ipf, open(args.split_support_out, 'w') as opf:
        for record in SeqIO.parse(ipf, 'fasta'):
            name_parts = record.id.split(":")[0].split("/")
            seq_name = "/".join(name_parts[0:3])
            if seq_name not in name_count:
                name_count[seq_name] = 0
            else:
                name_count[seq_name] += 1
            new_id = seq_name + str(name_count[seq_name])
            id_dict[record.id] = new_id
            new_id += "/" + name_parts[3]
            record.id = new_id
            record.description = ""
            SeqIO.write(record, opf, 'fasta')

    with open(args.split_patch_in, 'r') as ipf, open(args.split_patch_out, 'w') as opf:
        for record in SeqIO.parse(ipf, 'fasta'):
            if record.id not in id_dict:
                continue
            record.id = id_dict[record.id]
            record.description = ""
            SeqIO.write(record, opf, 'fasta')

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit()
