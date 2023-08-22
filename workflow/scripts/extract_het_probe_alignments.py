#!/usr/bin/env python

"""
Restrict microarray probe flanking sequence alignments to heterozygous positions.
"""

import argparse
import pysam
import sys

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('het_probe_list', help='Path to heterozygous probe list, first column should be SNP IDs that '
                                               'match the header in the probe fasta files')
    parser.add_argument('output_bam', help='Path to filtered probe alignments')
    parser.add_argument('--input_bam', help='Path to probe alignments in bam format. Default is to '
                                            'read from stdin, use this flag to specify existing file', default='-')
    args = parser.parse_args()

    het_probe_set = load_het_probes(args.het_probe_list)
    filter_bam(args.input_bam, args.output_bam, het_probe_set)


def filter_bam(input_bam, output_bam, het_probe_set):
    """
    Read through probe flanking sequence alignments in bam format and write out only the heterozygous positions

    :param input_bam: Probe alignment file in BAM format
    :param output_bam: Heterozygous probe aligment file output
    :param het_probe_set: Set of heterozygous probe names
    """

    with pysam.AlignmentFile(input_bam, 'rb') as ipb, pysam.AlignmentFile(output_bam, 'wb', template=ipb) as opb:
        for align in ipb:
            if align.query_name in het_probe_set:
                opb.write(align)
    pysam.index(output_bam)

def load_het_probes(probe_file):
    """

    :param probe_file:
    :return:
    """
    het_probe_set = set()
    with open(probe_file, 'r') as ipf:
        for line in ipf:
            items = line.strip().split("\t")
            het_probe_set.add(items[0])
    return het_probe_set


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("USer interrupted, exiting")
        sys.exit(1)