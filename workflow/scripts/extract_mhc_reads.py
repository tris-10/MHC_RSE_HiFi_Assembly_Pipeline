#!/usr/bin/env python

"""
The read names from all primary alignments in the input bam file are stored.  The input fastq file is scanned and
matching reads are written to the output file in fasta format.  It is assumed that the input bam contains only
MHC-specific alignments.
"""

from Bio import SeqIO
import pysam
import sys
import argparse
import gzip


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('alignment_file', help='MHC reference panel alignment file in bam format')
    parser.add_argument('fastq_file', help='Gzipped-compressed fastq file containing the pacbio reads used '
                                                   'for alignment.')
    parser.add_argument('output_file', help='MHC-specific reads in fasta format.')
    args = parser.parse_args()

    read_set = load_read_names(args.alignment_file)
    extract_reads(args.fastq_file, args.output_file, read_set)


def extract_reads(input_fastq, output_fasta, read_set):
    """
    Reads through a fastq file and writes out any reads found within read_set

    :param input_fastq: Path to PacBio fastq file
    :param output_fasta: Path to MHC-specific fasta output file
    :param read_set: Set of MHC-specific read names
    """

    print("Writing MHC reads ")
    with gzip.open(input_fastq, "rt") as ipf, open(output_fasta, "w") as opf:
        for record in SeqIO.parse(ipf, "fastq"):
            if record.id in read_set:
                SeqIO.write(record, opf, "fasta")


def load_read_names(input_bam):
    """
    Reads in an alignment file in bam format and stores all non-secondary, non-supplementary reads.

    :param input_bam:
    :return: set of read names.
    """

    read_set = set()
    print("Loading MHC read names")
    with pysam.AlignmentFile(input_bam, "rb") as ipb:
        for record in ipb:
            if not record.is_secondary and not record.is_supplementary:
                read_set.add(record.query_name)
    return read_set


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit(1)
