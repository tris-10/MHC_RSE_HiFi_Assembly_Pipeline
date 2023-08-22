#!/usr/bin/env python

"""
Target reads are extracted from the pool of split support reads.  Support reads fully overlap target reads to
maximize chimera detection.  Target reads are restricted by length and limited to a maximum number per
unedited PacBio read.
"""

from Bio import SeqIO
import sys
import argparse
import paf_io


class EditedRead:
    """
    Container for split read information
    """
    def __init__(self, read):
        self.read_data = read
        self.read_length = len(read)


class RawRead:
    """
    Container for unedited read information
    """
    def __init__(self, name):
        self.read_name = name
        self.edited_reads = []

    def add_edited_read(self, read_data):
        split_read = EditedRead(read_data)
        self.edited_reads.append(split_read)

    def get_top(self, read_cap, min_length):
        """
        Return read_cap split reads that are above min_length.

        :param read_cap: Maximum number of reads to return
        :param min_length: Minimum read length to return
        :return: list of EditedRead retained for RawRead
        """

        sort_reads = sorted(self.edited_reads, key=lambda x: x.read_length, reverse=True)
        return [x.read_data for x in sort_reads[0: read_cap] if x.read_length >= min_length]

    def get_count(self):
        return len(self.edited_reads)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('edited_support_reads', help='Path to edited support reads in fasta format')
    parser.add_argument('target_reads', help='Path to target reads in fasta format')
    parser.add_argument('edited_target_reads', help='Path to edited target read output')
    parser.add_argument('-c', '--max_count', help='Max number of edited reads per raw read', default=4, type=int)
    parser.add_argument('-l', '--min_length', help='Minimum edited read length', default=1000, type=int)
    args = parser.parse_args()

    target_reads = load_target_reads(args.target_reads)
    load_split_reads(args.edited_support_reads, target_reads)
    write_target_reads(args.edited_target_reads, target_reads, args.max_count, args.min_length)


def load_target_reads(target_read_file):
    """
    Load target read names and return a dictionary of RawRead objects keyed by the name.

    :param target_read_file: Path to unedited target reads in fasta format
    :return: dictionary of RawRead objects keyed by the name
    """

    target_reads = {}
    with open(target_read_file, 'r') as ip_fasta:
        for read in SeqIO.parse(ip_fasta, 'fasta'):
            read_name = paf_io.get_original_name(read.id)
            target_reads[read_name] = RawRead(read_name)
    print('Loaded {0} raw target reads'.format(len(target_reads.keys())))
    return target_reads


def load_split_reads(split_read_file, target_reads):
    """
    Load split reads and add to raw read container.

    :param split_read_file: Path to edited support read file in fasta format
    :param target_reads: dictionary of RawRead objects keyed by the name
    :return: number of split target reads
    """

    with open(split_read_file, 'r') as ip_fasta:
        for read in SeqIO.parse(ip_fasta, 'fasta'):
            read_name = paf_io.get_original_name(read.id)
            if read_name not in target_reads:
                continue

            target_reads[read_name].add_edited_read(read)


def write_target_reads(target_read_output, target_reads, max_count, min_length):
    """
    Write out edited target reads

    :param target_read_output: Path to edited target read file in fasta file
    :param target_reads: dictionary of target read objects
    :param max_count: Max number of edited reads per raw read
    :param min_length: Minimum edited read length
    """

    no_split = one_split = three_split = more_split = new_count = 0

    with open(target_read_output, 'w') as op_fasta:
        for raw_read in target_reads.values():
            count = raw_read.get_count()
            if count == 1:
                no_split += 1
            elif count == 2:
                one_split += 1
            elif count <= 4:
                three_split += 1
            else:
                more_split += 1
            for read in raw_read.get_top(max_count, min_length):
                new_count += 1
                SeqIO.write(read, op_fasta, 'fasta')

    print('Unedited target reads', no_split)
    print('1 split', one_split)
    print('2-3 split', three_split)
    print('>3 split', more_split)
    print('Edited target reads',new_count)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit(1)
