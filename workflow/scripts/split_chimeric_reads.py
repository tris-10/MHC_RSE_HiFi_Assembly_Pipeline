#!/usr/bin/env python

"""
The MDA reaction introduces chimeras into the PacBio read data (Warris, 2018) (Lasken, 2007).  Most of the chimeras are
palindromic, meaning the second segment of the chimera is an inverted repeat of the first.  Non-palindromic chimeras are
also observed, where the segments are nearby each other in the genome, but either non-adjacent, out of order or inverse
with respect to the genome.  If chimeras are not removed prior to assembly, they could both be incorporated into
individual contig sequences and can cause confusion in the assembly, resulting in highly fragmented contigs.
Reads were cleaned by identifying potential chimeric boundaries in each PacBio read and then splitting the read at
these breakpoints using a method adapted from SACRA (Kiguchi 2021).

K-mer counts across the reads were generated using jellyfish (v2.3.0 -L 2 -C -m 51)

Overlaps between MHC-specific reads were found using minimap2 (v2.17 -x ava-pb -I 20G -m 200 –dual=yes).
Reads with more than 16 self-overlaps longer than 1KB were removed from the analysis.  All overlaps to a given query
read are checked for early termination, which is defined by starting and/or ending coordinates more than 100bp from
the boundaries of both reads in the overlap.  Query read coordinates that are observed as early termination endpoints
3 or more times are retained as potential breakpoints.  A 50 bp window is centered on each potential breakpoint and
the count-1 of each early termination endpoint within the window is added to the count of the potential breakpoint.
All potential breakpoints within the read are ranked by window count.  Potential breakpoints are removed if the
window overlaps with a higher-ranking breakpoint window.  A breakpoint percentage is calculated by dividing the window
count by the sum of the window count and the number of overlaps that pass through the breakpoint position with more
than 100bp of overhang.

The median kmer depth across each high-quality read was calculated using the kmer counts generated by jellyfish.
If the median kmer depth was less than 20, breakpoints with a breakpoint percentage lower than 20% were removed, if the
kmer depth is higher than 20, the breakpoint percentage cutoff of 10%.  The count of the kmer centered on each
high-quality read breakpoint was compared to the median count, if the ratio was greater than 0.15, the breakpoint was
removed.

Breakpoints for low quality reads were removed if the breakpoint percentage was lower than 20%.

Each raw PacBio read was split using the list of filtered breakpoints and each subread was written to file if longer
than 1KB.
"""

from Bio import SeqIO
from collections import OrderedDict
from collections import defaultdict
from multiprocessing import Pool
import sys
import os
import argparse
import pandas as pd
import numpy as np
import functools
import paf_io


class BreakPoint:
    """This class stores information about potential breakpoints within in PacBio read."""

    def __init__(self, pos, window):
        """
        Initialize the BreakPoint object

        :param pos: Overlap termination position within the query read
        :param window: Window size used to group breakpoints
        """

        self.pos = pos
        self.window_iv = pd.Interval(pos-int(window/2), pos+int(window/2))
        self.pos_count = 0
        self.window_count = 0
        self.pass_through_count = 0
        self.break_per = 0
        self.stronger_break_found = False

    def increment_pos_count(self):
        """
        Increment the count of premature overlap terminiations at the read position.
        """

        self.pos_count += 1

    def generate_window_count(self, break_list):
        """
        Calculate the number of premature overlap terminations across the window centered at the read position.
        All potential read breaks are compared to current position. Window count is increased by position count for
        the exact match, window count is increased by position count - 1 for positions within the window.

        :param break_list: List of all potential breakpoints within the read
        """

        for b in break_list:
            if b.pos == self.pos:
                self.window_count += b.pos_count
            elif b.pos in self.window_iv:
                self.window_count += (b.pos_count-1)

    def calculate_break_percentage(self, align_iv_list):
        """
        Calculate the percentage of premature overlap terminations at the position. The full list of read overlaps
        are compared to the current position, identifying the number of overlaps that pass through the position with
        at least 100bp (default) of overhang.  The number of early terminations in the window centered on the position
        is divided by the number of overlaps passing through the position plus the number of early terminations to
        generate the breakpoint percentage.

        :param align_iv_list: List of all overlaps intervals for the read, reduced by 100bp on each end
        """

        for align_iv in align_iv_list:
            if self.pos in align_iv:
                self.pass_through_count += 1

        self.break_per = self.window_count / (self.pass_through_count+self.window_count) * 100


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('align_file', help='Path to all-to-all overlap file in PAF format')
    parser.add_argument('fasta_file', help='Path to raw PacBio reads in fasta format')
    parser.add_argument('kmer_file', help='Path to kmer count file generated by jellyfish, if not specified, '
                                                  'the kmer analysis is not run.')
    parser.add_argument('output_dir', help='Path to output files')
    parser.add_argument('output_prefix', help='Prefix for output files')
    parser.add_argument('-s', '--max_self', help='Maximum number of self overlaps allowed before the read is '
                                                 'dropped', type=int, default=16)
    parser.add_argument('-l', '--min_length', help='Minimum overlap length retained for analysis',
                        type=int, default=1000)
    parser.add_argument('-d', '--end_distance', help='Minimum distance from the end of the read to be considered as a '
                                                     'breakpoint', type=int, default=100)
    parser.add_argument('-w', '--window_size', help='All potential breakpoints within the window size are '
                                                    'grouped when calculating final depth and percentage', type=int,
                        default=50)
    parser.add_argument('-p', '--min_pos_count', help='Minimum breakpoint depth at a given coordinate that will be '
                                                     'retained', type=int, default=3)
    parser.add_argument('-q', '--min_window_count', help='Minimum breakpoint depth across a window that will be '
                                                        'retained', type=int, default=4)
    parser.add_argument('-t', '--cpu', help='Number of cpus to use for breakpoint processing', default=20, type=int)
    parser.add_argument('-r', '--lq_reads', help='File containing a list of low-quality reads in the dataset, one per line.'
                                                  'All reads in the file are considered low quality. If not specified, '
                                                  'all reads are assumed to be high-quality', default=None)
    parser.add_argument('-e', '--min_hq_per', help='High quality reads are split at breakpoints with this percentage or '
                                               'higher, unless the kmer centered on the breakpoint is commonly observed.',
                                                default=10.0, type=float)
    parser.add_argument('-f', '--min_lq_per', help='Low quality reads are split at breakpoints with this percentage or '
                                               'higher', default=20.0, type=float)
    parser.add_argument('-m', '--min_kmer_depth', help='If the median kmer count across a high quality read above this '
                                                     'threshold, the min_hq_per cutoff is used, otherwise min_lq_per is '
                                                     'used', default=20.0, type=float)
    parser.add_argument('-n', '--min_kmer_ratio', help='If the median kmer count across a high quality read is above min_kmer_cov '
                                                       'and the observation count of the kmer spanning a breakpoint '
                                                       'divided by the median kmer count across the read * 100 is greater to '
                                                       'or equal the kmer_ratio, the read is not split at the breakpoint.  Only '
                                                       'used for high quality reads', type=float, default=15.0)

    args = parser.parse_args()

    print('>Starting chimeric read splitting')
    # Containers
    filtered_reads = OrderedDict()
    bp_dict = {}
    kmer_flank = 0
    lq_reads = set()

    if args.lq_reads is not None:
        print('>Read in low quality read names')
        lq_reads = load_lq_read_names(args.lq_reads)

    print('>Read in kmer counts')
    kmer_dict = load_kmer_counts(args.kmer_file)
    if len(kmer_dict.keys()) > 0:
        kmer_flank = int(len(list(kmer_dict.keys())[0])/2)

    print('>Identifying reads with excessive self-overlaps')
    for query_align_list in paf_io.load_query_alignments(args.align_file):
        if len([x for x in query_align_list if x.query_name == x.target_name and x.align_length >= args.min_length]) > args.max_self:
            filtered_reads[query_align_list[0].query_name] = None

    print('>Processing valid overlapss')
    with Pool(processes=args.cpu) as pool:
        for query_name, bp_list in pool.imap_unordered(functools.partial(identify_breakpoints, args.end_distance, args.min_pos_count,
                                                       args.min_window_count, args.window_size), (align_info for align_info
                                                       in paf_io.load_query_alignments(args.align_file, length=args.min_length,
                                                                                       filtered_read_dict=filtered_reads))):
            bp_dict[query_name] = bp_list

    print('>Writing out breakpoint data')
    with open(os.path.join(args.output_dir, args.output_prefix) + "_bp.txt", "w") as of_bp:
        for query_name, bp_list in bp_dict.items():
            for b in sorted(bp_list, key=lambda x: x.pos):
                of_bp.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(query_name, b.pos, b.window_count,
                                                                    b.pass_through_count, b.break_per))

    print('>Filter breakpoints by percentage')
    bp_dict = filter_breakpoints_by_percentage(bp_dict, lq_reads, args.min_hq_per, args.min_lq_per)

    print('>Splitting reads by breakpoint')
    split_reads_by_breakpoint(args.fasta_file, os.path.join(args.output_dir, args.output_prefix), args.min_lq_per, args.min_kmer_depth,
                              args.min_kmer_ratio, args.min_length, bp_dict, kmer_dict, kmer_flank, lq_reads, filtered_reads)

    print('>Writing out filtered read names')
    write_high_self_reads(filtered_reads, os.path.join(args.output_dir, args.output_prefix))

    print('>Finished chimeric read splitting')


def identify_breakpoints(end_distance, min_pos_count, min_window_count, window_size, align_info):
    """
    Identify premature overlaps terminations and store as potential breakpoints.  If there is enough
    evidence for a breakpoint, re-calculate using all breakpoints in the window. Sort potential breakpoints
    by window evidence and return all non-overlapping points.

    :param end_distance: Minimum distance from the end of the read to be considered as a breakpoint
    :param min_pos_count: Minimum breakpoint depth at a given coordinate that will be retained
    :param min_window_count: Minimum breakpoint depth across a window that will be retained
    :param window_size: All potential breakpoints within the window size are grouped when calculating final depth and percentage
    :param align_info: list of PAFAlignment objects for the query sequence
    :return: query name and a list of BreakPoint objects for the given query.
    """
    align_iv_list = []
    break_dict = {}

    for align in align_info:
        # Overlap interval is reduced by end_distance from both ends. If breakpoint intersects with interval,
        # it is considered a pass through read
        align_iv_list.append(pd.Interval(align.query_start+end_distance, align.query_end-end_distance))

        # If the overlap start is more than end_distance from both the target start and query start, it is
        # considered a premature termination and added to the potential breakpoint list
        if align.query_start > end_distance and align.target_start > end_distance:
            if align.query_start not in break_dict:
                break_dict[align.query_start] = BreakPoint(align.query_start, window_size)
            break_dict[align.query_start].increment_pos_count()

        # If the overlap end is more than end_distance from both the target end and query end, it is considered
        # a premature termination and added to the potential breakpoint list
        if align.query_end < (align.query_length-end_distance) and align.target_end < (align.target_length-end_distance):
            if align.query_end not in break_dict:
                break_dict[align.query_end] = BreakPoint(align.query_end, window_size)
            break_dict[align.query_end].increment_pos_count()

    # Only generate window counts if position has at least min_pos_count observations
    for b in break_dict.values():
        if b.pos_count >= min_pos_count:
            b.generate_window_count(break_dict.values())

    retained_breaks = []
    # Sort breaks by decreasing window_count order, if break has not been masked, report and mask all other potential
    # breakpoints in the window
    for b in sorted(break_dict.values(), key=lambda x: x.window_count, reverse=True):
        # TODO: Breakpoints are prioritized based on raw count, try prioritizing based on break percentage
        if not b.stronger_break_found and b.window_count >= min_window_count:
            b.calculate_break_percentage(align_iv_list)
            retained_breaks.append(b)
            for b2 in break_dict.values():
                # TODO: Breakpoints are removed if windows overlap, try removing if point itself within window
                if not b2.stronger_break_found and b.window_iv.overlaps(b2.window_iv):
                    b2.stronger_break_found = True

    return align_info[0].query_name, sorted(retained_breaks, key=lambda x: x.pos)


def split_reads_by_breakpoint(fasta_file, output_prefix, min_lq_per, min_kmer_depth, min_kmer_ratio, min_length,
                              bp_dict, kmer_dict, kmer_flank, lq_read_names, filtered_reads):
    """
    Split up reads based on the predicted breakpoints.  Low quality reads are split using all breakpoints in the
    dictionary. High quality reads are treated like low quality reads if the median kmer depth <= min_kmer_depth. If the
    median kmer depth > min_kmer_depth, the kmer centered at the predicted breakpoint is counted across all reads.
    If the kmer is infrequent, (kmer_count / median depth < min_kmer_ratio), the breakpoint is retained. If the kmer
    is frequent, the breakpoint is not used.  Once breakpoints are finalized, the method writes out the split reads
    in fasta format if longer than min_length.  The script also writes out all raw PacBio reads that pass filtering
    for use later in the pipeline.

    :param fasta_file: Path to raw PacBio reads
    :param output_prefix: Prefix for outpt files
    :param min_lq_per: Low quality reads are split at breakpoints with this percentage or higher
    :param min_kmer_depth: If the median kmer count across a high quality read above this threshold, it is treated
                           like a low quality read
    :param min_kmer_ratio: If the ratio of bp kmer count to median count is above this threshold, the kmer is considered
                           frequent and the bp is dropped. Only used for hiqh quality reads
    :param min_length: Minimum read length written to file
    :param bp_dict: List of predicted breakpoints for each read
    :param kmer_dict: dictionary containing the count of everykmer in the dataset
    :param kmer_flank: kmer sizes are odd, flank is kmer size / 2 rounded down the nearest integer
    :param lq_read_names: Set of low quality reads in the dataset
    :param filtered_reads: Ordered dictionary with the list of removed reads as the keys
    """

    with open(fasta_file, 'r') as if_fasta, open(output_prefix + '_filt.fasta', 'w') as of_filt, \
            open(output_prefix + '_split.fasta', 'w') as of_clean,  open(output_prefix + '_split_log.txt', 'w') as of_log:
        for record in SeqIO.parse(if_fasta, 'fasta'):
            orig_name = paf_io.get_original_name(record.id)
            record.description = ''

            # If the sequence has a high level of self overlap, skip, otherwise write it out and continue
            if orig_name in filtered_reads:
                continue
            SeqIO.write(record, of_filt, 'fasta')

            # If there is no recorded breakpoint, write original fasta, otherwise split reads
            if orig_name not in bp_dict:
                write_split_seq(record, 0, len(record), of_clean, min_length)
            else:
                seq = str(record.seq)
                kmer_self_dict = kmer_counts_across_read(seq, kmer_flank)
                med_kmer_depth = calculate_kmer_depth(kmer_self_dict, kmer_dict)

                retained_breaks = []
                for bp in bp_dict[orig_name]:
                    bp_kmer_count, kmer = generate_bp_kmer_count(bp, seq, kmer_flank, kmer_self_dict, kmer_dict)

                    used = False
                    if orig_name in lq_read_names or \
                            (med_kmer_depth > min_kmer_depth and (bp_kmer_count == 0 or bp_kmer_count / med_kmer_depth * 100 < min_kmer_ratio)) or \
                            (med_kmer_depth <= min_kmer_depth and bp.break_per >= min_lq_per):
                        retained_breaks.append(bp)
                        used = True

                    # Write out information to log.
                    of_log.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t\t{6}\t{7}\t{8}\t{9}\n".format(record.id, bp.pos,
                                 bp.break_per, kmer, bp_kmer_count, kmer_self_dict[kmer], med_kmer_depth,
                                 round(bp_kmer_count / med_kmer_depth * 100, 2), used, orig_name in lq_read_names))

                start_pos = 0
                for bp in retained_breaks:
                    write_split_seq(record, start_pos, bp.pos, of_clean, min_length)
                    start_pos = bp.pos
                write_split_seq(record, start_pos, len(record), of_clean, min_length)


def load_kmer_counts(kmer_count_file):
    """
    Load in kmer counts generated by jellyfish. Each line contains the kmer sequence followed by the count of the kmer
    in the dateset. Low count kmers are not included in the output to reduce filesize.

    :param kmer_count_file: Path to a kmer count file generated by jellyfish
    :return: Dictionary of kmer counts, key is the kmer sequence, value is the count of the kmer in the dataset.
    """

    kmer_count_dict = defaultdict(int)
    with open(kmer_count_file, "r") as ipf:
        for line in ipf:
            items = line.strip().split(" ")
            kmer_count_dict[items[0]] = int(items[1])
    return kmer_count_dict


def load_lq_read_names(lq_read_name_file):
    """
    Load all low quality read names, return as a set.

    :param lq_read_name_file: Path to a file with all low quality read names, one per line
    :return: set of low quality read names
    """

    lq_read_names = set()
    with open(lq_read_name_file, "r") as ipf:
        for line in ipf:
            lq_read_names.add(line.strip())
    return lq_read_names


def write_split_seq(record, start, end, opf, min_length):
    """
    Extracts subsequence from original PacBio sequence and writes to file with the sequence boundaries in the header

    :param record: Original PacBio sequence as a SeqRecord object
    :param start: Starting coordinate of the subsequence
    :param end:  Ending coordinate of the subsequence
    :param opf: Handle to the output file
    :param min_length: Minimum read length to write to file
    """

    sub_seq_record = record[start:end]
    if len(sub_seq_record) >= min_length:
        sub_seq_record.id = "{0}:{1}-{2}".format(record.id, start, end)
        sub_seq_record.description = ""
        SeqIO.write(sub_seq_record, opf, "fasta")


def write_high_self_reads(filtered_reads, output_prefix):
    """
    Write out reads with high self-overlaps

    :param filtered_reads: Ordered dictionary with the list of removed reads as the keys
    :param output_prefix: Prefix for output files
    """

    with open(output_prefix + '_high.txt', "w") as of_filt:
        for read_name in filtered_reads.keys():
            of_filt.write(read_name + "\n")


def filter_breakpoints_by_percentage(bp_dict, lq_read_names, min_hq_per, min_lq_per):
    """
    Remove breakpoints where premature overlap ends make up a small portion of the reads aligning to the position.
    The threshold is set lower for high quality reads because the overlap should be more accurate and because kmer
    counts can be used to protect false breakpoints.

    :param bp_dict: List of all potential breakpoints for each read
    :param lq_read_names: Set of low quality reads in the dataset
    :param min_hq_per: Breakpoints in high quality reads with a break_per lower than min_hq_per are removed
    :param min_lq_per: Breakpoints in low quality reads with a break_per lower than min_lq_per are removed
    :return: List of filtered breakpoints for each read.
    """

    # TODO: This was added later, maybe do this filtering as part of the breakpoint detection
    filt_bp_dict = defaultdict(list)
    for read, bp_list in bp_dict.items():
        hq = False if read in lq_read_names else True
        for bp in bp_list:
            if (hq and bp.break_per >= min_hq_per) or (not hq and bp.break_per >= min_lq_per):
                filt_bp_dict[read].append(bp)
    return filt_bp_dict


def generate_kmer(seq, pos, kmer_flank):
    """
    A kmer is generated by retrieving the sequence kmer_flank upstream and downstream of the position.  The
    reverse complement sequence is generated and the two are compared, retaining the lowest lexicographically.

    :param seq: Raw PacBio sequence
    :param pos: Potential breakpoint position
    :param kmer_flank: kmer sizes are odd, flank is kmer size / 2 rounded down the nearest integer
    :return:
    """
    kmer_for = seq[pos-kmer_flank: pos+kmer_flank+1]
    kmer_rev = reverse_complement(kmer_for)
    return kmer_for if kmer_for < kmer_rev else kmer_rev


def kmer_counts_across_read(seq, kmer_flank):
    """
    Generate a dictionary containing the count of every kmer across a read

    :param seq: Raw PacBio sequence
    :param kmer_flank: kmer sizes are odd, flank is kmer size / 2 rounded down the nearest integer
    :return: dictionary containing the count of every kmer in the read
    """
    kmer_self = defaultdict(int)
    for i in range(0, len(seq) - (2*kmer_flank + 1)):
        kmer_self[generate_kmer(seq, i, kmer_flank)] += 1
    return kmer_self


def calculate_kmer_depth(kmer_self, kmer_dict):
    """
    Calculate the median kmer count across the read. If the kmer is not in the dictionary, it is assumed to contain
    a sequencing error and not used for calculating the median.  The kmer count is reduced by the number of times
    it appears in the read itself.

    :param kmer_self: dictionary containing the count of every kmer in the read
    :param kmer_dict: dictionary containing the count of everykmer in the dataset
    :return: median kmer count across the read, or 1 if kmers are never observed.
    """

    kmer_counts = []
    for kmer, self_count in kmer_self.items():
        non_self_cov = kmer_dict[kmer] - self_count
        if non_self_cov > 0:
            kmer_counts.append(non_self_cov)

    # TODO: A median could be calculated on just a handful of positions in the read, maybe require minimum
    if len(kmer_counts) > 0:
        return np.median(kmer_counts)
    else:
        return 1


def generate_bp_kmer_count(bp, seq, kmer_flank, kmer_self, kmer_dict):
    """
    Calculate the number of times the kmer centered on the breakpoint is observed in other reads in the dataset

    :param bp: BreakPoint object
    :param seq: Raw PacBio sequence
    :param kmer_flank: kmer sizes are odd, flank is kmer size / 2 rounded down the nearest integer
    :param kmer_self: dictionary containing the count of every kmer in the read
    :param kmer_dict: dictionary containing the count of every kmer in the dataset.
    :return: kmer observation count in the data and the kmer sequence as a tuple.
    """
    kmer = generate_kmer(seq, bp.pos, kmer_flank)
    kmer_count = 0
    if kmer_dict[kmer] > 0:
        kmer_count = kmer_dict[kmer] - kmer_self[kmer]
    return kmer_count, kmer


def reverse_complement(seq):
    """
    Reverse complement the input sequence.  Assumes no ambiguous bases and uppercase letters.

    :param seq: DNA sequence
    :return: reverse complement sequence
    """

    trans_map = str.maketrans('ACGT', 'TGCA')
    return seq.translate(trans_map)[::-1]


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit(1)

