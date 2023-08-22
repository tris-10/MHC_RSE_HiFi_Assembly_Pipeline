#!/usr/bin/env python

"""
Identify and remove any chimeras at the end of contig sequences
"""

from Bio import SeqIO
import sys
import argparse
import re
import subprocess
import paf_io as pio
import contig_io as cio
import os_utils as osu

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_contig_path', help='Path to input contig fasta file')
    parser.add_argument('trim_contig_path', help='Path to output contig fasta file')
    parser.add_argument('minimap2_path', help='Path to minimap2 binary')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Minimap2 threads')
    parser.add_argument('-m', '--min_contig_length', type=int, default=20000,
                        help='Minimum contig length, shorter contigs are not written to output.')
    parser.add_argument('-i', '--min_identity', type=float, default=0.98,
                        help='Minimum alignment identity, alignments below this identity are ignored.')
    parser.add_argument('-o', '--min_overlap', type=int, default=1000,
                        help='Minimum alignment overlap length, alignments below this length are ignored.')
    parser.add_argument('-d', '--max_edge_distance', type=int, default=5000,
                        help='Maximum alignment distance from contig edge, alignments further from the edge are ignored')
    parser.add_argument('-f', '--full_overlap_frac', type=float, default=0.9,
                        help='Amount of overlap to consider an interval completely overlapping')
    args = parser.parse_args()

    print('>Starting Contig Trimming')
    self_align_path = args.trim_contig_path + '_self.paf'

    # Check inputs
    for path in [args.input_contig_path]:
        if not osu.check_path(path):
            print('Could not find data at specified path: {0}'.format(path))
            sys.exit(1)

    print('>Loading Contig information')
    tig_data = cio.create_tiginfo_dict(args.input_contig_path)

    print('>Create self alignments')
    align_contigs_to_self(args.input_contig_path, self_align_path, args.minimap2_path, args.threads)

    print('>Trim contigs by overlap')
    trim_contigs_by_overlap(self_align_path, tig_data, args.full_overlap_frac, args.min_identity, args.min_overlap,
                            args.max_edge_distance)

    print('>Write out trimmed contigs')
    cio.write_trimmed_contigs(args.trim_contig_path, tig_data, args.min_contig_length)

    osu.remove_file(self_align_path)
    print(">Finished Contig Trimming")


def trim_contigs_by_overlap(self_alignment_file, tig_data, full_overlap_frac, min_identity, min_overlap, max_edge_distance):
    for align_group in pio.load_query_alignments(self_alignment_file):
        for a in align_group:
            full_overlap = True if max(0, (min(a.query_end, a.target_end) - max(a.query_start, a.target_start))) / a.query_length > full_overlap_frac else False

            if (a.query_name == a.target_name and a.strand == '-' and a.identityNM >= min_identity and
                    a.align_length > min_overlap):

                if a.query_start < max_edge_distance and a.target_start < max_edge_distance * 2:
                    if not full_overlap:
                        pred_start = min(a.query_start, a.target_start)
                    else:
                        pred_start = a.query_start + int(a.match_count / 2)

                    tig_data[a.query_name].trim_start(pred_start)
                if (a.query_length - a.query_end) < max_edge_distance and (a.target_length - a.target_end) < max_edge_distance * 2:
                    if not full_overlap:

                        pred_end = min(a.query_length - a.query_start, a.target_length - a.target_start)
                        print('not_full', a.query_start, a.target_start, a.query_length-pred_end)
                    else:
                        print('full')
                        pred_end = (a.query_length - a.query_end) + int(a.match_count / 2)
                    tig_data[a.query_name].trim_end(pred_end)


def align_contigs_to_self(contig_file, align_file, minimap_path, threads):
    """
    Align contigs to themselves using minimap2

    :param contig_file: Path to the location/read count filtered contig file
    :param align_file: Path to output alignment file
    :param minimap_path: Path the minimap2 executable
    :param threads: Number of threads to use for alignment
    :return: Path to the alignment file in PAF format.
    """
    cproc = subprocess.run([minimap_path, '-x', 'map-hifi', '-D', '-c', '-t', str(threads), '-o',
                             align_file, contig_file, contig_file])

    if cproc.returncode != 0:
        print('Error running self-alignment')
        sys.exit(1)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("user interrupted, exiting")
        sys.exit(1)
