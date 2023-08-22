#!/usr/bin/env python

"""
Find overlaps between contigs and drop any that are redundant or have high self alignment. Trim alignment edges
if they appear to be chimeric
"""

import sys
import argparse
import subprocess
import paf_io as pio
import os_utils as osu
import contig_io as cio


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_contig_path', help='Path to canu directory, as specified when running canu.')
    parser.add_argument('filtered_contig_path', help='Path to filtered contig fasta output file.')
    parser.add_argument('minimap2_path', help='Path to minimap2 executable')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of minimap2 threads')
    parser.add_argument('-c', '--min_contig_length', type=int, default=0, help='Minimum contig length')
    g1 = parser.add_argument_group('Overlap Filtering')
    g1.add_argument('-i', '--min_identity', type=float, default=0.96,
                    help='Minimum identity allowed to retain alignment')
    g1.add_argument('-l', '--max_align_length', type=int, default=40000,
                    help='Maximum overlap length allowed to retain alignment')
    g1.add_argument('-o', '--min_self_overlap', type=float, default=0.25,
                    help='Minimum allowed fraction self-overlap')
    g1.add_argument('-p', '--min_other_overlap', type=float, default=0.90,
                    help='Minimum allowed fraction overlap to a single longer contig')
    g1.add_argument('-q', '--min_combined_overlap', type=float, default=0.98,
                    help='Minimum allowed fraction overlap to all longer contigs')
    g2 = parser.add_argument_group('Edge Trimming')
    g2.add_argument('-d', '--max_edge_distance', type=int, default=5000,
                    help='Maximum distance from contig edge to check for terminal chimera')
    g2.add_argument('-f', '--max_edge_overlap', type=int, default=1000,
                    help='Maximum allowed self-overlap bases at contig edge')
    args = parser.parse_args()

    print('>Starting contig filtration by overlap')

    # Check inputs
    for path in [args.input_contig_path]:
        if not osu.check_path(path):
            print('Could not find data at specified path: {0}'.format(path))
            sys.exit(1)

    # Temporary files
    self_align_paf = args.filtered_contig_path + '_self.paf'

    print('>Loading contig information')
    tig_dict = cio.create_tiginfo_dict(args.input_contig_path)

    print('>Finding contig overlaps')
    align_contigs_to_self(args.input_contig_path, self_align_paf, args.minimap2_path, args.threads)

    print('>Filter contigs by overlap')
    filter_by_overlap(self_align_paf, tig_dict, args.min_identity, args.max_align_length, args.min_self_overlap,
                      args.min_other_overlap, args.min_combined_overlap, args.max_edge_overlap, args.max_edge_distance)

    print('>Writing updated contig list')
    cio.write_trimmed_contigs(args.filtered_contig_path, tig_dict, args.min_contig_length)

    print('>Finished contig filtration by overlap')


def align_contigs_to_self(contig_file, align_file, minimap_path, threads):
    """
    Align contigs to themselves using minimap2

    :param contig_file: Path to the location/read count filtered contig file
    :param align_file: Path to output alignment file
    :param minimap_path: Path the minimap2 executable
    :param threads: Number of threads to use for alignment
    :return: Path to the alignment file in PAF format.
    """
    cproc = subprocess.run([minimap_path, '-x', 'map-hifi', '-D', '-t', str(threads), '-o',
                             align_file, contig_file, contig_file])

    if cproc.returncode != 0:
        print('Error running self-alignment')
        sys.exit(1)


def filter_by_overlap(self_align, tig_dict, min_identity, max_align_length, min_self_overlap, min_other_overlap,
                      min_combined_overlap, min_end_overlap, max_end_distance):
    """
    Check contigs for high self-overlap and high overlap with other contigs.  These contigs are set to be filtered
    and not written to output.  Contig with high self-overlap at just the ends are trimmed.

    :param self_align: All-to-all contig alignment file in PAF format
    :param tig_dict: Dictionary of all contig objects reported by Canu
    :param min_identity: Minimum identity allowed to retain alignment
    :param max_align_length: Maximum overlap length allowed to retain alignment
    :param min_self_overlap: Minimum allowed fraction self-overlap
    :param min_other_overlap: Minimum allowed fraction overlap to a single longer contig
    :param min_combined_overlap: Minimum allowed fraction overlap to all longer contigs
    :param min_end_overlap: Maximum allowed self-overlap bases at contig end
    :param max_end_distance: Maximum distance from contig end to check for terminal chimera
    """

    self_paf = pio.load_query_alignments(self_align)

    for align_list in self_paf:
        for a in align_list:
            if a.identity < min_identity:
                continue

            if a.query_name == a.target_name:
                if a.query_overlap >= min_self_overlap:
                    tig_dict[a.query_name].filter_contig()
                    print('Removing {0}, self-aligned with overlap {1} and identity {2}'.format(a.query_name, round(a.query_overlap, 4),
                                                                                                round(a.identity, 4)))
                elif a.match_count > min_end_overlap and a.align_length > min_end_overlap:
                    if a.query_start < max_end_distance and a.target_start < max_end_distance:
                        pred_start = a.query_start + int(a.match_count / 2)
                        tig_dict[a.query_name].trim_start(pred_start)
                    if ((a.query_length - a.query_end) < max_end_distance and
                            (a.target_length - a.target_end) < max_end_distance):
                        pred_end = (a.query_length - a.query_end) + int(a.match_count / 2)
                        tig_dict[a.query_name].trim_end(pred_end)
            elif a.query_length < a.target_length and a.query_length < max_align_length:
                if a.query_overlap >= min_other_overlap:
                    tig_dict[a.query_name].filter_contig()
                    print("Removing {0}, high overlap {1} and identity {2} with {3}".format(a.query_name, round(a.query_overlap, 4),
                                                                                            round(a.identity, 4), a.target_name))
                else:
                    tig_dict[a.query_name].add_contig_overlap(a.query_start, a.query_end)

    for tig in tig_dict.values():
        combined_overlap = tig.check_contig_overlap()
        if combined_overlap >= min_combined_overlap:
            print('Removing {0}, {1} overlap with other contigs'.format(tig.tig_name, round(combined_overlap, 4)))
            tig.filter_contig()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit(1)
