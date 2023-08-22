#!/usr/bin/env python

"""
Identify patch reads
"""

from Bio import SeqIO
from collections import defaultdict
import sys
import argparse
import identify_edge_reads as ier
import paf_io

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_alignment', help='Alignment of unedited PacBio reads to the haploid contigs, bam format')
    parser.add_argument('self_alignment', help='Path to self-alignment alignment file, paf format')
    parser.add_argument('matching_hap_fasta',
                        help='Fasta file containing the reads used to generate target haplotype contigs')
    parser.add_argument('opposite_hap_fasta',
                        help='Fasta file containing the reads used to generate the opposite haplotype contigs')
    parser.add_argument('input_fasta', help='Unedited PacBio read fasta')
    parser.add_argument('output_edge_fasta', help='Reads that align to the contig edges in fasta format')
    parser.add_argument('output_support_fasta', help='Reads that align to the contig edges + padding in fasta format')
    parser.add_argument('-e', '--edge_dist', help='Max distance from contig end to retain', type=int, default=500)
    parser.add_argument('-s', '--support_dist', help='Max distance from edge_dist to retain', type=int, default=2000)
    parser.add_argument('-l', '--min_length', help='Minimum contig length', type=int, default=5000)
    parser.add_argument('-a', '--min_align_score', help='Minimum alignment score', type=int, default=500)
    parser.add_argument('-m', '--overlap_distance', help='Max self-align distance from contig edge', type=int, default=2000)
    args = parser.parse_args()

    print('>Starting patch read identification')
    print('>Loading contig overlaps')
    contig_overlaps = identify_contig_overlaps(args.self_alignment, args.min_length, args.overlap_distance)

    print('>Loading edge reads')
    edge_dict = ier.parse_contig_edge_data(args.input_alignment, args.edge_dist, args.support_dist,
                                                           args.min_length, args.min_align_score, contig_overlaps)

    for edge in edge_dict.values():
        print(edge.name, len(edge.edge_read_names), len(edge.support_read_names))

    print('>Grouping shared edges')
    ier.group_edges(edge_dict)

    # TODO: This step isn't required if minimap2 was used for alignment with -Y option. In that case full seq in bam file
    print('>Loading sequence information')
    ier.load_read_data(args.input_fasta, edge_dict)
    match_reads = load_read_names(args.matching_hap_fasta)
    opp_reads = load_read_names(args.opposite_hap_fasta)
    cross_assign_reads = identify_cross_assignment(match_reads, opp_reads)

    print('>Remove non-haplotype reads')
    filter_non_haplotpe_reads(edge_dict, match_reads, cross_assign_reads)

    print('>Writing edge and support reads')
    ier.write_edge_reads(args.output_edge_fasta, args.output_support_fasta, edge_dict)

    print('>Patch read identification finished')


def filter_non_haplotpe_reads(edge_dict, match_reads, cross_reads):
    for edge in edge_dict.values():
        updated_edge_names = set()
        updated_support_names = set()
        for read_name in edge.edge_read_names:
            if read_name in match_reads and read_name not in cross_reads:
                updated_edge_names.add(read_name)

        for read_name in edge.support_read_names:
            if read_name in match_reads and read_name not in cross_reads:
                updated_support_names.add(read_name)

        edge.support_read_names = updated_support_names
        edge.edge_read_names = updated_edge_names


def load_read_names(input_fasta):
    read_name_dict = defaultdict(set)
    with open(input_fasta, 'r') as ipf:
        for record in SeqIO.parse(ipf, 'fasta'):
            read_name_dict[paf_io.get_original_name(record.id)].add(record.id)
    return read_name_dict


def identify_cross_assignment(hap1_reads, hap2_reads):
    cross_assign_reads = set()
    for read in hap1_reads.keys():
        if read in hap2_reads:
            if len(hap1_reads[read]) < len(hap2_reads[read]) or hap1_reads[read].isdisjoint(hap2_reads[read]):
                cross_assign_reads.add(read)
    return cross_assign_reads


def identify_contig_overlaps(self_align_paf, min_contig_length, max_edge_distance):
    """

    :param self_align_paf:
    :param min_contig_length:
    :param max_edge_distance:
    :return:
    """
    custom_edge_dist = {}
    for align_group in paf_io.load_query_alignments(self_align_paf):
        start_overlap = None
        end_overlap = None
        for align in align_group:
            if align.query_length >= min_contig_length:
                if align.query_start < max_edge_distance and (start_overlap is None or start_overlap < align.query_end):
                    start_overlap = align.query_end
                if ((align.query_length - align.query_end) < max_edge_distance
                        and (end_overlap is None or end_overlap < (align.query_length - align.query_start))):
                    end_overlap = (align.query_length - align.query_start)
        if start_overlap is not None:
            custom_edge_dist["{0}_{1}_{2}".format(align_group[0].query_name, "start", align_group[0].query_length)] = start_overlap
        if end_overlap is not None:
            custom_edge_dist["{0}_{1}_{2}".format(align_group[0].query_name, "end", align_group[0].query_length)] = end_overlap
    return custom_edge_dist

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit(1)

