#!/usr/bin/env python

"""
Identify reads that align to the edges of contigs.  These reads will not be cleaned properly by alignment, since the
alignment will terminate at the edge of the contig vs the edge of a chimera.  All reads within a certain distance
of the edge are written to file, along with 'support' reads, which are the set of reads that overlap the edge reads.
The support reads will be used to identify breakpoints and only edge reads will be retained for assembly.
The original cleaned reads are not used because the breakpoint finding was done with diploid data, which makes
breakpoint finding more difficult.

During chimera detection, edges are treated separately. The only exception to this is when a read aligns to two
different edges.  In this case, the edges are considered grouped and processed together.  This can happen if there
is only a small gap between two contigs in the assembly.
"""

from Bio import SeqIO
from collections import defaultdict
import sys
import argparse
import pysam
import copy
import fileinput
import paf_io


class ContigEdge:
    """
    Container for contig edge-associated information.  Stores the names and sequences of both the edge and supporting
    reads associated with the contig edge. Also stores a list of other edges that are linked due to shared reads.
    """
    def __init__(self, name):
        self.name = name
        self.grouped_name = None
        self.linked_edges = set()
        self.edge_read_names = set()
        self.support_read_names = set()
        self.all_reads = {}

    def add_support_read_name(self, read_name):
        self.support_read_names.add(read_name)

    def add_edge_read_name(self, read_name):
        self.edge_read_names.add(read_name)

    def add_read(self, read):
        read = copy.copy(read)
        read_name = read.id
        read.id = '{0}/{1}'.format(read.id, self.grouped_name)
        read.description = ''
        self.all_reads[read_name] = read

    def add_linked_edge(self, linked_edge):
        self.linked_edges.add(linked_edge)


def main():
    # Parse command line argumnets
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_align', help='Alignment of unedited and haplotype-specific PacBio reads to contigs in '
                                            'BAM format')
    parser.add_argument('input_fasta', help='Unedited PacBio reads')
    parser.add_argument('output_edge_fasta', help='Reads that align to the contig edges in fasta format')
    parser.add_argument('output_support_fasta', help='Reads that align to the contig edges + padding in fasta format')
    parser.add_argument('-e', '--edge_dist', help='Max distance from contig end to retain', type=int, default=500)
    parser.add_argument('-s', '--support_dist', help='Max distance from edge_dist to retain', type=int, default=4000)
    parser.add_argument('-l', '--min_length', help='Minimum contig length', type=int, default=5000)
    parser.add_argument('-a', '--min_align_score', help='Minimum alignment score', type=int, default=0)
    args = parser.parse_args()

    print('>Edge read identification starting')
    print('>Loading edge reads')
    edge_dict = parse_contig_edge_data(args.input_align, args.edge_dist, args.support_dist, args.min_length,
                                       args.min_align_score)

    print('>Grouping shared edges')
    group_edges(edge_dict)

    # TODO: This step isn't required if minimap2 was used for alignment with -Y option. In that case full seq in bam file
    print('>Loading sequence information')
    load_read_data(args.input_fasta, edge_dict)

    print('>Writing edge and support reads')
    write_edge_reads(args.output_edge_fasta, args.output_support_fasta, edge_dict)

    print('>Edge read identification finished')


def parse_contig_edge_data(align_file, edge_dist, support_dist, min_length, min_align_score, custom_edge_dist_dict={}):
    """
    Read through an alignment file and collect the reads names of edge reads and supporting reads.  If the same
    read aligns to two different edges, track the connected edges.

    :param align_file: Alignment of unedited and haplotype-specific PacBio reads to contigs
    :param edge_dist: Reads that align within edge_dist from the end of the contigs are considered edge reads.
    :param support_dist: Reads that align within edge_dist + support dist from the end of the contigs are support reads.
    :param min_length: Minimum contig length to process
    :param min_align_score: Minimum alignment score
    :param custom_edge_dist_dict: Edge distance can be overridden using values in this dictionary
    :return: Dictionary of ContigEdge objects, keyed by their names [ contig_name_start|end_contig_length ]
    """
    edge_dict = {}

    with pysam.AlignmentFile(align_file, "rb") as ip_bam:
        for contig_name, contig_length in zip(ip_bam.references, ip_bam.lengths):
            if contig_length < min_length or not contig_name.startswith('tig'):
                continue

            # load contig start information
            edge_name = '{0}_start_{1}'.format(contig_name, contig_length)
            edge_dict[edge_name] = ContigEdge(edge_name)
            if edge_name in custom_edge_dist_dict:
                load_edge_reads(ip_bam, contig_name, 0, custom_edge_dist_dict[edge_name], edge_dict, edge_name,
                                min_align_score)
                load_support_reads(ip_bam, contig_name, 0, custom_edge_dist_dict[edge_name] + support_dist, edge_dict,
                                   edge_name, min_align_score)
            else:
                load_edge_reads(ip_bam, contig_name, 0, edge_dist, edge_dict, edge_name, min_align_score)
                load_support_reads(ip_bam, contig_name, 0, edge_dist + support_dist, edge_dict, edge_name,
                                   min_align_score)

            # load contig end information
            edge_name = '{0}_end_{1}'.format(contig_name, contig_length)
            edge_dict[edge_name] = ContigEdge(edge_name)
            if edge_name in custom_edge_dist_dict:
                load_edge_reads(ip_bam, contig_name, contig_length - custom_edge_dist_dict[edge_name], contig_length,
                                edge_dict, edge_name, min_align_score)
                load_support_reads(ip_bam, contig_name, contig_length - (custom_edge_dist_dict[edge_name]+support_dist),
                                   contig_length, edge_dict, edge_name, min_align_score)
            else:
                load_edge_reads(ip_bam, contig_name, contig_length - edge_dist, contig_length, edge_dict, edge_name,
                                min_align_score)
                load_support_reads(ip_bam, contig_name, contig_length - (edge_dist + support_dist), contig_length,
                                   edge_dict, edge_name, min_align_score)
    return edge_dict


def load_edge_reads(ip_bam, contig_name, start, end, edge_dict, edge_name, min_align_score):
    """
    Identify reads that align to contig edges.  If a read is already associated with an edge, link the two edge
    objects.

    :param ip_bam: handle to open alignment file
    :param contig_name: name of the contig
    :param start: starting coordinate of contig edge region
    :param end: ending coordinate of contig edge region
    :param edge_dict: Dictionary of ContigEdge objects, keyed by their names
    :param edge_name: Name of contig edge [ contig_name_start|end_contig_length ]
    :param min_align_score: Minimum alignment score
    """

    for align in ip_bam.fetch(contig_name, start, end):
        if align.get_tag('AS') < min_align_score:
            continue
        read_name = paf_io.get_original_name(align.query_name)

        # If the edge read aligns to a different edge, note the linkage so they can be grouped.
        for current_edge in edge_dict.values():
            if read_name in current_edge.edge_read_names and edge_name != current_edge.name:
                edge_dict[current_edge.name].add_linked_edge(edge_name)
                edge_dict[edge_name].add_linked_edge(current_edge.name)
        edge_dict[edge_name].add_edge_read_name(read_name)


def load_support_reads(ip_bam, ref_name, start, end, edge_dict, edge_name, min_align_score):
    """
    Identify reads supporting edge reads, which is defined as the region edge_dist + support_dist from contig edge

    :param ip_bam: handle to open alignment file
    :param ref_name: name of the contig
    :param start: starting coordinate of the contig support region
    :param end: ending coordinate of the contig support region
    :param edge_dict: Dictionary of ContigEdge objects, keyed by their names
    :param edge_name: Name of contig edge [ contig_name_start|end_contig_length ]
    :param min_align_score: Minimum alignment score
    """
    for align in ip_bam.fetch(ref_name, start, end):
        if align.get_tag('AS') < min_align_score:
            continue
        edge_dict[edge_name].add_support_read_name(paf_io.get_original_name(align.query_name))


def group_edges(edge_dict):
    """
    Group edges together and label with a shared name

    :param edge_dict: Dictionary of ContigEdge objects, keyed by their names
    """

    e_idx = 0
    for current_edge in edge_dict.values():
        # If the current edge shares reads with another edge, label them all with the same name
        if len(current_edge.linked_edges) > 0:
            if current_edge.grouped_name is None:
                current_edge.grouped_name = 'JG{0}'.format(e_idx)
                print('Creating edge group {0} {1}'.format(current_edge.name, current_edge.grouped_name))
                e_idx += 1

            # Label each associated edge
            for linked_edge_name in current_edge.linked_edges:
                linked_edge = edge_dict[linked_edge_name]
                if linked_edge.grouped_name is None:
                    linked_edge.grouped_name = current_edge.grouped_name
                    print('Adding to edge group {0} {1}'.format(linked_edge.name, linked_edge.grouped_name))

        # Current edge is not linked to anything else
        else:
            current_edge.grouped_name = 'JG{0}'.format(e_idx)
            print('Solo edge {0} {1}'.format(current_edge.name, current_edge.grouped_name))
            e_idx += 1


def load_read_data(input_fasta, edge_dict):
    """
    Load fasta data and add to edge objects

    :param input_fasta: Unedited PacBio reads
    :param edge_dict: Dictionary of ContigEdge objects, keyed by their names
    """
    with fileinput.hook_compressed(input_fasta, 'rt') as ip_fasta:
        for read in SeqIO.parse(ip_fasta, 'fasta'):
            read_name = paf_io.get_original_name(read.id)
            for edge in edge_dict.values():
                if read_name in edge.support_read_names:
                    edge.add_read(read)


def write_edge_reads(output_edge, output_support, edge_dict):
    """
    Write out edge and support reads to file.

    :param output_edge: Path to edge read output file
    :param output_support: Path to support read output file
    :param edge_dict: Dictionary of ContigEdge objects, keyed by their names
    :return:
    """
    with open(output_edge, 'w') as of_edge, open(output_support, 'w') as of_support:
        written_edges = defaultdict(set)
        written_support = defaultdict(set)

        for edge in edge_dict.values():
            print('Contig: {0} has {1} edge reads and {2} supporting'.format(edge.name, len(edge.edge_read_names),
                                                                             len(edge.support_read_names)))
            for read_name in edge.edge_read_names:
                if read_name not in written_edges[edge.grouped_name]:
                    SeqIO.write(edge.all_reads[read_name], of_edge, 'fasta')
                    written_edges[edge.grouped_name].add(read_name)

            for read_name in edge.support_read_names:
                if read_name not in written_support[edge.grouped_name]:
                    SeqIO.write(edge.all_reads[read_name], of_support, 'fasta')
                    written_support[edge.grouped_name].add(read_name)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit(1)


