#!/usr/bin/env python

"""
Contigs are patched nby aligning the edges of each contig in the current assembly (r2) to an earlier assembly (r1).
If two r2 contig edges align close to each other on the same r1 contig, a patched contig is created by replacing the
r2 contig edge sequences with the corresponding r1 sequence and adding the remaining r2 contig sequences on either end.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import argparse
import subprocess
import os_utils as osu
import paf_io as pio


class ContigEdgeAlignment:
    """
    Container for contig edge alignment information
    """

    def __init__(self, ce_label, orient, align_start):
        """
        Initialize contig edge alignment object with the name, location and orientation of the r2 contig
        aligment relative to the corresponding r1 contig.

        :param ce_label: Name of the contig edge, which is comprised of the r2 contig name and 'start' or 'end'
        :param orient: Orientation of the r2 contig relative to r1 contig
        :param align_start: Starting coordinate of the alignment
        """
        self.ce_label = ce_label
        self.name, self.contig_end = ce_label.split("_")
        self.orient = orient
        self.align_start = align_start

    def get_seq(self, contig_dict, length):
        """
        Return r2 contig sequence from the start/end of the contig up to the patch location. If the r1 and r2 assemblies
        have differing orientation, reverse complement r2 contig sequence.

        :param contig_dict: Dictionary with r2 contig sequences
        :param length: Amount of sequence from contig end to return
        :return: r2 contig sequence from start/end to patch location
        """
        if self.contig_end == "start":
            sub_seq = contig_dict[self.name][length:]
        else:
            sub_seq = contig_dict[self.name][:-length]
        if self.orient == "-":
            return str(sub_seq.reverse_complement().seq)
        else:
            return str(sub_seq.seq)


class PatchedContig:
    """
    Container for the r2 contig edge alignment locations on a r1 contig.  If two r2 contig edges align nearby on
    the same r1 contig, the r1 sequence across the alignment region are used to join the corresponding r2 contigs.
    """

    def __init__(self, name, start, end, ce_label, orient):
        """
        Initialize with the name of the r1 contig, the alignment coordinates of the r2 contig edge and the
        r2 alignment edge object

        :param name: name of r1 contig
        :param start: alignment starting coordinate for r2 contig edge
        :param end: alignment ending coordinate for r2 contig edge
        :param ce: ContigEdgeAlignment object
        """
        self.name = name
        self.start = start
        self.end = end
        self.ce_list = [ContigEdgeAlignment(ce_label, orient, start)]

    def is_nearby(self, pc2, min_distance):
        """
        If a second contig edge aligns close to the initial contig edge, expand to coordinates to encompass
        both alignments.

        :param pc2: Comparison PatchedContig object
        :param min_distance: If the two PatchedContig objects are within min_distance, they are merged.
        :return: true if PatchedContig objects are nearby and merged, false if not.
        """
        if self.name == pc2.name and (abs(self.start - pc2.start) < min_distance or abs(self.start - pc2.end) < min_distance or
                                      abs(self.end - pc2.end) < min_distance or abs(self.end - pc2.start) < min_distance):
            self.ce_list.extend(pc2.ce_list)
            self.start = min(self.start, pc2.start)
            self.end = max(self.end, pc2.end)
            return True
        else:
            return False

    def generate_patched_contig(self, r1_contig_dict, r2_contig_dict, edge_length):
        """
        Merge the r1 contig sequence and the two r2 contig sequences

        :param r1_contig_dict: Dictionary with r1 contig sequences
        :param r2_contig_dict: Dictionary with r2 contig sequences
        :param edge_length: Number of bp from r2 contig edge to use for alignment to r1 contigs
        """

        ce_list = sorted(self.ce_list, key=lambda x: x.align_start)

        print("Merging {0} and {1} using sequence from {2}".format(ce_list[0].name, ce_list[1].name, self.name))
        final_seq = ce_list[0].get_seq(r2_contig_dict, edge_length) + \
                    str(r1_contig_dict[self.name][self.start:self.end].seq) + \
                    ce_list[1].get_seq(r2_contig_dict, edge_length)

        record = SeqRecord(Seq(final_seq), id="{0}".format(r2_contig_dict[ce_list[0].name].id, name="", description=""))

        # for every tig involved in the new merge, set to the updated record
        for ce in ce_list:
            for tig in r2_contig_dict[ce.name].id.split("_"):
                r2_contig_dict[tig] = record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('r2_contig_path', help='Path to second haplotype-specific assembly')
    parser.add_argument('r1_contig_path', help='Path to first haplotype-specific assembly')
    parser.add_argument('patched_contig_path', help='Path to patched and trimmed assembly')
    parser.add_argument('minimap2_path', help='Path to minimap2 binary')
    parser.add_argument('-t', '--threads', type=int, default=10, help='Number of minimap2 threads')
    parser.add_argument('-l', '--min_contig_length', type=int, default=20000, help='Minimum contig length')
    parser.add_argument('-e', '--edge_length', type=int, default=5000,
                        help='Number of bp from r2 contig edge to use for alignment to r1 contigs')
    parser.add_argument('-a', '--min_edge_align_length', type=int, default=4000,
                        help='Minimum r2 contig edge to r1 contig alignment length. Shorter alignments are discarded')
    parser.add_argument('-d', '--min_edge_align_dist', type=int, default=2000,
                        help='Minimum r2 contig edge to r1 contig edge aligment distance. Alignments close to the r1 '
                             'contig edge are dropped')
    args = parser.parse_args()

    print('> Starting Contig Patching')

    # check inputs
    for path in [args.r2_contig_path, args.r1_contig_path]:
        if not osu.check_path(path):
            print('Could not find data at specified path: {0}', path)
            sys.exit(1)

    # Create temporary file paths.
    contig_edge_fasta = args.patched_contig_path + '_edges.fasta'
    contig_edge_align = args.patched_contig_path + '_edges.paf'

    print('> Loading contigs')
    r2_contig_dict = SeqIO.to_dict(SeqIO.parse(args.r2_contig_path, "fasta"))
    r1_contig_dict = SeqIO.to_dict(SeqIO.parse(args.r1_contig_path, "fasta"))

    print('> Writing r2 contig edges to file')
    write_contig_ends(r2_contig_dict, contig_edge_fasta, args.min_contig_length, args.edge_length)

    print('> Aligning r2 contig edges to r1 contigs')
    align_contig_ends(args.minimap2_path, contig_edge_fasta, contig_edge_align, args.r1_contig_path)

    print('> Process r2 contig edge alignments')
    potential_patches = load_contig_edge_assembly(contig_edge_align, args.min_contig_length, args.min_edge_align_length,
                                                  args.min_edge_align_dist, args.edge_length)

    print('> Generate patched contigs')
    process_patches(potential_patches, r1_contig_dict, r2_contig_dict, args.edge_length)

    print('> Writing patched contigs to file')
    write_patched_contigs(args.patched_contig_path, r2_contig_dict)

    for path in [contig_edge_fasta, contig_edge_align]:
        osu.remove_file(path)

    print('> Finshed Contig Patching')


def write_trimmed_contigs(trimmed_assembly, tig_dict, min_length):
    """
    Write out trimmed sequences to file

    :param trimmed_assembly: Path to trimmed contig objects
    :param tig_dict: Dictionary of TigInfo objects
    :param min_length: Minimum allow contig length
    """
    with open(trimmed_assembly, 'w') as op_fasta:
        for tig in tig_dict.values():
            if tig.filtered or len(tig.tig_sequence) < min_length:
                continue
            SeqIO.write(tig.get_updated_sequence(), op_fasta, 'fasta')


def write_patched_contigs(patched_assembly, r2_assembly_dict):
    """
    Write out patched contigs to a new file. Patched contig sequences are stored under the names of the joined r2
    contigs, so contig IDs are checked before writing to ensure they are only written once.

    :param patched_assembly: Path to patched assembly fasta
    :param r2_assembly_dict: Dictionary of r2 contigs. Patched contig sequence will replace originals
    :return: Dictionary of patched reads
    """
    patched_contig_set = set()
    with open(patched_assembly, 'w') as opf:
        for contig in r2_assembly_dict.values():
            if contig.id not in patched_contig_set:
                SeqIO.write(contig, opf, 'fasta')
                patched_contig_set.add(contig.id)


def write_contig_ends(r2_contig_dict, contig_edge_fasta, min_contig_length, edge_length):
    """
    Write out r2 contig edge sequences to a fasta file

    :param r2_contig_dict: Dictionary of r2 contig sequences
    :param contig_edge_fasta: Path to contig edge fasta output
    :param min_contig_length: Minimum contig length to process
    :param edge_length: Contig edge length
    """

    with open(contig_edge_fasta, 'w') as op_fasta:
        for record in r2_contig_dict.values():
            if len(record) <= min_contig_length:
                continue
            op_fasta.write(">{0}_start\n{1}\n".format(record.id, str(record.seq)[0:edge_length]))
            op_fasta.write(">{0}_end\n{1}\n".format(record.id, str(record.seq)[-edge_length:]))


def process_patches(potential_patches, r1_contig_dict, r2_contig_dict, edge_length):
    """
    Check each potential patch to see if two different contig ends align nearby. If so, combine the two r2 sequences
    by replacing the alignment end sequences with the corresponding r1 sequence and appending the remaining
    r2 sequences on either end.

    :param potential_patches: List of potential patch locations
    :param r1_contig_dict: Dictionary of r1 contigs
    :param r2_contig_dict: Dictionary of r2 contigs
    :param edge_length: Number of bp from r2 contig edge to use for alignment to r1 contigs
    """
    usable_ends = set()
    for pc in potential_patches:

        # if we find exacly two patch contigs end aligning in the same vicinity of an original contig
        # we can join the patch contigs using the original contig sequence
        failed = False
        if len(pc.ce_list) == 2:

            # Throw an error if the same end of a patch contig is used twice, that would suggest
            # duplicate original contigs
            for ce in pc.ce_list:
                if ce.ce_label in usable_ends:
                    print("Patched contig end overlaps multiple original contigs, can't process", ce.ce_label)
                    failed = True
                usable_ends.add(ce.ce_label)

            # Create the combined sequence and update seq records
            if not failed:
                pc.generate_patched_contig(r1_contig_dict, r2_contig_dict, edge_length)


def load_contig_edge_assembly(contig_edge_align, min_contig_length, min_align_length, min_align_dist, edge_length):
    """
    Load the alignment of each r2 contig edge to the r1 contig sequences.  If the alignment is longer than min_align_length
    and at least min_align_dist away from the r1 contig end, it is treated as a potential patch location.  As each
    potential patch is found, it is compared to previously found patch locations to see if they can be merged.

    :param contig_edge_align: Path to contig edge alignment file in PAF format
    :param min_contig_length: Minimum contig length to process
    :param min_align_length: Minimum alignment length to process
    :param min_align_dist: Mininum alignment distance from r1 contig end to process
    :param edge_length: Contig edge length
    :return: List of PatchedContig objects.
    """
    potential_patches = []

    for align_group in pio.load_query_alignments(contig_edge_align):
        for a in align_group:
            if (a.target_length >= min_contig_length and a.align_length >= min_align_length and
                    a.target_start > min_align_dist and a.target_end < (a.target_length - min_align_dist)):
                pc1 = PatchedContig(a.target_name, a.target_start, a.target_end, a.query_name, a.strand)
                for pc2 in potential_patches:
                    if pc2.is_nearby(pc1, edge_length):
                        break
                else:
                    potential_patches.append(pc1)
    return potential_patches


def align_contig_ends(minimap2_path, contig_edge_fasta, contig_edge_align, r1_assembly):
    """
    Align r2 contig edges to r1 contigs

    :param minimap_path: Path to minimap2 binary
    :param contig_edge_fasta: Path to r2 contig edge fasta file
    :param contig_edge_align: Path to r2 contig edge to r1 contig alignment file
    :param r1_assembly: Path to r1 contig fasta file
    :return:
    """
    result = subprocess.run(
        [minimap2_path, '-x', 'map-hifi', '-c', '-o', contig_edge_align, r1_assembly, contig_edge_fasta])
    if result.returncode != 0:
        print('Minimap2 failed, exiting')
        sys.exit(1)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit(1)

