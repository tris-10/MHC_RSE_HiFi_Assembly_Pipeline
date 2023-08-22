#!/usr/bin/env python

"""
Contigs are removed if they are supported by a limited number of raw PacBio reads or if they don't align anywhere
inside the targeted region.
"""

from Bio import SeqIO
import argparse
import pysam
import os
import sys
import subprocess
import pandas as pd
import paf_io as pio
import os_utils as osu
import contig_io as cio


class TigInfoCov(cio.TigInfo):
    def __init__(self, tig_name, tig_sequence):
        """
        Initialize TigInfo object.

        :param tig_name: Contig name, which is the ID padded with leading zeros out to 8 positions
        :param tig_sequence: Contig sequence as SeqRecord object
        """

        super().__init__(tig_name, tig_sequence)
        self.read_ids = set()
        self.orig_names = set()

    def add_read(self, read_id, orig_name):
        """
        Store the read_id Canu uses and the raw PacBio read name, prior to splitting.

        :param read_id: Canu read ID
        :param orig_name: raw PacBio read name
        """
        self.orig_names.add(orig_name)
        self.read_ids.add(read_id)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('canu_directory', help='Path to canu directory, as specified when running canu.')
    parser.add_argument('canu_prefix', help='Path to canu output prefix, as sepecified when running canu')
    parser.add_argument('minimap2_path', help='Path to minimap2 executable')
    parser.add_argument('ref_path', help='Path to reference file')
    parser.add_argument('filtered_contig_path', help='Path to filtered contig fasta output file.')
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of minimap2 threads")
    parser.add_argument('-c', '--chrom', default='chr6', help='Target region chromosome')
    parser.add_argument('-s', '--target_start', type=int, default=29654000, help='Target region starting coordinate')
    parser.add_argument('-e', '--target_end', type=int, default=33323000, help='Target region ending coordinate')
    parser.add_argument('-r', '--min_read_support', type=int, default=10,
                        help='Minimum number of raw pacbio reads supporting contig')

    args = parser.parse_args()

    print('>Starting dipoid Filtering')

    # Setup file handles
    canu_path = os.path.join(args.canu_directory, args.canu_prefix)
    contig_path = canu_path + '.contigs.fasta'
    ref_align_paf = canu_path + '_ref.paf'

    # Check inputs
    for path in [contig_path, args.ref_path]:
        if not osu.check_path(path):
            print('Could not find data at specified path: {0}'.format(path))
            sys.exit(1)

    print('>Aligning contigs to reference')
    align_contigs_to_ref(contig_path, ref_align_paf, args.ref_path, args.minimap2_path, args.threads)

    print('>Loading contig sequences')
    tig_dict = create_tiginfocov_dict(contig_path)

    print('>Filter contig by read support')
    filter_by_read_support(canu_path, args.min_read_support, tig_dict)

    print('>Filter contig by location')
    filter_by_location(ref_align_paf, tig_dict, args.chrom, args.target_start, args.target_end)

    print('>Write out passing contigs')
    write_passing_contigs(args.filtered_contig_path, tig_dict)

    # Cleanup
    for path in [ref_align_paf]:
        osu.remove_file(path)

    print('>Ending Diploid filtering')


def create_tiginfocov_dict(tig_fasta):
    """
    Load contig information from contig fasta file, return dictionary of TigInfo objects

    :param tig_fasta: Path to contig file in fasta format
    """
    contig_dict = {}
    with open(tig_fasta, 'r') as ip_fasta:
        for contig in SeqIO.parse(ip_fasta, 'fasta'):
            contig_dict[contig.id] = TigInfoCov(contig.id, contig)
    return contig_dict


def filter_by_read_support(canu_path, min_reads, tig_dict):
    """
    Load the read IDs that make up each contig as reported by Canu.  Look up the original read name, prior
    to splitting, for each ID and count up the number of unique names. If the number of unique names is lower than
    min_reads, the contig is set to be filtered and not written to output.

    :param canu_path: Path to the assembly data
    :param min_reads: Minimum number of raw pacbio reads supporting contig
    :param tig_dict: Dictionary of all contig objects reported by Canu
    """
    read_name_dict = {}
    with open(os.path.join(canu_path + '.seqStore', 'readNames.txt'), 'r') as ipf:
        for line in ipf:
            read_id, read_name = line.strip().split("\t")
            read_name_dict[int(read_id)] = pio.get_original_name(read_name)

    read_to_tig_df = pd.read_csv(os.path.join(canu_path + '.contigs.layout.readToTig'), sep='\t')
    for idx in read_to_tig_df.index:
        tig_name = TigInfoCov.generate_name(read_to_tig_df.iloc[idx]['tigID'])

        if tig_name in tig_dict:
            tig_dict[tig_name].add_read(idx, read_name_dict[read_to_tig_df.iloc[idx]['#readID']])

    for tig in tig_dict.values():
        if len(tig.orig_names) < min_reads:
            print('Removing {0}, only {1} unique reads and {2} total'.format(tig.tig_name, len(tig.orig_names),
                                                                             len(tig.read_ids)))
            tig.filter_contig()


def filter_by_location(ref_align, tig_dict, chrom, target_start, target_end):
    """
    Check contigs to see if they align somewhere in the targeted region. If no part of the contig aligns to the
    targeted region, it is set to be filtered and not written to output.

    :param ref_align: Alignment file in PAF format of contigs to target chromosome
    :param tig_dict: Dictionary of all contig objects reported by Canu
    :param chrom: target region chromosome
    :param target_start: target region starting coordinate
    :param target_end: target region ending coordinate
    """
    valid_tigs = set()
    ref_paf = pio.load_query_alignments(ref_align)
    for align_list in ref_paf:
        for align in align_list:
            if align.target_start < target_end and align.target_end > target_start and align.target_name == chrom:
                valid_tigs.add(align.query_name)

    for tig in tig_dict.values():
        if tig.tig_name not in valid_tigs:
            print('Removing {0}, does not align within the target region'.format(tig.tig_name))
            tig.filter_contig()


def align_contigs_to_ref(contig_path, align_path, ref_path, minimap_path, threads):
    """
    Align contigs to target chromosome using minimap2.

    :param contig_path: Path to the contig fasta
    :param align_path: Path to alignment output file
    :param ref_path: Path to the reference chromosome
    :param minimap_path: Path to the minimap2 executable
    :param threads: Number of threads to use for alignment
    :return: Path to the alignment file in PAF format
    """

    cproc = subprocess.run([minimap_path, '-x', 'map-hifi', '--secondary=no', '-t', str(threads), '-o', align_path,
                            ref_path, contig_path])
    if cproc.returncode != 0:
        print('Error running reference alignment')
        sys.exit(1)


def write_passing_contigs(output_fasta, tig_dict):
    """
    Write trimmed contigs to file is the 'filtered' flag is false.

    :param output_fasta: Path to the filtered fasta output file
    :param tig_dict: Dictionary of all contig objects reported by Canu
    """
    with open(output_fasta, 'w') as opf_fasta:
        for tig_id in sorted(tig_dict.keys()):
            tig = tig_dict[tig_id]
            if not tig.filtered:
                SeqIO.write(tig.get_updated_sequence(), opf_fasta, 'fasta')

    pysam.faidx(output_fasta)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('User interrupted, exiting')
        sys.exit(1)

