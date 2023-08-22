#!/usr/bin/env python

"""
Caps read coverage across contigs to speed up assembly.  Alignments are sorted by decreasing alignment score. The
coverage across the contig is evaluated if the read is included, calculating the number of positions in the read that
get the coverage closer to the cap without going over. If the number of positions is higher than the contribution
threshold, the read is written to a fasta file.
"""
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import argparse
import pysam
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_alignment', help='Haplotype-specific edited reads to contig alignment file, bam format')
    parser.add_argument('output_fasta', help='Reads retained after coverage capping, fasta format')
    parser.add_argument('-d', '--depth', help='Target depth of coverage', type=int, default=80)
    parser.add_argument('-q', '--min_qual', help='Minimum alignment quality required to retain read', type=int, default=1)
    parser.add_argument('-c', '--min_contrib', help='Minimum bp under depth required to retain read', type=int,
                        default=100)
    parser.add_argument('-t', '--sort_threads', help='Number of sorting threads used by samtools', type=int, default=5)
    args = parser.parse_args()

    nas_bam = args.input_alignment + '.n_as.bam'
    sas_bam = args.input_alignment + '.s_as.bam'

    print('>Negating AS')
    negate_alignment_score(args.input_alignment, nas_bam)

    print('>Sorting alignments by decreasing AS')
    sort_alignment_by_align_score(nas_bam, sas_bam, args.sort_threads)

    print('>Capping coverage')
    cap_bam_coverage(sas_bam, args.output_fasta, args.depth, args.min_qual, args.min_contrib)

    print('>Deleting temp files')
    for temp in [nas_bam, sas_bam]:
        if os.path.exists(temp):
            os.remove(temp)

def negate_alignment_score(input_align, output_align):
    """
    Negate the alignment score so that samtools will sort reads in decreasing alignment score order

    :param input_align: Haplotype-specific edited reads to contig alignment file, bam format
    :param output_align: Path to temporary alignment file with negated alignment scores
    """
    with pysam.AlignmentFile(input_align, 'rb') as ip_bam, pysam.AlignmentFile(output_align, 'wb', template=ip_bam) as op_bam:
        for align in ip_bam:
            if align.is_unmapped:
                continue
            align.set_tag('AS', -align.get_tag('AS'))
            op_bam.write(align)


def sort_alignment_by_align_score(input_align, output_align, threads):
    """
    Sort by AS tag using samtools sort

    :param input_align: Path to temporary alignment file with negated alignment scores
    :param output_align: Path to temporary alignment file sorted by alignment score
    :param threads: number of threads to use when sorting
    """

    cproc = subprocess.run(['samtools', 'sort', '-@', str(threads), '-t', 'AS', '-o', output_align, input_align])

    if cproc.returncode != 0:
        print('Samtools sort failed, exiting')
        sys.exit(1)


def cap_bam_coverage(input_align, output_fasta, max_depth, min_qual, min_contrib):
    """
    Count how many positions of the alignment are under the coverage cap. If the alignment has meaningful contribution
    write to file and increment coverage.

    :param input_align: Path to temporary alignment file sorted by alignment score
    :param output_fasta: Path to capped read output fasta file
    :param max_depth: Target read depth across contig
    :param min_qual: Minimum alignment quality
    :param min_contrib: Minimum bp under depth required to retain read
    """
    depth_dict = {}

    with pysam.AlignmentFile(input_align, 'rb') as ip_bam, open(output_fasta, 'w') as op_fasta:
        print(">Creating depth array")
        for ref_name in ip_bam.references:
            ref_tid = ip_bam.get_tid(ref_name)
            length = ip_bam.get_reference_length(ref_name)
            depth_dict[ref_tid] = [0] * length

        total_reads = 0
        retained_reads = 0
        for align in ip_bam:
            if align.mapping_quality < min_qual or align.is_unmapped or align.is_supplementary or align.is_secondary:
                continue

            if total_reads % 10000 == 0 and total_reads != 0:
                print('Processing: ', total_reads, retained_reads, round(retained_reads / total_reads * 100, 2))

            # Check to see how many positions are useful in the read
            contrib = 0
            for i in range(align.reference_start, align.reference_end):
                if depth_dict[align.tid][i] < max_depth:
                    contrib += 1

            # Write out read if it contributes meaningfully to coverage
            if contrib > min_contrib:
                for i in range(align.reference_start, align.reference_end):
                    depth_dict[align.tid][i] += 1

                seq = Seq(align.query_sequence)
                if align.is_reverse:
                    seq = seq.reverse_complement()

                record = SeqRecord(seq, id=align.query_name, description='')
                SeqIO.write(record, op_fasta, 'fasta')
                retained_reads += 1
            total_reads += 1

        print('Finished: ', total_reads, retained_reads, round(retained_reads / total_reads * 100, 2))


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit(1)