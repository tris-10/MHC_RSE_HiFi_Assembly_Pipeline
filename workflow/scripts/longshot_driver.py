#!/usr/bin/env python

"""
Splits up contigs into chunks and runs the Longshot variant caller in parallel. After all chunks are completed, the
resulting VCF files are merged into a single file and then deleted from the filesystem. It is assumed that there is
only one sample in the bam file and the sample name is used in the VCF output. A limited number of Longshot
options are supported on the command line.
"""

from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
import os
import sys
import pysam
import argparse
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bam_path", help='Path to alignment file in BAM format')
    parser.add_argument('ref_path', help='Path contig fasta')
    parser.add_argument('output_prefix', help='Prefix of the output variant file in VCF format')
    parser.add_argument('-p', '--num_proc', help='Number of processes to launch', type=int, default=10)
    parser.add_argument('-c', '--min_alt_count', help='Require a potential SNV to have at least this many alternate '
                                                      'allele observations', type=int, default=10)
    parser.add_argument('-f', '--min_alt_frac', help='Require a potential SNV to have at least this fraction of '
                                                     'alternate allele observations', type=float, default=0.125)
    parser.add_argument('-s', '--chunk_size', help='Variant calling region size', type=int, default=100000)
    parser.add_argument('-r', '--het_rate', help='Specify the heterozygous SNV Rate for genotype prior estimation',
                        type=float, default=0.005)
    parser.add_argument('-q', '--hom_rate', help='Specify the homozygous SNV Rate for genotype prior estimation',
                        type=float, default=0.0001)
    args = parser.parse_args()

    # container for vcf chunks
    vcf_file_dict = {}

    sample_name = extract_sample_name_from_bam(args.bam_path)
    ref_list = reference_chunker(args.ref_path, args.chunk_size)

    with ProcessPoolExecutor(max_workers=args.num_proc) as e:
        for ref in ref_list:
            vcf_file_dict[ref] = e.submit(longshot_launcher, ref, args.bam_path, args.ref_path, args.output_prefix, sample_name,
                              args.min_alt_count, args.min_alt_frac, args.het_rate, args.hom_rate)

    vcf_paths, status = check_longshot_results(ref_list, vcf_file_dict)
    if status == 0:
        print("Concatenating VCF files")
        vcf_concat(vcf_paths, args.output_prefix)
    else:
        print("Longshot command failed")

    print("Cleaning up individual VCF files")
    vcf_cleanup(vcf_paths)
    print("Done")


def extract_sample_name_from_bam(bam_path):
    """
    Extract the first sample name from a bam header. It is assumed that there is only one sample per bam and the
    RG:SM field is set.

    :param bam_path: Path to the alignment file in bam format.
    :return Sample name as a string
    """

    with pysam.AlignmentFile(bam_path, "rb") as ipb:
        return ipb.header.to_dict()['RG'][0]['SM']


def reference_chunker(ref_path, chunk_size):
    """
    Divides each contig into equal-sized regions based on chunk_size. Last chunk per contig will likely be shorter.
    Returns a list of strings in the format contig_name:start-end, which is how Longshot accepts regions on the command
    line.

    :param ref_path: Path the contig fasta file
    :param chunk_size: Size of each variant calling chunk
    :return: List of regions in the format: 'contig_name:start-end'
    """

    ref_list = []
    with open(ref_path, "r") as ipf:
        for record in SeqIO.parse(ipf, "fasta"):
            length = len(record.seq)
            for i in range(0, length, chunk_size):
                ref_list.append("{0}:{1}-{2}".format(record.id, i + 1, min(length, i + chunk_size)))
    return ref_list


def vcf_concat(vcf_paths, output_prefix):
    """
    System call to vcfcombine on an ordered list of VCF files

    :param vcf_paths: list of paths to the VCF file chunks
    :param output_prefix: Prefix of the final output file
    """

    if len(vcf_paths) == 0:
        print("No VCF files to combine")
        return

    out_file = output_prefix + '.vcf'

    cproc = subprocess.run(['bcftools', 'concat', '-o', out_file] + vcf_paths)
    if cproc.returncode == 1:
        print('VCF concat failed')
        sys.exit(1)
    pysam.tabix_index(output_prefix + '.vcf', preset='vcf', force=True)


def vcf_cleanup(vcf_paths):
    """
    Remove VCF file chunks from filesytem.

    :param vcf_paths:  List of paths to the VCF file chunks
    """

    if len(vcf_paths) == 0:
        print('No VCF files to remove')
        return

    cproc = subprocess.run(['rm'] + vcf_paths)
    if cproc.returncode != 0:
        print('VCF cleanup failed, exiting')
        sys.exit(1)


def check_longshot_results(ref_list, vcf_output):
    """
    Generate a list of paths to the VCF file chunks. Longshot calls that completed successfully are added to the
    file list. If any of the results

    :param ref_list: List of variant calling regions
    :param vcf_output: List of Longshot outputs
    :return: tuple with the list of VCF paths and the Longshot completion status 1=failed, 0=passed.
    """

    vcf_paths = []
    status = 0
    for ref in ref_list:
        if vcf_output[ref].result() is None:
            status = 1
        elif os.path.exists(vcf_output[ref].result()):
            vcf_paths.append(vcf_output[ref].result())

    return vcf_paths, status


def longshot_launcher(ref, bam_path, ref_path, out_prefix, sample_name, min_alt_count, min_alt_frac, het_rate, hom_rate):
    """
    Launch Longshot on a reference chunk.

    :param ref: Contig region used for variant calling
    :param bam_path: Path to the alignment file in BAM format
    :param ref_path: Path to the contig sequences in fasta format
    :param out_prefix: Prefix of the VCF chunks
    :param sample_name: Sample name to report in VCF header
    :param min_alt_count: minimum number of alt observations to call variant
    :param min_alt_frac: minimum alt allele fraction to call variant
    :param het_rate: Specify the heterozygous SNV Rate for genotype prior estimation
    :param hom_rate: Specify the homozygous SNV Rate for genotype prior estimation
    :return: Path to VCF output if Longshot passed, otherwise None
    """

    out_name = out_prefix + '_' + ref + '.vcf'
    cproc = subprocess.run(['longshot', '--bam', bam_path, '--ref', ref_path, '--out', out_name, '-F', '--sample_id',
                              sample_name, '--min_alt_count', str(min_alt_count), '--min_alt_frac', str(min_alt_frac),
                              '--het_snv_rate', str(het_rate), '--hom_snv_rate', str(hom_rate), '-r', ref,
                            '-q', '30'])

    if cproc.returncode == 0:
        return out_name
    else:
        return None


if __name__ == "__main__":
    try:

        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit(1)
