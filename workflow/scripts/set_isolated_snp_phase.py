#!/usr/bin/env python

"""
Set isolated het positions as phased, with only the SNP in the block. This allows the
haplotagging of reads if the SNP happens to intersect with a microarray het.
"""

import argparse
import pysam
import sys

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_vcf', help='Path to Whatshap phased VCF file')
    parser.add_argument('output_vcf', help='Path to output VCF file')
    parser.add_argument('-d', '--min_dist', help='Minimum distance to nearest SNP for a call to be considered isolated',
                        type=int, default=3000)
    args = parser.parse_args()

    isolated = identify_isolated_hets(args.input_vcf, args.min_dist)
    set_phase(args.input_vcf, args.output_vcf, isolated)
    print("Set phase for {0} isolated hets.".format(len(isolated)))


def set_phase(input_vcf, output_vcf, isolated_hets):
    """
    Read through VCF file and set isolated het positions to be phased and the phase set (PS) to be the
    position of the variant.  Write updated date to new VCF file

    :param input_vcf: Path original phased VCF file
    :param output_vcf: Path to updated VCF file with isolated hets phased
    :param isolated_hets: set of isolated het locations.
    """
    with pysam.VariantFile(input_vcf, 'r') as ipv, pysam.VariantFile(output_vcf, 'w', header=ipv.header) as opv:
        samples = list(ipv.header.samples)
        if len(samples) != 1:
            print('Only one sample is allowed in VCF, found: {0}'.format(len(samples)))
            sys.exit(1)

        for var in ipv:
            loc = "{0}:{1}".format(var.chrom, var.pos)
            if loc in isolated_hets:
                sample_data = var.samples[samples[0]]
                sample_data['GT'] = (0, 1)
                sample_data.phased = True
                sample_data['PS'] = var.pos
                print(loc)
            opv.write(var)

    pysam.tabix_index(output_vcf, preset='vcf', force=True)


def identify_isolated_hets(vcf_file, min_dist):
    """
    Read through phased VCF file and identify unphased heterozygous positions that are min_distance away from
    any other het. Return list of isolated postions as a set.

    :param vcf_file: Path to original phased VCF file
    :param min_dist: Minimum distance from closest het to be considered isolated.
    :return: Set of isolated het positions
    """

    isolated = set()

    with pysam.VariantFile(vcf_file, 'r') as ipv:
        samples = list(ipv.header.samples)
        if len(samples) != 1:
            print('Only one sample is allowed in VCF, found: {0}'.format(len(samples)))
            sys.exit(1)

        prev_pos = 0
        prev_chrom = None
        pot_iso = None

        for var in ipv:
            sample_data = var.samples[samples[0]]
            if not sample_data.phased and sample_data['GT'][0] == 0 and sample_data['GT'][1] == 1:
                if var.pos - prev_pos >= min_dist or prev_chrom != var.chrom:
                    pot_iso = '{0}:{1}'.format(var.chrom, var.pos)
            if pot_iso is not None:
                if var.pos - prev_pos >= min_dist or prev_chrom != var.chrom:
                    isolated.add(pot_iso)
                pot_iso = None
            prev_pos = var.pos
            prev_chrom = var.chrom
        if pot_iso is not None:
            isolated.add(pot_iso)

    return isolated


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("user interrupted, eiting")