#!/usr/bin/env python

"""
Split Longshot VCF by genotype, heterozygous variants into one file, homozygous variants into another.
"""

import pysam
import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_vcf', help='Path to Longshot VCF file')
    parser.add_argument('het_vcf', help='Path to heterozygous variant file output')
    parser.add_argument('hom_vcf', help='Path to homozygous variant file output')
    args = parser.parse_args()

    print('>Starting VCF splitting')

    total = het_count = hom_alt_count = other_count = 0

    with pysam.VariantFile(args.input_vcf, 'r') as ipv, pysam.VariantFile(args.het_vcf, 'w', header=ipv.header) as opv_het, \
            pysam.VariantFile(args.hom_vcf, 'w', header=ipv.header) as opv_hom:

        samples = list(ipv.header.samples)
        if len(samples) != 1:
            print('Only one sample is allowed in VCF, found: {0}'.format(len(samples)))
            sys.exit(1)

        for var in ipv:
            total += 1

            sample_data = var.samples[samples[0]]
            sep = '|' if sample_data.phased else '/'
            gt = sep.join([str(x) for x in sample_data['GT']])

            if gt == '0/1' or gt == '0|1' or gt == '1|0':
                sample_data.phased = False
                sample_data['PS'] = None
                sample_data['GT'] = (0, 1)
                opv_het.write(var)
                het_count += 1
            elif gt == '1/1' or gt == '1|1':
                sample_data.phased = False
                sample_data['PS'] = None
                sample_data['GT'] = (1, 1)
                opv_hom.write(var)
                hom_alt_count += 1
            else:
               other_count += 1

    pysam.tabix_index(args.hom_vcf, preset='vcf', force=True)
    pysam.tabix_index(args.het_vcf, preset='vcf', force=True)

    print('{0} variants found, {1} were hets, {2} were hom. '
          'alt and {3} were other.'.format(total, het_count, hom_alt_count, other_count))
    print('>Finished variant splitting')


if __name__ == "__main__":
    try:

        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit(1)

