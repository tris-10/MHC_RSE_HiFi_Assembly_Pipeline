#!/usr/bin/env python

import sys
import argparse
import pysam
import logging

"""
Restrict VCF to heterozygous variants that are found in a list of known variants and fall outside of the RCCX region.  
Optionally filter based on mapping quality or remove phasing information
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_vcf", help="Path to variant calls in VCF format.")
    parser.add_argument("known_variants", help="Path to a list of valid variants in format chr6-pos-ref-alt.")
    parser.add_argument("output_vcf", help="Path to filtered vcf file")
    parser.add_argument("--min_qual", type=float, default=50, help="Minimum variant quality value")
    parser.add_argument("--min_depth", type=int, default=0, help="Minimum depth of coverage across two alleles")
    parser.add_argument("--min_obs", type=int, default=6, help="Minimum observation of ref/alt allele")
    parser.add_argument("--min_frac", type=float, default=0.12, help="Minimum fraction of ref/alt allele")
    parser.add_argument("--retain_hom", default=False, action="store_true", help="Retain homozygous reference positions")
    parser.add_argument("--sample_name",default="SAMPLE",help="Sample name in VCF header")
    parser.add_argument("--retain_phase", default=False, action="store_true", help="Retain longshot phasing info")
    args = parser.parse_args()

    logging.basicConfig(format='%(levelname)s:%(asctime)s %(message)s', level=logging.INFO)

    comp_start = 31980291
    comp_end = 32046429

    valid_positions = {}
    with open(args.known_variants, "r") as ipf:
        for line in ipf:
            items = line.strip().split("\t")
            valid_positions[items[1]] = items[0]

    with pysam.VariantFile(args.input_vcf) as ip_vcf, pysam.VariantFile(args.output_vcf, "w", header=ip_vcf.header) as op_vcf:
        for r in ip_vcf:
            sd = r.samples[args.sample_name]
            gt = sd['GT']
            ac = r.info['AC']

            if not args.retain_phase: # Strip out phasing information
                sd.phased = False
                sd['PS'] = None

            ref = r.ref.upper()
            alt = r.alts[0].upper()

            key = "{0}-{1}-{2}-{3}".format(r.chrom, r.pos, ref, alt)
            if key in valid_positions: # Make sure variant is in known list
                r.id = valid_positions[key] # Add in variant Id from known list

                if (comp_start < r.pos < comp_end) or r.qual < args.min_qual or sum(ac) < args.min_depth: # Filter all
                    continue
                if gt[0] != gt[1]: # Filter heterozygous
                    if min(ac) < args.min_obs or min(ac) / sum(ac) < args.min_frac:
                        continue
                elif not args.retain_hom: 
                    continue
                op_vcf.write(r)

    pysam.tabix_index(args.output_vcf, preset='vcf', force=True)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logging.error("user interrupted, exiting")
        sys.exit(1)
