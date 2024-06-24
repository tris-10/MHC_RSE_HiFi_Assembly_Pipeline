import pysam
import os
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="Convert VCF file into lookup table")
    parser.add_argument("input_vcf")
    parser.add_argument("output_table")
    parser.add_argument("--sample",default="SAMPLE")
    args = parser.parse_args()

    with pysam.VariantFile(args.input_vcf, "r") as ipf, open(args.output_table, "w") as opf:
        for record in ipf:
            sd = record.samples[args.sample]
            gt = sd["GT"]

            if gt[0] == 0 and gt[1] == 1:
                base1 = record.ref
                base2 = record.alts[0]
            elif gt[0] == 1 and gt[1] == 0:
                base2 = record.ref
                base1 = record.alts[0]
            else:
                continue

            opf.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(record.id, record.chrom, record.pos, base1, base2))


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("user interrupted, exiting")
        sys.exit(1)
