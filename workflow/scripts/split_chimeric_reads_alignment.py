#!/usr/bin/env python

"""
Raw PacBio reads are split based on alignment to an assembly. Read data for each primary & supplementary alignment is
renamed and restricted to just the aligned portion of the read, effectively splitting reads within the alignment file.
Reads can be restricted to primary only or filtered based on length and alignment score. Finally, supplementary
alignments that align to a different reference than the primary can be filtered.
"""

from collections import defaultdict
import argparse
import sys
import pysam

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_alignment', help='Path to minimap2 alignment file in bam format')
    parser.add_argument('output_alignment', help='Path to split read output file in bam format')
    parser.add_argument('-l', '--min_read_length', help='Minimum alignment length to retain', type=int, default=1000)
    parser.add_argument('-s', '--min_alignment_score', help='Minimum alignment score to retain', type=int, default=1000)
    parser.add_argument('-p', '--primary_only', help='If true, only primary alignments are retained',
                        action="store_true", default=False)
    parser.add_argument('-m', '--match_ref', help='If true, supplementary alignments that align to a different reference '
                                                  'from the primary are dropped', action="store_true", default=False)
    parser.add_argument('-q', '--cap_quality', help='If set, mapping quality is set to this value', type=int, default=-1)
    args = parser.parse_args()

    read_idx = defaultdict(int)
    total = short = poor = written = mismatch_ref = 0

    with pysam.AlignmentFile(args.input_alignment, 'rb') as ipb, \
            pysam.AlignmentFile(args.output_alignment, 'wb', template=ipb) as opb:
        for a in ipb:
            # Skip unaligned, secondary reads, supplementary if set in options
            if a.is_unmapped or a.is_secondary or (args.primary_only and a.is_supplementary):
                continue

            # Skip reads that don't pass quality or length filters
            total += 1
            if a.query_alignment_length < args.min_read_length:
                short += 1
                continue
            if a.get_tag("AS") < args.min_alignment_score:
                poor += 1
                continue

            # Generate split read name
            if a.is_supplementary:
                ref_name = ipb.get_reference_name(a.tid)
                prim_ref_name = a.get_tag("SA").split(",")[0]
                if args.match_ref and ref_name != prim_ref_name:
                    mismatch_ref += 1
                    continue

                read_idx[a.query_name] += 1
                a.query_name = a.query_name + "/" + str(read_idx[a.query_name])
                a.is_supplementary = False
            else:
                a.query_name = a.query_name + "/0"

            # Retain all operations aside from soft-clipping, which represent unaligned regions.  Hard clipping
            # is not expected.
            new_cigar = []
            for op in a.cigartuples:
                if op[0] != 4:
                    new_cigar.append(op)

            # Updated alignment with updated cigar, sequence and qualities
            seq = a.query_alignment_sequence
            qual = a.query_alignment_qualities
            a.query_sequence = seq
            a.query_qualities = qual
            a.cigartuples = new_cigar
            if args.cap_quality != -1:
                a.mapping_quality = args.cap_quality

            opb.write(a)
            written += 1

    pysam.index(args.output_alignment)

    if args.match_ref:
        print("Alignment splitting: Found {0} alignments. {1} were retained, {2} were removed due to length, "
              "{3} were removed due to quality and {4} were removed "
              "for mismatched reference.".format(total, written, short, poor, mismatch_ref))
    else:
        print("Alignment splitting: Found {0} alignments. {1} were retained, {2} were removed due to length, "
              "and {3} were removed due to quality.".format(total, written, short, poor))


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting.")
        sys.exit(1)
