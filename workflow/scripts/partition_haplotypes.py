#!/usr/bin/env python

"""
Split reads by haplotype
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import argparse
import copy
import re
import sys
import pysam

import paf_io as pio
from partition_classes import PhaseBlock
from partition_classes import ProbeAlign
from partition_classes import ProbeInfo
from partition_classes import InfoLogger

LOGGER = None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('contig_fasta', help='Path to contig fasta file')
    parser.add_argument('tagged_bam', help='Path to the haplotagged bam file')
    parser.add_argument("primary_bam", help='Path to primary bam file')
    parser.add_argument('probe_file', help='Phased microarray information')
    parser.add_argument('phased_vcf', help='Path to whatshap phased VCF file')
    parser.add_argument('phase_blocks_gtf', help='Path to whatshap phased blocks')
    parser.add_argument('front_probe_bam', help='Path to front probe alignment file')
    parser.add_argument('back_probe_bam', help='Path to back probe alignment file')
    parser.add_argument('output_prefix', help="Prefix of the output files")
    parser.add_argument('-m', '--min_gt_count', type=int, default=1,
                        help='Minimum intersecting microarray positions required to assign phased block')
    parser.add_argument('-b', '--max_frac_contam', type=float, default=0.25,
                        help='Maximum allowed minor haplotype fraction to assign phased block')
    parser.add_argument('-c', '--max_count_contam', type=int, default=5,
                        help='Maximum number of minor haplotype observations to assign phased block')
    parser.add_argument('-e', '--edge_distance', type=int, default=1000,
                        help='Number of basepairs from contig edge to distrust when calling haploid blocks. '
                             'Genotypes are dropped if they are inconsistent with internal positions.')
    parser.add_argument('-d', '--min_haploid_depth', type=int, default=30,
                        help='Minimum depth required to call a haploid position')
    parser.add_argument('-r', '--min_recover_depth', type=int, default=20,
                        help='Minimum observations required for phased block recovery')
    parser.add_argument('-f', '--max_recover_frac', type=float, default=0.25,
                        help='Maximum allowed minor haplotype fraction to recover phased block')
    parser.add_argument('-l', '--log_file', default=sys.stdout)
    parser.add_argument('-q', '--silent', default=False, action='store_true',
                        help='Detailed partitioning information is not logged')
    args = parser.parse_args()

    global LOGGER
    LOGGER = InfoLogger(args.log_file, args.silent)

    print(">Starting")
    snp_coverage = defaultdict(int)

    print('>Loading contig info')
    contig_data = load_contig_data(args.contig_fasta)

    print('>Loading phased block coordinates')
    phased_blocks = load_phased_blocks(args.phase_blocks_gtf)

    print('>Loading unphased block coordinates')
    unphased_blocks = generate_unphased_blocks(phased_blocks, contig_data)

    print('>Loading phased microarray hets information')
    probe_info_dict = load_probe_data(args.probe_file)

    print('>Loading microarray probe alignment data')
    front_align_dict = load_probe_alignment(args.front_probe_bam, probe_info_dict, True)
    back_align_dict = load_probe_alignment(args.back_probe_bam, probe_info_dict, False)

    print('>Identifying valid alignments and determine contig haplotype at probe')
    valid_probe_dict = identify_valid_probe_locations(front_align_dict, back_align_dict, contig_data)

    print('>Intersect VCF with probe locations in phased blocks')
    unphased_hets = add_probes_to_phased_blocks(args.phased_vcf, args.primary_bam, phased_blocks, valid_probe_dict, snp_coverage)

    print('>Intersect VCF with probe locations in unphased blocks')
    add_probes_to_unphased_blocks(contig_data, args.primary_bam, unphased_blocks, valid_probe_dict, snp_coverage,
                                  args.min_haploid_depth, args.edge_distance)

    print('>Check for inconsistent phasing')
    phased_blocks = split_phased_blocks(phased_blocks, args.min_gt_count, args.max_count_contam, args.max_frac_contam)

    print('>Classifying reads in phased blocks')
    read_assignments, unknown_unphased_blocks = assign_reads_in_phased_blocks(args.tagged_bam, phased_blocks,
                                                                              args.min_gt_count, args.max_count_contam,
                                                                              args.max_frac_contam)

    print('>Classifying reads in unphased blocks')
    read_assignments = assign_reads_in_unphased_blocks(args.tagged_bam, unphased_blocks, unphased_hets,
                                                       read_assignments, contig_data, args.max_count_contam,
                                                                              args.max_frac_contam)

    print('>Recover unknown phase blocks')
    read_assignments = recover_unknown_phased_blocks(read_assignments, unknown_unphased_blocks, args.min_recover_depth,
                                                     args.max_recover_frac)

    print('>Write out fasta files')
    write_split_fasta(args.tagged_bam, read_assignments, args.output_prefix)

    print("Complete")


def load_contig_data(contig_file):
    """
    Load in contig sequences

    :param contig_file: Path to contig fasta
    :return: dictionary of contig sequence object, keyed by name
    """

    contig_data = SeqIO.to_dict(SeqIO.parse(contig_file, 'fasta'))
    LOGGER.log('{0} contigs found'.format(len(contig_data.keys())))
    return contig_data


def load_phased_blocks(block_file):
    """Load in phased block data from whatshap, gtf format"""

    phase_blocks = {}
    with open(block_file, 'r') as ipf:
        for line in ipf:
            items = line.strip().split('\t')
            m = re.match('gene_id "(\\w+)"', items[8])
            if m is None:
                print("Could not find gene id {0}".format(items[8]))
                sys.exit(1)
            phase_set_name = items[0] + ":" + m.group(1)
            phase_blocks[phase_set_name] = PhaseBlock(phase_set_name, items[0], int(items[3]), int(items[4]) + 1)

    LOGGER.log('{0} phased blocks found'.format(len(phase_blocks.keys())))
    return phase_blocks


def generate_unphased_blocks(phased_blocks, contig_file):
    """Generate unphased blocks that sit in-between phased blocks"""

    unphased_blocks = {}
    for contig in contig_file.keys():
        start = 0
        for block in sorted((b for b in phased_blocks.values() if b.contig == contig), key=lambda x: start):
            end = block.start
            name = "{0}:{1}-{2}".format(contig, start, end)
            unphased_blocks[name] = PhaseBlock(name, contig, start, end)
            start = block.end
        name = "{0}:{1}-{2}".format(contig, start, len(contig_file[contig]))
        unphased_blocks[name] = PhaseBlock(name, contig, start, len(contig_file[contig]))
    LOGGER.log('{0} unphased blocks found'.format(len(unphased_blocks.keys())))
    return unphased_blocks


def load_probe_data(probe_file):
    """
    Load phased microarray snp information from tab-delim text file
    ID\tCHR\tPOS\tBASE1\tBASE2\n

    :param probe_file: Path to phased het file
    :return: dictionary of ProbeInfo objects, keyed on SNP ID
    """

    probe_info = {}
    with open(probe_file, "r") as ipf:
        for line in ipf:
            items = line.strip().split('\t')
            probe_info[items[0]] = ProbeInfo(items[1], items[3], items[4])

    LOGGER.log('{0} total heterozygous microarray positions'.format(len(probe_info.keys())))
    return probe_info


def load_probe_alignment(alignment_file, probe_info, is_front):
    """Load probe alignment data """

    total_alignments = valid_alignments = 0
    probe_align_dict = {}

    with pysam.AlignmentFile(alignment_file, "r") as ip_bam:
        for a in ip_bam:
            total_alignments += 1

            # skip unaligned probes or probes with trimming
            if a.is_unmapped or len([cigar for cigar in a.cigartuples if cigar[0] == 5 or cigar[0] == 4]) > 0:
                continue

            # Drop probes with more than two mismatches, unless spike-in.
            if probe_info[a.query_name].spike:
                if a.get_tag("NM") > 1:
                    continue
            elif a.get_tag("NM") > 2:
                continue
            valid_alignments += 1

            # calculate the position of the SNP on the contig based on front/back probe and alignment orientation
            if is_front:
                pos = a.reference_start if a.is_reverse else a.reference_end + 1
            else:
                pos = a.reference_end + 1 if a.is_reverse else a.reference_start

            # create probe alignment location object
            contig = ip_bam.get_reference_name(a.reference_id)
            key = '{0}:{1}'.format(contig, pos)
            probe_align = ProbeAlign(a.query_name, pos, ip_bam.get_reference_name(a.reference_id), a.is_reverse,
                                     copy.copy(probe_info[a.query_name]))
            probe_align_dict[key] = probe_align

            if probe_info[a.query_name].spike:
                LOGGER.log('SPIKE: Contig: {0} Location: {1} label: {2}'.format(contig, pos, a.query_name))

    LOGGER.log('{0} probes alignments found, {1} passed filters'.format(total_alignments, valid_alignments))
    return probe_align_dict


def identify_valid_probe_locations(front_probe_dict, back_probe_dict, contig_data):
    """Compare locations of the front and back probe to determine if the probe alignment location is acceptable
    and also determines if the reference sequence matches hap1 or hap2"""

    hap1_cnt = hap2_cnt = unk_cnt = valid_cnt = 0
    valid_probe_dict = {}

    for key, probe in front_probe_dict.items():
        if not probe.probe_info.spike and (key not in back_probe_dict or probe.is_reverse != back_probe_dict[key].is_reverse):
            continue
        valid_cnt += 1

        base = contig_data[probe.contig_name][probe.contig_pos - 1]
        if base == probe.get_base(True):
            probe.set_haplotype('HAP1')
            hap1_cnt += 1
        elif base == probe.get_base(False):
            probe.set_haplotype('HAP2')
            hap2_cnt += 1
        else:
            probe.set_haplotype('NONE')
            unk_cnt += 1

        valid_probe_dict[key] = probe

    LOGGER.log_sect('Found {0} valid probe alignment locations on contigs'.format(valid_cnt))
    LOGGER.log_sect('{0} probe pairs identified hap1 base on contig'.format(hap1_cnt))
    LOGGER.log_sect('{0} probe pairs identified hap2 base on contig'.format(hap2_cnt))
    LOGGER.log_sect('{0} probe pairs identified unexpected base on contig'.format(unk_cnt))
    LOGGER.flush_sect('Merge probe data stats')

    return valid_probe_dict


def add_probes_to_phased_blocks(phased_vcf, primary_bam, phase_blocks, valid_probe_data, cov_dict):
    """Read through phased SNPs called using pacbio data aligned to CANU contigs. Store phasing information"""

    phased_cnt = usable_cnt = hap1_cnt = hap2_cnt = unk_cnt = 0
    unphased_hets = defaultdict(set)

    last_contig = None
    pileup = None

    with pysam.VariantFile(phased_vcf) as ip_vcf, pysam.AlignmentFile(primary_bam, 'rb') as ip_bam:
        samples = list(ip_vcf.header.samples)
        for v in ip_vcf:
            key = '{0}:{1}'.format(v.contig, v.pos)

            # Load pileup for each newly encountered contig
            if v.contig != last_contig:
                LOGGER.log('Loading pileups for {0}'.format(v.contig))
                pileup = ip_bam.pileup(v.contig)
                last_contig = v.contig

            # Extract genotype information
            sample_data = v.samples[samples[0]]
            sep = '|' if sample_data.phased else '/'
            gt = sep.join([str(g) for g in sample_data['GT']])

            # This is simply to annotate the number of variants in the block, no matter if there is microarray data
            if gt == '0|1' or gt == '1|0':
                phased_cnt += 1
                phase_set = '{0}:{1}'.format(v.contig, sample_data['PS'])
                phase_blocks[phase_set].add_vcf_pos(v.pos)

            # only process positions found in MA data
            if key in valid_probe_data:
                probe = valid_probe_data[key]
                cov_dict[probe.snp_id] += check_coverage(pileup, probe.contig_pos)

                # if a position is unphased, we move on.  However, we store unphased hets to try and avoid
                # overcalling hemizygous regions
                if not sample_data.phased:
                    unphased_hets[v.contig].add(v.pos)
                    continue
                usable_cnt += 1
                phase_set = '{0}:{1}'.format(v.contig, sample_data['PS'])
                ref, alt = (v.ref, v.alts[0]) if gt == '0|1' else (v.alts[0], v.ref)

                assign = None
                if ref == probe.get_base(True) and alt == probe.get_base(False):
                    assign = 'HAP1'
                    hap1_cnt += 1
                elif ref == probe.get_base(False) and alt == probe.get_base(True):
                    assign = 'HAP2'
                    hap2_cnt += 1
                else:
                    unk_cnt += 1

                if assign:
                    probe.set_haplotype(assign)
                    phase_blocks[phase_set].add_snp(probe)

    LOGGER.log_sect('{0} total variant calls in VCF that were phased by whatshap'.format(phased_cnt))
    LOGGER.log_sect('{0} intersecting variant calls in VCF that were phased by whatshap'.format(usable_cnt))
    LOGGER.log_sect('{0} intersecting variant calls phased with hap1 as reference base'.format(hap1_cnt))
    LOGGER.log_sect('{0} intersecting variant calls phased with hap2 as reference base'.format(hap2_cnt))
    LOGGER.log_sect('{0} intersecting variant calls phased with unknown reference base'.format(unk_cnt))
    LOGGER.flush_sect('Phased SNP Assignment Results')
    return unphased_hets


def add_probes_to_unphased_blocks(contig_dict, primary_bam, unphased_blocks, valid_probe_dict, cov_dict, min_depth,
                                  min_distance):
    """This script adds probes to unphased blocks"""

    last_contig = None
    pileup = None
    used = dropped = 0

    with pysam.AlignmentFile(primary_bam, 'rb') as ip_bam:
        for block in sorted(unphased_blocks.values(), key=lambda x: (x.contig, x.start)):
            if block.contig != last_contig:
                LOGGER.log('Loading pileups for {0}'.format(block.contig))
                pileup = ip_bam.pileup(block.contig)
                last_contig = block.contig

            for probe in valid_probe_dict.values():
                if probe.contig_name == block.contig and block.start <= probe.contig_pos < block.end:
                    # Check to see if position is close to contig edge
                    if (probe.contig_pos < min_distance or
                            probe.contig_pos > (len(contig_dict[probe.contig_name]) - min_distance)):
                        probe.set_edge()

                    # Check coverage
                    cov = cov_dict[probe.snp_id] + check_coverage(pileup, probe.contig_pos)
                    if cov >= min_depth:
                        block.add_snp(probe)
                        used += 1
                    else:
                        LOGGER.log('Uncalled het {0} removed from contig {1} position {2} with '
                                   'global coverage {3}'.format(probe.snp_id, probe.contig_name, probe.contig_pos, cov))
                        dropped += 1

    LOGGER.log_sect('{0} probes retained in unphased blocks'.format(used))
    LOGGER.log_sect('{0} probes dropped from unphased blocks'.format(dropped))
    LOGGER.flush_sect('Unphased SNP Assignment Results')


def split_phased_blocks(phased_blocks, min_count, max_contam_count, max_contam_frac):
    updated_blocks = []
    for block in phased_blocks.values():
        result = block.classify_block(min_count, max_contam_count, max_contam_frac)
        if result == 'MIXED':
            LOGGER.log("Splitting block: {0}".format(block.name))
            LOGGER.log(block.genotype_string)
            updated = block.split_block(10)
            if len(updated) > 0:
                LOGGER.log("Generated {0} blocks".format(len(updated)))
                updated_blocks.extend(updated)
            else:
                updated_blocks.extend(block)
        else:
            updated_blocks.append(block)
    return updated_blocks


def assign_reads_in_phased_blocks(tagged_bam, phased_blocks, min_count, max_count_contam, max_frac_contam):
    """This function iterates the alignments in a phased regions and assigns a haplotype to each read"""

    # For each phase block, grab member reads
    read_assignments = {}
    unknown_unphased_blocks = {}

    reads_in_phased_blocks = 0
    hap1_reads_in_phased_blocks = 0
    hap2_reads_in_phased_blocks = 0
    unknown_reads_in_phased_blocks = 0
    unknown_reads_at_single_snps = 0
    dropped_reads_in_phased_blocks = 0

    phased_snps = 0
    unphased_single_snps = 0
    unphased_block_snps = 0

    with pysam.AlignmentFile(tagged_bam, 'rb') as ip_bam:
        for block in sorted(phased_blocks, key=lambda x: (x.contig, x.start)):
            hap = block.classify_block(min_count, max_count_contam, max_frac_contam)
            phased_snp_count = len(block.vcf_pos)

            # Count number of assigned VCF variants
            if hap in ['HAP1_CLEAN', 'HAP2_CLEAN', 'HAP1_CONTAM', 'HAP2_CONTAM']:
                phased_snps += phased_snp_count
            elif len(block.vcf_pos) == 1:
                unphased_single_snps += 1
            else:
                unphased_block_snps += phased_snp_count

            # Create a container for blocks that don't intersect with MA data
            if hap == 'UNKNOWN':
                unknown_unphased_blocks[block.name] = []

            b_tot_cnt = b_assign_cnt = b_drop_cnt = b_unk_cnt = 0
            for read in ip_bam.fetch(block.contig, block.start, block.end):
                try:
                    # Skip if read has already been assigned or if read isn't primary
                    if read.query_name in read_assignments or read.is_supplementary or read.is_secondary:
                        continue

                    b_tot_cnt += 1
                    reads_in_phased_blocks += 1
                    tag = read.get_tag('HP')

                    if hap in ['HAP1_CLEAN', 'HAP1_CONTAM']:
                        b_assign_cnt += 1
                        if tag == 1:
                            read_assignments[read.query_name] = 'HAP1'
                            hap1_reads_in_phased_blocks += 1
                        else:
                            read_assignments[read.query_name] = 'HAP2'
                            hap2_reads_in_phased_blocks += 1
                    elif hap in ['HAP2_CLEAN', 'HAP2_CONTAM']:
                        b_assign_cnt += 1
                        if tag == 1:
                            read_assignments[read.query_name] = 'HAP2'
                            hap2_reads_in_phased_blocks += 1
                        else:
                            read_assignments[read.query_name] = 'HAP1'
                            hap1_reads_in_phased_blocks += 1
                    else:
                        read_assignments[read.query_name] = 'UNPHASED'
                        b_unk_cnt += 1
                        if phased_snp_count < 2:
                            unknown_reads_at_single_snps += 1
                        else:
                            unknown_reads_in_phased_blocks += 1

                        if hap == 'UNKNOWN':
                            unknown_unphased_blocks[block.name].append((read.query_name, tag))
                except KeyError:
                    b_drop_cnt += 1
                    read_assignments[read.query_name] = 'DROP'
                    dropped_reads_in_phased_blocks += 1

            if hap == 'HAP1_CLEAN':
                call = 'Assigning phased block HAP1 clean.'
            elif hap == 'HAP1_CONTAM':
                call = 'Assigning phased block HAP1 w/ contamination.'
            elif hap == 'HAP2_CLEAN':
                call = 'Assigning phased block HAP2 clean.'
            elif hap == 'HAP2_CONTAM':
                call = 'Assigning phased block HAP2 w/ contamination.'
            elif hap == 'UNKNOWN':
                call = 'Ignoring unknown phased block.'
            elif hap == 'MIXED':
                call = 'Ignoring mixed phased block.'

            call_string = ''
            if hap in ['HAP1_CONTAM', 'HAP2_CONTAM', 'MIXED']:
                call_string = block.genotype_string
            LOGGER.log(call + " {0} has {1} reads across {9} SNPs: {2} tagged {3} untagged and {4} unknown. "
                              "Call based on {5} genotypes ({6}/{7}). {8}".format(block.name, b_tot_cnt, b_assign_cnt,
                                                                                  b_drop_cnt, b_unk_cnt,
                                                                                  len(block.snp_list), block.hap1_count,
                                                                                  block.hap2_count,
                                                                                  call_string, phased_snp_count))

    LOGGER.log_sect('{0} reads phased blocks'.format(reads_in_phased_blocks))
    LOGGER.log_sect('{0} reads are assigned to hap1'.format(hap1_reads_in_phased_blocks))
    LOGGER.log_sect('{0} reads are assigned to hap2'.format(hap2_reads_in_phased_blocks))
    LOGGER.log_sect('{0} reads are within phased blocks, but are not tagged by whatshap'.format(dropped_reads_in_phased_blocks))
    LOGGER.log_sect('{0} reads are within phased blocks, but the block can not be assigned to '
                   'paternal/maternal using genotyping'.format(unknown_reads_in_phased_blocks))
    LOGGER.log_sect('{0} reads are at a single SNP, but the block cannot be assigned to '
                   'paternal/maternal using genotyping'.format(unknown_reads_at_single_snps))
    LOGGER.log_sect('{0} SNPs in phased blocks'.format(phased_snps + unphased_block_snps + unphased_single_snps))
    LOGGER.log_sect('{0} SNPs in assigned blocks'.format(phased_snps))
    LOGGER.log_sect('{0} SNPs in unassigned blocks'.format(unphased_block_snps))
    LOGGER.log_sect('{0} SNPs in single-snp blocks'.format(unphased_single_snps))
    LOGGER.flush_sect('Read Assignment Stats Phased Blocks')

    return read_assignments, unknown_unphased_blocks


def assign_reads_in_unphased_blocks(tagged_bam, unphased_blocks, unphased_hets, read_assignments, contig_data,
                                    max_count_contam, max_frac_contam):
    """This function iterates over the reads in unphased blocks and assigns haplotypes"""

    reads_in_unphased_blocks = 0
    hap1_reads_in_unphased_blocks = 0
    hap2_reads_in_unphased_blocks = 0
    unknown_reads_in_unphased_blocks = 0
    reads_in_short_blocks = 0

    with pysam.AlignmentFile(tagged_bam, 'rb') as ip_bam:
        for block in sorted(unphased_blocks.values(), key=lambda x: (x.contig, x.start)):
            contig_length = len(contig_data[block.contig])

            # Don't process short blocks
            if block.length < 1000:
                total = 0
                for read in ip_bam.fetch(block.contig, block.start, block.end):
                    if read.query_name in read_assignments and read_assignments[read.query_name] != "UNKNOWN":
                        continue
                    read_assignments[read.query_name] = 'UNKNOWN'
                    reads_in_unphased_blocks += 1
                    reads_in_short_blocks += 1
                    total += 1
                LOGGER.log('Skipping short unphased block. {0} with {1} reads.'.format(block.name, total))
                continue

            # edges of contigs are more likely to be haploid, so use more aggressive settings
            if block.length == contig_length or block.start == 0 or block.end == contig_length:
                location = 'edge'
                if block.length >= 50000:
                    hap = block.classify_block(2, max_count_contam, max_frac_contam, unphased_hets)
                else:
                    hap = block.classify_block(1, max_count_contam, max_frac_contam, unphased_hets)
            else:
                location = 'middle'
                hap = block.classify_block(5, max_count_contam, max_frac_contam, unphased_hets)

            assign_cnt = unk_cnt = tot_cnt = 0
            for read in ip_bam.fetch(block.contig, block.start, block.end):
                tot_cnt += 1

                # Skip over already assigned reads
                if read.query_name in read_assignments and (read_assignments[read.query_name] == 'HAP1'
                                                            or read_assignments[read.query_name] == 'HAP2'
                                                            or read_assignments[read.query_name] == 'HAP1_HAPLOID'
                                                            or read_assignments[read.query_name] == 'HAP2_HAPLOID'):
                    continue

                reads_in_unphased_blocks += 1
                if hap in ['HAP1_CLEAN', 'HAP1_CONTAM']:
                    assign_cnt += 1
                    read_assignments[read.query_name] = 'HAP1_HAPLOID'
                    hap1_reads_in_unphased_blocks += 1
                elif hap in ['HAP2_CLEAN', 'HAP2_CONTAM']:
                    assign_cnt += 1
                    read_assignments[read.query_name] = 'HAP2_HAPLOID'
                    hap2_reads_in_unphased_blocks += 1
                else:
                    unk_cnt += 1
                    unknown_reads_in_unphased_blocks += 1
                    read_assignments[read.query_name] = 'UNKNOWN'

            if hap == "UNPHASED":
                call = "Ignoring haploid het unphased block."
            elif hap == "HAP1_CLEAN":
                call = "Assigning haploid unphased block HAP1 clean"
            elif hap == "HAP1_CONTAM":
                call = "Assigning haploid unphased block HAP1 w/ contamination."
            elif hap == "HAP2_CLEAN":
                call = "Assigning haploid unphased block HAP2 clean."
            elif hap == "HAP2_CONTAM":
                call = "Assigning haploid unphased block HAP2 w/ contamination."
            elif hap == "UNKNOWN":
                call = "Ignoring haploid unknown unphased block."
            elif hap == "MIXED":
                call = "Ignoring haploid multi-mixed unphased block."

            genotype_string = ""
            edge_reduced = " "
            if hap in ['HAP1_CONTAM', 'HAP2_CONTAM', 'MIXED']:
                genotype_string = block.genotype_string
            if block.edge_reduced:
                edge_reduced = ' edge reduced '
            LOGGER.log(call + ' {0} loc: {8} length: {9} has {1} reads, {2} assigned and {3} unknown.  '
                              'Call based on {4}{10}genotypes ({5}/{6}). {7}'.format(block.name, tot_cnt, assign_cnt, unk_cnt,
                                                                                     len(block.snp_list), block.hap1_count,
                                                                                     block.hap2_count, genotype_string, location,
                                                                                     block.length, edge_reduced))

    LOGGER.log_sect('{0} reads in unphased blocks'.format(reads_in_unphased_blocks))
    LOGGER.log_sect('{0} reads are assigned to hap1'.format(hap1_reads_in_unphased_blocks))
    LOGGER.log_sect('{0} reads are assigned to hap2'.format(hap2_reads_in_unphased_blocks))
    LOGGER.log_sect('{0} reads can not be assigned >500bp'.format(unknown_reads_in_unphased_blocks))
    LOGGER.log_sect('{0} reads can nott be assigned <500bp'.format(reads_in_short_blocks))
    LOGGER.flush_sect('Read Assignment Stats Unphased Blocks')

    return read_assignments

def recover_unknown_phased_blocks(read_assignments, unknown_phased_blocks, min_recov_depth, max_recov_frac):
    for unknown_block_name in unknown_phased_blocks.keys():
        unknown_block_reads = unknown_phased_blocks[unknown_block_name]
        hap1_hp1 = hap1_hp2 = hap2_hp1 = hap2_hp2 = total = 0

        lookup = {}

        # For each unknown read, store HP under the original pacbio name
        for read in unknown_block_reads:
            if read[1] == 1 or read[1] == 2:
                lookup[pio.get_original_name(read[0])] = read[1]
                total += 1

        # look at assignments, see if the read was assigned elsewhere
        for read_name, read_class in read_assignments.items():
            orig_name = pio.get_original_name(read_name)
            if orig_name in lookup:
                if read_class in ['HAP1', 'HAP1_HAPLOID']:
                    if lookup[orig_name] == 1:
                        hap1_hp1 += 1
                    elif lookup[orig_name] == 2:
                        hap1_hp2 += 1
                elif read_class in ['HAP2', 'HAP2_HAPLOID']:
                    if lookup[orig_name] == 1:
                        hap2_hp1 += 1
                    elif lookup[orig_name] == 2:
                        hap2_hp2 += 1

        same = hap1_hp1 + hap2_hp2
        cross = hap1_hp2 + hap2_hp1
        per = max(same, cross) / (same + cross) if (same + cross) > min_recov_depth else 0

        # If there are enough assigned reads and they have consistent category
        if per > (1-max_recov_frac):
            for read in unknown_block_reads:
                if read[1] == 1:
                    read_assignments[read[0]] = 'HAP1' if same > cross else 'HAP2'
                elif read[1] == 2:
                    read_assignments[read[0]] = 'HAP2' if same > cross else 'HAP1'

            LOGGER.log('Recovering phased block {0} with {1} phased reads. Percent {2}: {3} {4} {5} {6}'.format(
                unknown_block_name, total, round(per * 100, 2), hap1_hp1, hap1_hp2, hap2_hp1, hap2_hp2))
        else:
            LOGGER.log('Could not recover phased block {0} with {1} phased reads. Percent {2}: {3} {4} {5} {6}'.format(
                unknown_block_name, total, round(per * 100, 2), hap1_hp1, hap1_hp2, hap2_hp1, hap2_hp2))
    return read_assignments


def check_coverage(pileups, position):
    for pc in pileups:
        if pc.pos == position:
            return pc.n
        elif pc.pos > position:
            return 0
    return 0

def write_split_fasta(tagged_bam, collapsed_read_assignments, output_prefix):
    """Given a list of PacBio read assignments, read through bam file and write to a new file based on assignment"""
    total_written = total_reads = 0
    hap1_diploid = hap1_haploid = 0
    hap2_diploid = hap2_haploid = 0
    dropped = unphased = unobserved = unknown = 0

    with pysam.AlignmentFile(tagged_bam, 'rb') as ip_bam, open(output_prefix + '_hap1_only.fasta', 'w') as hap1_out, \
            open(output_prefix + '_hap2_only.fasta', 'w') as hap2_out, open(output_prefix + '_unknown.fasta', 'w') as unknown_out:
        for a in ip_bam:
            total_reads += 1
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue
            total_written += 1
            record = SeqRecord(Seq(a.query_alignment_sequence), id=a.query_name)
            if a.query_name in collapsed_read_assignments:
                call = collapsed_read_assignments[a.query_name]
                if call == 'HAP1':
                    SeqIO.write(record, hap1_out, 'fasta')
                    hap1_diploid += 1
                elif call == 'HAP1_HAPLOID':
                    SeqIO.write(record, hap1_out, "fasta")
                    hap1_haploid += 1
                elif call == 'HAP2':
                    SeqIO.write(record, hap2_out, 'fasta')
                    hap2_diploid += 1
                elif call == 'HAP2_HAPLOID':
                    SeqIO.write(record, hap2_out, 'fasta')
                    hap2_haploid += 1
                elif call == 'UNKNOWN':
                    SeqIO.write(record, unknown_out, 'fasta')
                    unknown += 1
                elif call == 'DROP':
                    SeqIO.write(record, unknown_out, 'fasta')
                    dropped += 1
                elif call == 'UNPHASED':
                    SeqIO.write(record, unknown_out, 'fasta')
                    unphased += 1
            else:
                SeqIO.write(record, unknown_out, 'fasta')
                unobserved += 1

        LOGGER.log_sect('{0} reads in alignment file'.format(total_reads))
        LOGGER.log_sect('{0} reads written'.format(total_written))
        LOGGER.log_sect('{0} reads are assigned to hap1'.format(hap1_diploid + hap1_haploid))
        LOGGER.log_sect('  {0} diploid contigs'.format(hap1_diploid))
        LOGGER.log_sect('  {0} haploid contigs'.format(hap1_haploid))
        LOGGER.log_sect('{0} reads are assigned to hap2'.format(hap2_diploid + hap2_haploid))
        LOGGER.log_sect('  {0} diploid contigs'.format(hap2_diploid))
        LOGGER.log_sect('  {0} haploid contigs'.format(hap2_haploid))
        LOGGER.log_sect("{0} reads are can't be assigned to maternal/paternal with genotyping".format(unknown))
        LOGGER.log_sect("{0} reads are not observed, probably between split phased regions, and can't be assigned".format(unobserved))
        LOGGER.log_sect("{0} reads couldn't be assigned phase by whatshap and were removed.".format(dropped))
        LOGGER.log_sect('{0} reads are unphased'.format(unphased))
        LOGGER.flush_sect('Fasta Stats')


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit(1)

