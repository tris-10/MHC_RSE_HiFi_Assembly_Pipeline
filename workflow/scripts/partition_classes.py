"""
Haplotype Partioning classes
"""

from collections import defaultdict
import sys


class ProbeInfo:
    '''
    This class stores information about the heterozygous SNP-chip probes expected to be captured
    in the assembly.
    '''

    def __init__(self, hg19_pos, base1, base2):
        self.hg19_pos = hg19_pos
        self.base1 = base1
        self.base2 = base2
        self.spike = True if hg19_pos == 'X' else False


class ProbeAlign:
    '''
    Class handles information on where a SNP-CHIP probe aligns in the assembly.
    '''

    rev_comp = {"A": "T", "T": "A", "C": "G", "G": "C", "X": "X"}

    def __init__(self, snp_id, contig_pos, contig_name, is_reverse, probe_info):
        self.snp_id = snp_id
        self.contig_pos = contig_pos
        self.contig_name = contig_name
        self.key = '{0}:{1}'.format(contig_name, contig_pos)
        self.is_reverse = is_reverse
        self.probe_info = probe_info
        self.haplotype = None
        self.edge = False

    def set_edge(self):
        self.edge = True

    def set_haplotype(self, hap):
        self.haplotype = hap

    def get_base(self, is_base1):
        base = self.probe_info.base1 if is_base1 else self.probe_info.base2
        return ProbeAlign.rev_comp[base] if self.is_reverse else base


class PhaseBlock:
    '''
    Information about phased/unphased blocks
    '''

    prefix = ["DRB1", "DRB3", "DRB4", "DRB5", "DPA1", "DPB1", "DQA1", "DQB1", "A", "B", "C"]

    def __init__(self, phase_set_name, contig, start, end):
        self.contig = contig
        self.name = phase_set_name
        self.start = int(start)
        self.end = int(end)
        self.length = self.end - self.start
        self.snp_list = []
        self.vcf_pos = []
        self.genotype_string = None
        self.homozygous = False
        self.edge_reduced = False
        self.hap1_count = 0
        self.hap2_count = 0
        self.hap_string = 0

    def add_snp(self, snp):
        self.snp_list.append(snp)

    def add_vcf_pos(self, vcf_pos):
        self.vcf_pos.append(vcf_pos)

    def add_snps(self, snps):
        self.snp_list = snps

    def print_snps(self):
        for snp in self.snp_list:
            print(snp.rs_number, snp.cell_line)

    @staticmethod
    def generate_hap_percent(self, count1, count2):
        return count1 / (count1 + count2) if count1 > 0 else 0

    def classify_block(self, min_count, max_contam_count, max_contam_frac, unphased_snps=None):
        hap1_cnt = hap2_cnt = ehap1_cnt = ehap2_cnt = unphased_cnt = switches = 0
        hap_string = ""
        ehap_string = ""

        # If spike in, only use the most common locus, to avoid DRB misalignments
        loc_count = defaultdict(int)
        for snp in self.snp_list:
            if snp.probe_info.spike:
                locus, pos = snp.snp_id.split(':')
                loc_count[locus] += 1
        if len(loc_count.keys()) > 0:
            top_loc = max(loc_count, key=loc_count.get)
        else:
            top_loc = None

        for snp in self.snp_list:
            # Skip, if SNP doesn't come from the most common locus
            if snp.probe_info.spike:
                locus, pos = snp.snp_id.split(':')
                if locus != top_loc:
                    continue

            if unphased_snps is not None and self.contig in unphased_snps and snp.contig_pos in unphased_snps[self.contig]:
                unphased_cnt += 1
            elif snp.haplotype == 'HAP2':
                hap2_cnt += 1
                hap_string += '2'
                if not snp.edge:
                    ehap2_cnt += 1
                    ehap_string += '2'
            elif snp.haplotype == 'HAP1':
                hap1_cnt += 1
                hap_string += '1'
                if not snp.edge:
                    ehap1_cnt += 1
                    ehap_string += '1'

        hap1p = hap1_cnt / (hap1_cnt + hap2_cnt) if hap1_cnt > 0 else 0
        hap2p = hap2_cnt / (hap1_cnt + hap2_cnt) if hap2_cnt > 0 else 0
        ehap1p = ehap1_cnt / (ehap1_cnt + ehap2_cnt) if ehap1_cnt > 0 else 0
        ehap2p = ehap2_cnt / (ehap1_cnt + ehap2_cnt) if ehap2_cnt > 0 else 0

        #If we clean up the haplotype call by removing edge variants, do it
        if (hap1_cnt > hap2_cnt and ehap1p > hap1p) or (hap2_cnt > hap1_cnt and ehap2p > hap2p):
            hap1_cnt, hap2_cnt, hap1p, hap2p, h_string = ehap1_cnt, ehap2_cnt, ehap1p, ehap2p, ehap_string
            self.edge_reduced = True

        self.hap1_count = hap1_cnt
        self.hap2_count = hap2_cnt
        self.hap_string = hap_string

        # check for haplotype switches
        last = None
        for s in hap_string:
            if last is not None and last != s:
                switches += 1
            last = s

        if len(self.snp_list) == 1 and self.snp_list[0].edge:
            return 'UNKNOWN'
        if unphased_cnt > 0:
            return 'UNPHASED'
        elif hap1_cnt > hap2_cnt and hap1_cnt >= min_count:
            if switches == 0:
                return 'HAP1_CLEAN'
            elif hap1p >= (1-max_contam_frac) and hap2_cnt < max_contam_count:
                return 'HAP1_CONTAM'
            else:
                return 'MIXED'
        elif hap2_cnt > hap1_cnt and hap2_cnt >= min_count:
            if switches == 0:
                return 'HAP2_CLEAN'
            elif hap2p >= (1-max_contam_frac) and hap1_cnt < max_contam_count:
                return 'HAP2_CONTAM'
            else:
                return 'MIXED'
        else:
            return 'UNKNOWN'

    def split_block(self, limit):
        new_blocks = []

        start_pos = last_pos = self.start
        curr_snps = []
        last_hap = None
        for snp in self.snp_list:
            if last_hap is not None and last_hap != snp.haplotype:
                b = PhaseBlock(self.contig + ':' + str(start_pos) + '-' + str(last_pos) + '_gen_split', self.contig,
                               start_pos, last_pos)
                b.add_snps(curr_snps)
                new_blocks.append(b)
                start_pos = snp.contig_pos
                curr_snps = []

            last_hap = snp.haplotype
            last_pos = snp.contig_pos
            curr_snps.append(snp)

        b = PhaseBlock(self.contig + ':' + str(start_pos) + '-' + str(self.end) + '_gen_split', self.contig, start_pos,
                       self.end)
        b.add_snps(curr_snps)
        new_blocks.append(b)

        stranded = 0
        for pos in self.vcf_pos:
            for block in new_blocks:
                if block.start <= pos <= block.end:
                    block.vcf_pos.append(pos)
                    break
            else:
                stranded += 1
        print("Splitting stranded {0} VCF variants".format(stranded))

        if len(new_blocks) > limit:
            return None
        else:
            return new_blocks


class InfoLogger:
    """
    Write information to stdout or file
    """
    def __init__(self, output, silent):
        self.silent = silent
        self.handle = output if output == sys.stdout else open(output, 'w')
        self.lines = []

    def log_sect(self, line):
        self.lines.append(line)

    def log(self, line):
        if not self.silent:
            self.handle.write('{0}\n'.format(line))

    def flush_sect(self, section_header):
        if not self.silent:
            banner = ''.join(['*'] * 100)
            self.handle.write('{0}\n{1}\n{2}\n{3}\n'.format(banner, section_header, '\n'.join(self.lines), banner))
        self.lines.clear()

    def close(self):
        if self.handle != sys.stdout:
            self.handle.close()

