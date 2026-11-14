#!/usr/bin/env python3
#
# sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
# Copyright (C) 2019-2026 Thiago L. A. Miller <tmiller@mochsl.org.br
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Simulate the alignments of a retrocopy insertion for testing.

The simulation builds a random reference with two contigs: 'chr1',
where the retrocopy is inserted, and 'chr2', with the parental gene.
The inserted sequence is the spliced parental exons - truncated at
the 5' end by --truncation bases - followed by a poly-A tail,
flanked by a target site duplication (TSD):

    chr1:  left flank | TSD |                            | TSD | right flank
    alt:   left flank | TSD | exon1 exon2 exon3 poly-A   | TSD | right flank

No aligner is run. The reads are sampled from the alternative
haplotype and written with the SAM fields an aligner would produce:

  - Discordant pairs: one read at a chr1 flank and its mate at
    a chr2 exon.
  - Left junction split reads: the read crosses the junction
    'TSD | exon1'. The primary alignment is clipped at the right
    ('<m>M<c>S') and the supplementary alignment is at exon1.
  - Right junction reads: the read crosses the junction
    'poly-A | TSD'. The primary alignment is clipped at the left
    ('<c>S<m>M') and there is no supplementary alignment, because
    the poly-A does not map anywhere.
  - Exon-exon junction reads: discordant pairs whose chr2 mate
    crosses the junction 'exon1 | exon2' of the retrocopy. It is
    aligned in two pieces: The primary at the end of exon1,
    clipped at the right, and the supplementary at the beginning
    of exon2.
  - Noise clips (--n-noise-clips): discordant reads at the flanks
    soft clipped at random positions, as produced by adapters or
    quality trimming. They must not define a breakpoint.

The sequences are written at the forward strand of the reference,
as in the SAM format, whatever the strand of the read (flag 0x10).
The output is deterministic for a given seed.

Output files (into --output-dir):

  sim.sam     Alignments grouped by queryname (SO:queryname)
  ref.fa      Reference genome
  genes.gtf   Parental gene annotation
  truth.tsv   Expected values: key<TAB>value

Example:

  $ tests/sim/retrocopy_sim.py --seed 11 --output-dir sim
  $ sider process-sample -s -a sim/genes.gtf -o sim sim/sim.sam
  $ sider merge-call -r sim/ref.fa -o sim -p merged sim/out.db
  $ sider make-vcf -r sim/ref.fa -o sim -p merged sim/merged.db
  $ tests/sim/retrocopy_sim.py --check sim/merged.vcf sim/truth.tsv
"""

import argparse
import os
import random
import sys

READ_LEN = 100
QUAL = 'I'

# Illumina like read names, with colons, as
# 'A00123:8:HSIMFLOWX:1:<qname>'
READ_PREFIX = 'A00123:8:HSIMFLOWX:1:'

# Must mirror RETROCOPY_MIN_CLIP, RETROCOPY_MIN_POLYA and
# RETROCOPY_MIN_SEED at src/retrocopy.h
MIN_CLIP = 10
MIN_POLYA = 6
MIN_SEED = 12

# Each right junction read has this many more exon
# bases at its clip than the previous one
RIGHT_CLIP_STEP = 3


def random_seq(n):
    return ''.join(random.choice('ACGT') for _ in range(n))


def fasta_wrap(seq, width=60):
    return '\n'.join(seq[i:i + width] for i in range(0, len(seq), width))


class Sim:
    def __init__(self, args):
        self.args = args
        self.recs = []

        random.seed(args.seed)

        self.chr1 = random_seq(args.chr1_len)
        self.chr2 = random_seq(args.chr2_len)

        # Parental exons at chr2: 1-based inclusive, + strand
        self.exons = [(1001, 1400), (2001, 2400), (3001, 3400)]

        # Inserted sequence: spliced exons, 5' truncated, + poly-A
        self.spliced = ''.join(self.chr2[s - 1:e] for s, e in self.exons)
        self.ins = self.spliced[args.truncation:] + 'A' * args.polya_len

        # TSD at chr1: 1-based inclusive
        self.tsd_start = args.insertion_point
        self.tsd_end = args.insertion_point + args.tsd_len - 1
        self.tsd = self.chr1[self.tsd_start - 1:self.tsd_end]

        # Alternative haplotype
        left = self.chr1[:self.tsd_end]
        right = self.chr1[self.tsd_start - 1:]
        self.alt = left + self.ins + right

        # Alternative haplotype offsets (0-based) of the junctions
        self.left_junction = len(left)
        self.right_junction = len(left) + len(self.ins)

    def rec(self, qname, flag, chr_, pos, cigar, rnext, pnext, seq, mapq=60):
        self.recs.append('\t'.join(map(str, [
            READ_PREFIX + qname, flag, chr_, pos, mapq, cigar, rnext, pnext, 0,
            seq, QUAL * len(seq)
        ])))

    def exon_pos(self, exon):
        s, e = self.exons[exon]
        return random.randint(s, e - READ_LEN + 1)

    def exon_read(self, exon):
        pos = self.exon_pos(exon)
        return pos, self.chr2[pos - 1:pos - 1 + READ_LEN]

    def discordant(self):
        a = self.args

        for i in range(a.n_discordant):
            qname = 'disc%d' % i

            if i % 2 == 0:
                # Left flank, forward read. Mate at the first exon
                pos = random.randint(self.tsd_start - a.fragment_len,
                                     self.tsd_start - READ_LEN - 20)
                seq = self.chr1[pos - 1:pos - 1 + READ_LEN]
                mpos, mseq = self.exon_read(0)

                self.rec(qname, 97, 'chr1', pos, '%dM' % READ_LEN, 'chr2', mpos, seq)
                self.rec(qname, 145, 'chr2', mpos, '%dM' % READ_LEN, 'chr1', pos,
                         mseq)
            else:
                # Right flank, reverse read. Mate at the last exon
                pos = random.randint(self.tsd_end + 20,
                                     self.tsd_end + a.fragment_len - READ_LEN)
                seq = self.chr1[pos - 1:pos - 1 + READ_LEN]
                mpos, mseq = self.exon_read(len(self.exons) - 1)

                self.rec(qname, 81, 'chr1', pos, '%dM' % READ_LEN, 'chr2', mpos,
                         seq)
                self.rec(qname, 161, 'chr2', mpos, '%dM' % READ_LEN, 'chr1', pos, mseq)

    def junction_reads(self):
        a = self.args
        exon1_len = self.exons[0][1] - self.exons[0][0] + 1

        for i in range(a.n_junction_reads):
            qname = 'junc%d' % i

            # Matched bases at exon1 before the junction
            m = 40 + i
            c = READ_LEN - m

            start = exon1_len - a.truncation - m
            seq = self.ins[start:start + READ_LEN]

            # Left flank read, forward. The mate spans the junction
            pos = random.randint(self.tsd_start - a.fragment_len,
                                 self.tsd_start - READ_LEN - 20)
            fseq = self.chr1[pos - 1:pos - 1 + READ_LEN]
            mpos = self.exons[0][1] - m + 1

            self.rec(qname, 97, 'chr1', pos, '%dM' % READ_LEN, 'chr2', mpos, fseq)
            self.rec(qname, 145, 'chr2', mpos, '%dM%dS' % (m, c), 'chr1', pos, seq)
            self.rec(qname, 2193, 'chr2', self.exons[1][0], '%dH%dM' % (m, c),
                     'chr1', pos, seq[m:])

    def noise_clips(self):
        a = self.args

        for i in range(a.n_noise_clips):
            qname = 'noise%d' % i

            # Left flank, forward read, clipped at the right
            # as if the end of the read was trimmed
            pos = random.randint(self.tsd_start - a.fragment_len,
                                 self.tsd_start - READ_LEN - 20)
            c = random.randint(10, 30)
            m = READ_LEN - c
            seq = self.chr1[pos - 1:pos - 1 + READ_LEN]
            mpos, mseq = self.exon_read(0)

            self.rec(qname, 97, 'chr1', pos, '%dM%dS' % (m, c), 'chr2', mpos, seq)
            self.rec(qname, 145, 'chr2', mpos, '%dM' % READ_LEN, 'chr1', pos,
                     mseq)

    def split_left(self):
        a = self.args

        for i in range(a.n_split_left):
            qname = 'splitL%d' % i

            # Matched bases before the junction
            m = a.min_match + i
            c = READ_LEN - m

            start = self.left_junction - m
            seq = self.alt[start:start + READ_LEN]
            mpos, mseq = self.exon_read(0)

            # Primary clipped at the right and the
            # supplementary at the beginning of exon1
            self.rec(qname, 97, 'chr1', start + 1, '%dM%dS' % (m, c), 'chr2', mpos, seq)
            self.rec(qname, 2145, 'chr2', self.exons[0][0], '%dH%dM' % (m, c),
                     'chr2', mpos, seq[m:])
            self.rec(qname, 145, 'chr2', mpos, '%dM' % READ_LEN, 'chr1', start + 1,
                     mseq)

    def split_right(self):
        a = self.args

        for i in range(a.n_split_right):
            qname = 'splitR%d' % i

            # Clipped bases: the poly-A and some exon bases
            c = a.polya_len + RIGHT_CLIP_STEP * i
            m = READ_LEN - c

            start = self.right_junction - c
            seq = self.alt[start:start + READ_LEN]
            mpos, mseq = self.exon_read(len(self.exons) - 1)

            # Primary clipped at the left. No supplementary
            self.rec(qname, 81, 'chr1', self.tsd_start, '%dS%dM' % (c, m), 'chr2', mpos,
                     seq)
            self.rec(qname, 161, 'chr2', mpos, '%dM' % READ_LEN, 'chr1', self.tsd_start, mseq)

    def truth(self):
        """
        Expected VCF values. 'NA' means the INFO field must be
        absent. The left breakpoint needs one true split read,
        with a supplementary alignment. The right breakpoint
        reads have no supplementary alignment, so at least two
        of them must agree. Without breakpoints the retrocopy
        is IMPRECISE and POS is the window mean (not checked).
        """
        a = self.args

        # Right junction reads clipped less than MIN_CLIP
        # bases are ignored
        clips = [a.polya_len + RIGHT_CLIP_STEP * i for i in range(a.n_split_right)
                 if a.polya_len + RIGHT_CLIP_STEP * i >= MIN_CLIP]

        left = a.n_split_left >= 1
        right = len(clips) >= 2
        both = left and right

        # The poly-A run observed in the reads includes the
        # A bases at the end of the last exon. A read gives
        # the exact length only if its clip is longer than
        # the run, otherwise just a lower bound. Once the
        # breakpoint is known, every read clipped there with
        # a run of at least MIN_POLYA supports the poly-A,
        # even those clipped less than MIN_CLIP
        run = len(self.ins) - len(self.ins.rstrip('A'))
        all_clips = [a.polya_len + RIGHT_CLIP_STEP * i for i in range(a.n_split_right)]
        polya_reads = [c for c in all_clips if min(c, run) >= MIN_POLYA]
        polya = right and len(polya_reads) > 0
        exact = any(c > run for c in polya_reads)

        if not polya:
            polya_len, polya_min, polya_sr = 'NA', 'NA', 'NA'
        elif exact:
            polya_len, polya_min, polya_sr = run, 'NA', len(polya_reads)
        else:
            polya_len, polya_min, polya_sr = max(polya_reads), '', len(polya_reads)

        # Retrocopy ends into the parental gene (needs -r at
        # merge-call). The 5' end comes from the left junction
        # clips, which begin with the insertion: exon1 start +
        # truncation. The 3' end comes from the right junction
        # clips with at least MIN_SEED exon bases before the
        # poly-A run: the last exon end minus its trailing A.
        # The predicted length is the insertion length, but the
        # exon trailing A are counted twice when the 3' end is
        # assumed at the exon end
        exon_trailing = run - a.polya_len
        five_pos = self.exons[0][0] + a.truncation
        three_pos = self.exons[-1][1] - exon_trailing
        three_reads = [c for c in all_clips if c - run >= MIN_SEED]

        pg5p = 'ENSEP0:%d:%d' % (five_pos, a.n_split_left) if left else 'NA'
        pg3p = ('ENSEP%d:%d:%d' % (len(self.exons) - 1, three_pos, len(three_reads))
                if right and three_reads else 'NA')

        # The predicted sequence and its length. 'SKIP' means the
        # value is not checked: Without the poly-A the 3' end is
        # not well defined by this simulation
        if not left or not polya:
            svlen, insseq = ('NA', 'NA') if not left else ('SKIP', 'SKIP')
        elif right and three_reads:
            insseq = self.ins
            svlen = len(insseq)
        else:
            insseq = self.spliced[a.truncation:] + 'A' * run
            svlen = len(insseq)

        if not polya:
            pg3p = 'SKIP'

        # Observed sequences at the junctions: The consensus of
        # the clips, which are prefixes (left) or suffixes (right)
        # of the insertion
        left_seq = self.ins[:READ_LEN - a.min_match] if left else 'NA'
        right_seq = self.ins[-max(all_clips):] if right else 'NA'

        if right:
            pos = self.tsd_start - 1
        elif left:
            pos = self.tsd_end
        else:
            pos = 'NA'

        return [
            ('CHROM', 'chr1'),
            ('ID', 'RTC000001'),
            ('POS', pos),
            ('IMPRECISE', 'NA' if left or right else ''),
            ('LBKP', self.tsd_end if left else 'NA'),
            ('LSR', a.n_split_left if left else 'NA'),
            ('RBKP', self.tsd_start if right else 'NA'),
            ('RSR', len(clips) if right else 'NA'),
            ('TSDLEN', a.tsd_len if both else 'NA'),
            ('TSD', self.tsd if both else 'NA'),
            ('POLYA', polya_len),
            ('POLYAMIN', polya_min),
            ('POLYASR', polya_sr),
            ('POLARITY', '+'),
            ('PG', 'PARENT'),
            ('PGEXONS', 'PARENT:%d' % len(self.exons)),
            ('PG5P', pg5p),
            ('PG3P', pg3p),
            ('SVLEN', svlen),
            ('LEFT_SVINSSEQ', left_seq),
            ('RIGHT_SVINSSEQ', right_seq),
            ('INSSEQ', insseq),
            ('RNAMES', ','.join(sorted(set(r.split('\t')[0] for r in self.recs)))),
            ('ROLES', self.roles(left, right, run, all_clips, polya)),
            ('FILTER', self.filters(both, polya)),
            ('HALLMARKS', ','.join(h for h, ok in (
                ('TSD', both and a.tsd_len > 0),
                ('POLYA', polya),
                ('EEJ', a.n_junction_reads > 0)) if ok) or 'NA'),
            ('EEJ', a.n_junction_reads if a.n_junction_reads else 'NA'),
            ('EEJLIST', 'ENSEP0-ENSEP1' if a.n_junction_reads else 'NA'),
            ('INS', self.ins),
        ]

    def filters(self, both, polya):
        """
        Expected FILTER: The simulated retrocopy has a single,
        multi-exon parental gene, reads with MAPQ 60, the poly-A
        at '+' and a 5' end at the first exon. A check that cannot
        be evaluated (no breakpoints) does not fail
        """
        a = self.args
        f = []

        if both and a.tsd_len <= 0:
            f.append('NoTSD')
        if both and (a.tsd_len > 50 or a.tsd_len < -20):
            f.append('TSDOutlier')
        if both and not polya:
            f.append('NoPolyA')
        if a.n_junction_reads == 0:
            f.append('NoJunction')

        return ';'.join(f) if f else 'PASS'

    def roles(self, left, right, run, all_clips, polya):
        """
        Expected role of each read: 'name=role1,role2' separated
        by ';'. Only the reads supporting an accepted breakpoint
        (clipped at least MIN_CLIP bases) are split reads; a right
        junction read supports the poly-A if its run is long enough
        """
        a = self.args
        roles = {}

        for rec in self.recs:
            name = rec.split('\t')[0]
            kind = name[len(READ_PREFIX):]
            r = []

            if kind.startswith('splitL') and left:
                r.append('split_left')
            elif kind.startswith('splitR') and right:
                c = all_clips[int(kind[len('splitR'):])]
                if c >= MIN_CLIP:
                    r.append('split_right')
                if polya and min(c, run) >= MIN_POLYA:
                    r.append('polya')
            elif kind.startswith('junc'):
                r.append('junction')

            roles[name] = ','.join(r) if r else 'discordant'

        return ';'.join('%s=%s' % kv for kv in sorted(roles.items()))

    def write(self):
        a = self.args
        os.makedirs(a.output_dir, exist_ok=True)

        self.discordant()
        self.junction_reads()
        self.noise_clips()
        self.split_left()
        self.split_right()

        with open(os.path.join(a.output_dir, 'sim.sam'), 'w') as fp:
            fp.write('@HD\tVN:1.0\tSO:queryname\n')
            fp.write('@SQ\tSN:chr1\tLN:%d\n' % len(self.chr1))
            fp.write('@SQ\tSN:chr2\tLN:%d\n' % len(self.chr2))
            fp.write('\n'.join(self.recs) + '\n')

        with open(os.path.join(a.output_dir, 'ref.fa'), 'w') as fp:
            fp.write('>chr1\n%s\n' % fasta_wrap(self.chr1))
            fp.write('>chr2\n%s\n' % fasta_wrap(self.chr2))

        with open(os.path.join(a.output_dir, 'genes.gtf'), 'w') as fp:
            for k, (s, e) in enumerate(self.exons):
                fp.write('chr2\t.\texon\t%d\t%d\t.\t+\t.\t'
                         'gene_name "PARENT"; gene_id "ENSGP"; transcript_id "TP"; '
                         'transcript_type "protein_coding"; exon_id "ENSEP%d";\n'
                         % (s, e, k))

        with open(os.path.join(a.output_dir, 'truth.tsv'), 'w') as fp:
            for key, value in self.truth():
                fp.write('%s\t%s\n' % (key, value))


def check(vcf_file, truth_file, no_ref=False):
    """
    Compare the VCF body with truth.tsv. Return the number of
    errors. With no_ref, merge-call ran without the reference
    genome, so the retrocopy ends and the predicted sequence
    must be absent
    """
    truth = {}
    with open(truth_file) as fp:
        for line in fp:
            key, value = line.rstrip('\n').split('\t')
            truth[key] = value

    if no_ref:
        for key in ('PG5P', 'PG3P', 'SVLEN'):
            truth[key] = 'NA'
        truth['INSSEQ'] = 'NA'

    lines = [l for l in open(vcf_file) if not l.startswith('#')]
    if len(lines) != 1:
        print('FAIL: expected 1 retrocopy, found %d' % len(lines))
        return 1

    fields = lines[0].rstrip('\n').split('\t')
    info = dict(kv.split('=', 1) if '=' in kv else (kv, '')
                for kv in fields[7].split(';'))

    errors = 0

    def expect(what, got, exp):
        nonlocal errors
        if str(got) != str(exp):
            print('FAIL: %s = %s, expected %s' % (what, got, exp))
            errors += 1

    expect('CHROM', fields[0], truth['CHROM'])
    expect('ID', fields[2], truth['ID'])
    expect('FILTER', fields[6], truth['FILTER'])

    # The single sample carries the insertion: High QUAL
    if fields[5] == '.' or float(fields[5]) < 20:
        print('FAIL: QUAL = %s, expected >= 20' % fields[5])
        errors += 1

    if truth['POS'] != 'NA':
        expect('POS', fields[1], truth['POS'])

    for key in ('IMPRECISE', 'LBKP', 'LSR', 'RBKP', 'RSR', 'TSDLEN', 'TSD',
                'POLYA', 'POLYAMIN', 'POLYASR', 'POLARITY', 'PG',
                'PGEXONS', 'EEJ', 'EEJLIST', 'PG5P', 'PG3P', 'SVLEN',
                'LEFT_SVINSSEQ', 'RIGHT_SVINSSEQ', 'HALLMARKS'):
        if truth[key] == 'SKIP':
            continue
        if truth[key] == 'NA':
            if key in info:
                print('FAIL: %s = %s, expected absent' % (key, info[key]))
                errors += 1
        else:
            expect(key, info.get(key, '<missing>'), truth[key])

    # Side files written by make-vcf next to the VCF
    base = vcf_file[:-4] if vcf_file.endswith('.vcf') else vcf_file
    rnames_file = base + '.rnames.tsv'
    ins_file = base + '.insertions.fa'

    if os.path.exists(rnames_file):
        lines = open(rnames_file).read().splitlines()
        rows = [l.split('\t') for l in lines[1:]]
        if not lines or lines[0] != '#retrocopy_id\tsample\tread_name\trole':
            print('FAIL: rnames file header: %s' % (lines[0] if lines else '<empty>'))
            errors += 1
        elif any(len(r) != 4 or r[0] != truth['ID'] or r[1] != 'sim' for r in rows):
            print('FAIL: rnames rows must be <%s, sim, read_name, role>' % truth['ID'])
            errors += 1
        else:
            names = [r[2] for r in rows]
            if names != sorted(names):
                print('FAIL: rnames read names are not sorted')
                errors += 1
            missing = set(truth['RNAMES'].split(',')) - set(names)
            if missing:
                print('FAIL: rnames missing reads: %s' % ','.join(sorted(missing)))
                errors += 1

            # Roles: Without the reference, nothing changes
            expected = dict(kv.split('=') for kv in truth['ROLES'].split(';'))
            for name, role in ((r[2], r[3]) for r in rows):
                if expected.get(name) != role:
                    print('FAIL: rnames role of %s = %s, expected %s'
                          % (name, role, expected.get(name)))
                    errors += 1

    if os.path.exists(ins_file) and truth['INSSEQ'] != 'SKIP':
        lines = open(ins_file).read().splitlines()
        seq = ''.join(l for l in lines if not l.startswith('>'))
        if truth['INSSEQ'] == 'NA':
            if lines:
                print('FAIL: insertions file must be empty')
                errors += 1
        else:
            if not lines or not lines[0].startswith('>' + truth['ID'] + ' '):
                print('FAIL: insertions file must start with >%s' % truth['ID'])
                errors += 1
            expect('INSSEQ', seq, truth['INSSEQ'])

    if not errors:
        print('OK: %s matches %s' % (vcf_file, truth_file))

    return errors


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--check', nargs=2, metavar=('VCF', 'TRUTH'),
                   help='compare a VCF with a truth.tsv instead of simulating')
    p.add_argument('--no-ref', action='store_true',
                   help='with --check: merge-call ran without --reference-file')
    p.add_argument('-o', '--output-dir', default='.')
    p.add_argument('-s', '--seed', type=int, default=11)
    p.add_argument('--chr1-len', type=int, default=20000)
    p.add_argument('--chr2-len', type=int, default=10000)
    p.add_argument('--insertion-point', type=int, default=10000,
                   help='first base of the TSD at chr1 [10000]')
    p.add_argument('--tsd-len', type=int, default=12)
    p.add_argument('--polya-len', type=int, default=25)
    p.add_argument('--truncation', type=int, default=0,
                   help="bases missing at the 5' end of the retrocopy [0]")
    p.add_argument('--fragment-len', type=int, default=700)
    p.add_argument('--min-match', type=int, default=60,
                   help='matched bases of the first left split read [60]')
    p.add_argument('--n-discordant', type=int, default=20)
    p.add_argument('--n-split-left', type=int, default=8)
    p.add_argument('--n-split-right', type=int, default=8)
    p.add_argument('--n-junction-reads', type=int, default=6,
                   help='reads spanning the exon1-exon2 junction [6]')
    p.add_argument('--n-noise-clips', type=int, default=0,
                   help='discordant reads with random soft clips [0]')
    args = p.parse_args()

    if args.check:
        return check(*args.check, no_ref=args.no_ref)

    Sim(args).write()
    return 0


if __name__ == '__main__':
    sys.exit(main())
