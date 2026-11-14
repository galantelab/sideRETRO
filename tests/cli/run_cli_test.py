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
Test the command line of sider: Help, version, citation, invalid
options and full runs with the options. Used by 'meson test':

  $ run_cli_test.py SIDER RETROCOPY_SIM

Each case runs sider with some arguments and checks the exit code
and that some texts are present at its output (stdout and stderr).
"""

import os
import shutil
import subprocess
import sys
import tempfile

OK = 0
FAIL = 1


class Runner:
    def __init__(self, sider, d):
        self.sider = sider
        self.d = d
        self.failed = []
        self.n = 0

    def path(self, name):
        return os.path.join(self.d, name)

    def run(self, name, args, rc=OK, expect=(), absent=(), env=None):
        self.n += 1
        e = dict(os.environ)
        e.update(env or {})
        p = subprocess.run([self.sider] + args, capture_output=True,
                           text=True, cwd=self.d, env=e)
        out = p.stdout + p.stderr
        errors = []

        if p.returncode != rc:
            errors.append('exit code %d, expected %d' % (p.returncode, rc))
        for text in expect:
            if text not in out:
                errors.append('missing %r' % text)
        for text in absent:
            if text in out:
                errors.append('unexpected %r' % text)

        print('%-50s %s' % (name, 'ok' if not errors else 'FAIL'))
        if errors:
            self.failed.append(name)
            print('    $ sider ' + ' '.join(args))
            for err in errors:
                print('    ' + err)
            print('    ' + out[-600:].replace('\n', '\n    '))

        return out


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 2

    sider = os.path.abspath(sys.argv[1])
    sim = os.path.abspath(sys.argv[2])
    d = tempfile.mkdtemp(prefix='sider_cli_')
    t = Runner(sider, d)

    try:
        cases(t, sim)
    finally:
        shutil.rmtree(d)

    if t.failed:
        print('\n%d of %d cases failed' % (len(t.failed), t.n))
        return 1

    print('\nAll %d cases passed' % t.n)
    return 0


def cases(t, sim):
    P = t.path

    # Data: a simulated retrocopy
    subprocess.run([sim, '-o', t.d], check=True, capture_output=True)
    with open(P('list.txt'), 'w') as fp:
        fp.write(P('sim.sam') + '\n\n')
    with open(P('black.bed'), 'w') as fp:
        fp.write('chr1\t0\t100\tblack\n')
    with open(P('genes.gtf')) as gtf, open(P('genes.bed'), 'w') as bed:
        for line in gtf:
            f = line.rstrip('\n').split('\t')
            gene = f[8].split('gene_name "')[1].split('"')[0]
            bed.write('%s\t%d\t%s\t%s\t0\t%s\n'
                      % (f[0], int(f[3]) - 1, f[4], gene, f[6]))
    with open(P('black.gtf'), 'w') as fp:
        fp.write('chr1\t.\tgene\t1\t100\t.\t+\t.\tgene_name "black"; '
                 'gene_type "processed_pseudogene";\n')

    # ================= sider
    t.run('sider: usage', [], expect=['Usage: sider [-hv]', 'Commands:'])
    t.run('sider: help', ['--help'], expect=['Commands:'])
    t.run('sider: version', ['-v'], expect=['sideRETRO'])
    t.run('sider: citation', ['--cite'], expect=['@article', 'btaa689'])
    t.run('sider: invalid option', ['--ponga'], rc=FAIL, expect=['--help'])
    t.run('sider: no option', ['--'], expect=['Commands:'])
    t.run('sider: invalid command', ['ponga'], rc=FAIL,
          expect=["'ponga' is not a valid command"])

    # ================= process-sample
    ps = ['process-sample']
    ann = ['-a', P('genes.gtf')]
    sam = [P('sim.sam')]

    t.run('ps: usage', ps, expect=['Usage: sider process-sample'])
    t.run('ps: help', ps + ['-h'], expect=['Annotation Options:'])
    t.run('ps: invalid option', ps + ['--ponga'], rc=FAIL)
    t.run('ps: no alignment file', ps + ann, rc=FAIL,
          expect=['Missing alignment files'])
    t.run('ps: alignment file not found', ps + ann + ['ponga.bam'], rc=FAIL,
          expect=["alignment file 'ponga.bam': No such file"])
    t.run('ps: no annotation', ps + sam, rc=FAIL,
          expect=['Missing annotation file'])
    t.run('ps: annotation not found', ps + ['-a', 'ponga.gtf'] + sam, rc=FAIL,
          expect=["annotation file 'ponga.gtf': No such file"])
    t.run('ps: reciprocal with -F', ps + ann + ['-r', '-F', '0.5'] + sam,
          rc=FAIL, expect=['--reciprocal must be set solely with'])
    t.run('ps: exon-frac', ps + ann + ['-f', '0'] + sam, rc=FAIL,
          expect=['--exon-frac must be in the interval'])
    t.run('ps: alignment-frac', ps + ann + ['-F', '2'] + sam, rc=FAIL,
          expect=['--alignment-frac must be in the interval'])
    t.run('ps: threads', ps + ann + ['-t', '0'] + sam, rc=FAIL,
          expect=['--threads must be greater or equal to 1'])
    t.run('ps: cache-size', ps + ann + ['-c', '1'] + sam, rc=FAIL,
          expect=['--cache-size must be greater or equal to'])
    t.run('ps: max-distance', ps + ann + ['-m', '-1'] + sam, rc=FAIL,
          expect=['--max-distance must be'])
    t.run('ps: max-base-freq', ps + ann + ['-M', '0.1'] + sam, rc=FAIL,
          expect=['--max-base-freq must be in the range'])
    t.run('ps: phred-quality', ps + ann + ['-Q', '-1'] + sam, rc=FAIL,
          expect=['--phred-quality must be'])
    t.run('ps: gff-attribute', ps + ann + ['-A', 'ponga'] + sam, rc=FAIL,
          expect=["--gff-attribute 'ponga'"])

    t.run('ps: run with options',
          ps + ['--annotation-file', P('genes.gtf'), '-i', P('list.txt'),
                '-o', P('ps'), '-p', 'sim', '-l', P('ps.log'), '-d', '-s',
                '-D', '-t', '2', '-c', '300000', '-Q', '1', '-M', '0.9',
                '-m', '9000', '-f', '0.1', '-e', '-T', 'exon',
                '-A', 'transcript_type=protein_coding', '-N', 'gene_name',
                '-E', 'exon_id', P('sim.sam')],
          expect=['process-sample summary', 'After deduplication',
                  '--gff-attribute=', '--log-file=', '--deduplicate',
                  '--sorted', '--either'])
    t.run('ps: run with -r', ps + ann + ['-o', P('ps2'), '-f', '0.5', '-r'] + sam,
          expect=['--exon-frac=5.000000e-01', '--alignment-frac=5.000000e-01'])
    t.run('ps: run quiet with log', ps + ann + ['-o', P('ps3'), '-q',
          '-l', P('ps3.log')] + sam, absent=['summary'])
    t.run('ps: run BED annotation', ps + ['-a', P('genes.bed'), '-o', P('ps4')] + sam,
          expect=['Exons indexed            3'], absent=['--gff-feature'])

    # ================= merge-call
    mc = ['merge-call']
    db = [P('ps/sim.db')]

    t.run('mc: usage', mc, expect=['Usage: sider merge-call'])
    t.run('mc: help', mc + ['-h'], expect=['Blacklist Options:'])
    t.run('mc: invalid option', mc + ['--ponga'], rc=FAIL)
    t.run('mc: no database', mc + ['-o', P('x')], rc=FAIL,
          expect=['Missing SQLite3 databases'])
    t.run('mc: database not found', mc + ['ponga.db'], rc=FAIL,
          expect=["SQLite3 database 'ponga.db': No such file"])
    t.run('mc: blacklist not found', mc + ['-B', 'ponga.bed'] + db, rc=FAIL,
          expect=["--blacklist-region 'ponga.bed': No such file"])
    t.run('mc: cache-size', mc + ['-c', '1'] + db, rc=FAIL,
          expect=['--cache-size must be'])
    t.run('mc: epsilon', mc + ['-e', '-1'] + db, rc=FAIL,
          expect=['--epsilon must be'])
    t.run('mc: min-pts', mc + ['-m', '2'] + db, rc=FAIL,
          expect=['--min-pts must be greater than 2'])
    t.run('mc: parental-distance', mc + ['-x', '-1'] + db, rc=FAIL,
          expect=['--parental-distance must be'])
    t.run('mc: genotype-support', mc + ['-g', '-1'] + db, rc=FAIL,
          expect=['--genotype-support must be'])
    t.run('mc: blacklist-padding', mc + ['-P', '-1'] + db, rc=FAIL,
          expect=['--blacklist-padding must be'])
    t.run('mc: near-gene-rank', mc + ['-n', '0'] + db, rc=FAIL,
          expect=['--near-gene-rank must be'])
    t.run('mc: reference not found', mc + ['-r', 'ponga.fa'] + db, rc=FAIL,
          expect=["--reference-file 'ponga.fa': No such file"])
    t.run('mc: phred-quality', mc + ['-Q', '-1'] + db, rc=FAIL,
          expect=['--phred-quality must be'])
    t.run('mc: threads', mc + ['-t', '0'] + db, rc=FAIL,
          expect=['--threads must be'])
    t.run('mc: gff-attribute', mc + ['-A', '=ponga'] + db, rc=FAIL,
          expect=["--gff-attribute '=ponga'"])

    t.run('mc: run with GFF blacklist and options',
          mc + ['-o', P('mc'), '-p', 'merged', '-l', P('mc.log'), '-d',
                '-c', '300000', '-e', '300', '-m', '3', '-b', 'chrM',
                '-b', 'chrY', '-B', P('black.gtf'), '-P', '10', '-T', 'gene',
                '-A', 'gene_type=processed_pseudogene', '-N', 'gene_name',
                '-x', '100000', '-g', '1', '-n', '3', '-r', P('ref.fa'),
                '-t', '2', '-Q', '1'] + db,
          expect=['merge-call summary', '--blacklist-region=',
                  '--gff-attribute=', '--reference-file=', 'Predicted sequence'])
    t.run('mc: run with BED blacklist, two databases',
          mc + ['-o', P('mc2'), '-B', P('black.bed'), P('ps/sim.db'),
                P('ps2/out.db')],
          expect=['Databases                2'])
    shutil.copy(P('ps/sim.db'), P('inplace.db'))
    with open(P('dbs.txt'), 'w') as fp:
        fp.write(P('inplace.db') + '\n')
    t.run('mc: run in place from a list, quiet',
          mc + ['-I', '-q', '-l', P('mc3.log'), '-i', P('dbs.txt')],
          absent=['summary'])
    shutil.copy(P('ps/sim.db'), P('inplace1.db'))
    t.run('mc: run in place, one database',
          mc + ['-I', P('inplace1.db')],
          expect=['Databases                1'])
    shutil.copy(P('ps/sim.db'), P('inplace2.db'))
    t.run('mc: run in place, two databases',
          mc + ['-I', P('inplace2.db'), P('ps2/out.db')],
          expect=['Databases                2'])
    t.run('mc: no cluster', mc + ['-o', P('mc4'), '-m', '1000'] + db,
          expect=['No cluster found'])

    # ================= make-vcf
    vcf = ['make-vcf']
    mdb = [P('mc/merged.db')]

    t.run('vcf: usage', vcf, expect=['Usage: sider make-vcf'])
    t.run('vcf: help', vcf + ['-h'], expect=['Filter Options:'])
    t.run('vcf: invalid option', vcf + ['--ponga'], rc=FAIL)
    t.run('vcf: database not found', vcf + ['ponga.db'], rc=FAIL,
          expect=["SQLite3 database file 'ponga.db': No such file"])
    t.run('vcf: reference not found', vcf + ['-r', 'ponga.fa'] + mdb, rc=FAIL,
          expect=["--reference-file 'ponga.fa': No such file"])
    t.run('vcf: near-gene-dist', vcf + ['-n', '-1'] + mdb, rc=FAIL,
          expect=['--near-gene-dist must be'])
    t.run('vcf: orientation-error', vcf + ['-e', '-1'] + mdb, rc=FAIL,
          expect=['--orientation-error must be'])
    t.run('vcf: min-qual', vcf + ['-L', '-1'] + mdb, rc=FAIL,
          expect=['--min-qual must be'])
    t.run('vcf: min-mapq', vcf + ['-Q', '-1'] + mdb, rc=FAIL,
          expect=['--min-mapq must be'])
    t.run('vcf: unknown filter', vcf + ['-x', 'Ponga'] + mdb, rc=FAIL,
          expect=["--no-filter 'Ponga': Unknown filter"])

    t.run('vcf: run with options',
          vcf + ['-o', P('vcf'), '-p', 'sim', '-l', P('vcf.log'), '-d',
                 '-r', P('ref.fa'), '-n', '5000', '-e', '0.1', '-L', '10',
                 '-Q', '10', '-x', 'NoJunction', '-x', 'LowQual'] + mdb,
          expect=['make-vcf summary', 'PASS', '--no-filter=NoJunction',
                  '--reference-file='])
    t.run('vcf: run without reference, quiet',
          vcf + ['-o', P('vcf2'), '-q', '-l', P('vcf2.log')] + mdb,
          absent=['summary'])
    t.run('vcf: run without reference',
          vcf + ['-o', P('vcf3')] + mdb,
          expect=["REF is 'N'"])

    with open(P('normal.sam'), 'w') as fp:
        fp.write('@SQ\tSN:chr1\tLN:20000\n'
                 'N1\t99\tchr1\t1\t60\t10M\t=\t20\t29\tATCGATCGAT\t~~~~~~~~~~\n'
                 'N1\t147\tchr1\t20\t60\t10M\t=\t1\t-29\tAAAGGGCCCT\t~~~~~~~~~~\n')
    t.run('ps: two files', ps + ann + ['-o', P('ps6'), P('sim.sam'), P('normal.sam')],
          expect=['Total', 'No abnormal alignment found out of 2'])
    t.run('ps: no abnormal alignment', ps + ann + ['-o', P('ps7'), P('normal.sam')],
          expect=['No abnormal alignment was found: No retrocopy can be called'])
    with open(P('intergenic.sam'), 'w') as fp:
        fp.write('@SQ\tSN:chr1\tLN:20000\n@SQ\tSN:chr2\tLN:10000\n'
                 'X1\t97\tchr1\t500\t60\t10M\tchr2\t9000\t0\tATCGATCGAT\t~~~~~~~~~~\n'
                 'X1\t145\tchr2\t9000\t60\t10M\tchr1\t500\t0\tAAAGGGCCCT\t~~~~~~~~~~\n')
    t.run('ps: abnormal, but no exon', ps + ann + ['-o', P('ps8'), P('intergenic.sam')],
          expect=['Stored                   0 alignments',
                  'No abnormal fragment overlaps an exon: No retrocopy can be called'])

    # ================= quiet without log: --silent
    t.run('ps: silent', ps + ann + ['-o', P('ps5'), '--silent'] + sam,
          absent=['INFO'])


if __name__ == '__main__':
    sys.exit(main())
