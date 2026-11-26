#!/usr/bin/env python3
#
# sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
# Copyright (C) 2019-2026 Thiago L. A. Miller <tmiller@mochsl.org.br>
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
Run the whole pipeline over simulated retrocopies and compare the
results with the simulated truth. Used by 'meson test':

  $ run_sim_test.py SIDER RETROCOPY_SIM

Each scenario runs retrocopy_sim.py, process-sample, merge-call
and make-vcf into a temporary directory, and checks the VCF, the read names and the predicted
insertions files.

A last scenario checks that the calls depend only on the reads: The
same database, with its alignment ids shuffled - as the merge of
other databases does - must give the same results.
"""

import hashlib
import os
import random
import shutil
import sqlite3
import subprocess
import sys
import tempfile

# Scenarios: name => simulator arguments
SCENARIOS = [
    ('default',          []),
    ('truncated',        ['--seed', '2', '--truncation', '150']),
    ('long_polya_noise', ['--seed', '3', '--polya-len', '40',
                          '--n-split-right', '12', '--n-noise-clips', '5']),
    ('short_polya',      ['--seed', '9', '--polya-len', '7']),
    ('no_polya',         ['--seed', '9', '--polya-len', '3', '--n-split-right', '12']),
    ('no_split_reads',   ['--seed', '5', '--n-noise-clips', '8',
                          '--n-split-left', '0', '--n-split-right', '0']),
    ('left_only',        ['--seed', '5', '--n-split-left', '1', '--n-split-right', '0']),
    ('right_single',     ['--seed', '7', '--n-split-left', '0', '--n-split-right', '1']),
    ('right_only',       ['--seed', '7', '--n-split-left', '0', '--n-split-right', '2']),
    ('no_junction',      ['--seed', '4', '--n-junction-reads', '0']),
    ('short_tsd',        ['--seed', '6', '--tsd-len', '4', '--truncation', '37']),
    # Target site: A blunt insertion and a deletion pass as well
    ('blunt',            ['--seed', '6', '--tsd-len', '0']),
    ('target_deletion',  ['--seed', '6', '--tsd-len', '-15']),
    # A deletion longer than the cluster window: The retrocopies
    # at its two flanks are merged into one
    ('long_deletion',    ['--seed', '6', '--tsd-len', '-10000',
                          '--chr1-len', '40000', '--insertion-point', '25000']),
    # Filters: A TSD longer than the maximum (TSDOutlier), and many
    # split reads at just one side (not LowSR, but OneBreakpoint)
    ('long_tsd',         ['--seed', '8', '--tsd-len', '30']),
    ('left_many',        ['--seed', '7', '--n-split-left', '6', '--n-split-right', '0']),
    # No hallmark: A blunt insertion without the poly-A and the
    # exon-exon junctions (NoHallmark), and the same with a TSD
    ('no_hallmark',      ['--seed', '6', '--tsd-len', '0', '--polya-len', '3',
                          '--n-junction-reads', '0']),
    ('tsd_hallmark',     ['--seed', '6', '--polya-len', '3',
                          '--n-junction-reads', '0']),
    # A tandem repeat at the target site (TSITETR): The TSD is not
    # a hallmark there, but the poly-A still is
    ('repeat_site',      ['--seed', '6', '--site-repeat', 'CA']),
    ('repeat_site_only', ['--seed', '6', '--site-repeat', 'CA', '--polya-len', '3',
                          '--n-junction-reads', '0']),
    # A short piece of an exon at the site: No mate inside the
    # parental exons, only the pieces of the split reads
    ('short_insert',     ['--seed', '5', '--n-discordant', '0', '--n-split-right', '0',
                          '--split-mates-flank']),
    # Intron retention: Mates inside an intron of the parental
    # gene. A single one is not enough to report the intron
    ('intron_retention', ['--seed', '4', '--n-intron-mates', '4']),
    ('intron_single',    ['--seed', '4', '--n-intron-mates', '1']),
    # Genotype: Proper pairs clipped at the junctions are from the
    # alternate allele, with both or just one breakpoint known
    ('proper_clips',     ['--seed', '6', '--n-proper-clips', '6']),
    ('proper_clips_left', ['--seed', '5', '--n-split-left', '1',
                           '--n-split-right', '0', '--n-proper-clips', '6']),
    # The poly-A of the right junction reads aligned as a
    # supplementary piece of low complexity: It is left out
    ('polya_supplementary', ['--seed', '3', '--polya-len', '40',
                             '--polya-supplementary']),
]


def run(cmd, log):
    log.write('$ %s\n' % ' '.join(cmd))
    log.flush()
    return subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode


def scenario(sider, sim, name, args, log):
    with tempfile.TemporaryDirectory(prefix='sider_sim_') as d:
        steps = [
            [sim, '-o', d] + args,
            [sider, 'process-sample', '-q', '-s', '-a', os.path.join(d, 'genes.gtf'),
             '-o', d, os.path.join(d, 'sim.sam')],
            [sider, 'merge-call', '-q', '-o', d, '-p', 'merged',
             '-r', os.path.join(d, 'ref.fa'), os.path.join(d, 'out.db')],
            [sider, 'make-vcf', '-q', '-r', os.path.join(d, 'ref.fa'),
             '-o', d, '-p', 'merged', os.path.join(d, 'merged.db')],
        ]

        for cmd in steps:
            if run(cmd, log) != 0:
                return False

        check = [sim, '--check', os.path.join(d, 'merged.vcf'),
                 os.path.join(d, 'truth.tsv')]

        return run(check, log) == 0


def shuffle_ids(db_file, seed):
    """Renumber the alignments at random, keeping the reads."""
    rnd = random.Random(seed)
    db = sqlite3.connect(db_file)
    ids = [r[0] for r in db.execute('SELECT id FROM alignment')]
    new = list(range(100000, 100000 + 7 * len(ids), 7))
    rnd.shuffle(new)
    m = dict(zip(ids, new))
    # Negative first: Avoid colliding with the ids not renumbered yet
    for old, n in m.items():
        db.execute('UPDATE alignment SET id = ? WHERE id = ?', (-n, old))
        db.execute('UPDATE overlapping SET alignment_id = ? WHERE alignment_id = ?',
                   (-n, old))
    db.execute('UPDATE alignment SET id = -id')
    db.execute('UPDATE overlapping SET alignment_id = -alignment_id')
    db.commit()
    db.close()


def results(d, prefix):
    """The VCF records and the retrocopy results stored."""
    h = hashlib.md5()
    with open(os.path.join(d, prefix + '.vcf')) as fp:
        h.update(''.join(l for l in fp if not l.startswith('##')).encode())
    db = sqlite3.connect(os.path.join(d, prefix + '.db'))
    for t in ('retrocopy', 'insertion', 'cluster', 'genotype'):
        for row in db.execute('SELECT * FROM %s ORDER BY 1, 2' % t):
            h.update(repr(row).encode())
    db.close()
    return h.hexdigest()


def order_scenario(sider, sim, log):
    """The calls depend only on the reads, not on the alignment ids.

    The ids change with the databases merged, so the same
    database with its alignments renumbered must give the
    same results
    """
    with tempfile.TemporaryDirectory(prefix='sider_sim_') as d:
        steps = [
            [sim, '-o', d, '-s', '4', '--n-discordant', '12', '--n-split-left', '2',
             '--n-split-right', '2', '--n-noise-clips', '30'],
            [sider, 'process-sample', '-q', '-s', '-a', os.path.join(d, 'genes.gtf'),
             '-o', d, '-p', 'ps', os.path.join(d, 'sim.sam')],
        ]
        for cmd in steps:
            if run(cmd, log) != 0:
                return False

        digests = set()
        for k in range(4):
            prefix = 'run%d' % k
            db = os.path.join(d, prefix + '.db')
            shutil.copy(os.path.join(d, 'ps.db'), db)
            if k:
                shuffle_ids(db, k)
            for cmd in [
                [sider, 'merge-call', '-q', '-I', '-r', os.path.join(d, 'ref.fa'), db],
                [sider, 'make-vcf', '-q', '-r', os.path.join(d, 'ref.fa'), '-o', d,
                 '-p', prefix, db],
            ]:
                if run(cmd, log) != 0:
                    return False
            digests.add(results(d, prefix))

        log.write('Distinct results with shuffled ids: %d\n' % len(digests))
        return len(digests) == 1


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 2

    sider, sim = sys.argv[1:3]
    failed = []

    with tempfile.NamedTemporaryFile('w', prefix='sider_sim_', suffix='.log',
                                     delete=False) as log:
        for name, args in SCENARIOS:
            log.write('==== %s\n' % name)
            ok = scenario(sider, sim, name, args, log)
            print('%-32s %s' % (name, 'ok' if ok else 'FAIL'))
            if not ok:
                failed.append(name)

        name = 'alignment ids shuffled'
        log.write('==== %s\n' % name)
        ok = order_scenario(sider, sim, log)
        print('%-32s %s' % (name, 'ok' if ok else 'FAIL'))
        if not ok:
            failed.append(name)

    if failed:
        print('\n%d scenario(s) failed. See %s' % (len(failed), log.name))
        return 1

    os.unlink(log.name)
    print('\nAll %d scenarios passed' % (len(SCENARIOS) + 1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
