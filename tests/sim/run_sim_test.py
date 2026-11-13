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
Run the whole pipeline over simulated retrocopies and compare the
results with the simulated truth. Used by 'meson test':

  $ run_sim_test.py SIDER RETROCOPY_SIM

Each scenario runs retrocopy_sim.py, process-sample, merge-call
(with and without the reference) and make-vcf into a temporary
directory, and checks the VCF, the read names and the predicted
insertions files.
"""

import os
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
]


def run(cmd, log):
    log.write('$ %s\n' % ' '.join(cmd))
    log.flush()
    return subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode


def scenario(sider, sim, name, args, with_ref, log):
    with tempfile.TemporaryDirectory(prefix='sider_sim_') as d:
        steps = [
            [sim, '-o', d] + args,
            [sider, 'process-sample', '-q', '-s', '-a', os.path.join(d, 'genes.gtf'),
             '-o', d, os.path.join(d, 'sim.sam')],
            [sider, 'merge-call', '-q', '-o', d, '-p', 'merged']
            + (['-r', os.path.join(d, 'ref.fa')] if with_ref else [])
            + [os.path.join(d, 'out.db')],
            [sider, 'make-vcf', '-q', '-r', os.path.join(d, 'ref.fa'),
             '-o', d, '-p', 'merged', os.path.join(d, 'merged.db')],
        ]

        for cmd in steps:
            if run(cmd, log) != 0:
                return False

        check = [sim, '--check', os.path.join(d, 'merged.vcf'),
                 os.path.join(d, 'truth.tsv')]

        if not with_ref:
            check.append('--no-ref')

        return run(check, log) == 0


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 2

    sider, sim = sys.argv[1:3]
    failed = []

    with tempfile.NamedTemporaryFile('w', prefix='sider_sim_', suffix='.log',
                                     delete=False) as log:
        for name, args in SCENARIOS:
            for with_ref in (True, False):
                full = '%s%s' % (name, '' if with_ref else ' (no reference)')
                log.write('==== %s\n' % full)
                ok = scenario(sider, sim, name, args, with_ref, log)
                print('%-32s %s' % (full, 'ok' if ok else 'FAIL'))
                if not ok:
                    failed.append(full)

    if failed:
        print('\n%d scenario(s) failed. See %s' % (len(failed), log.name))
        return 1

    os.unlink(log.name)
    print('\nAll %d scenarios passed' % (2 * len(SCENARIOS)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
