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

"""Compare the results of two sideRETRO runs.

A change that must not alter the calls (e.g. a smaller database) is
checked by running the old and the new version on the same data and
comparing their results:

  - VCF: All records and header lines, but the date and the version
  - Read names (order does not matter) and predicted insertions
  - The summaries of process-sample, merge-call and make-vcf, but
    the elapsed time and the paths
  - The number of rows of each database table (informative: Some
    changes are expected to store less)
  - The database size

Usage:

  Two runs already done (e.g. with the docker images of each version):

    compare_runs.py OLD_DIR NEW_DIR

  Each directory holds the files of one run: <prefix>.vcf,
  <prefix>.rnames.tsv, <prefix>.insertions.fa, <prefix>.db and the
  logs ps.log, mc.log and vcf.log (the names given to --log-file).

  Simulated data, run with two sider executables:

    compare_runs.py --sim OLD_SIDER NEW_SIDER tests/sim/retrocopy_sim.py

Exit status: 0 if the calls are identical, 1 otherwise.
"""

import argparse
import os
import re
import sqlite3
import subprocess
import sys
import tempfile

# Header lines that change from run to run
VCF_VOLATILE = ('##fileDate=', '##source=')

# Summary lines that change from run to run
SUMMARY_VOLATILE = re.compile(r'Elapsed time|File |Database |Next step|VCF |'
                              r'Read names|Predicted insertions')

LOGS = ('ps.log', 'mc.log', 'vcf.log')


def read_lines(path):
    if not os.path.exists(path):
        return None
    with open(path) as fp:
        return fp.read().splitlines()


def vcf_lines(path):
    lines = read_lines(path)
    if lines is None:
        return None
    return [l for l in lines if not l.startswith(VCF_VOLATILE)]


def summary_lines(path):
    # The summary block: from its title to the closing rule
    lines = read_lines(path)
    if lines is None:
        return None

    out = []
    inside = False

    for l in lines:
        if re.match(r'^=+ .* summary =+$', l):
            inside = True
        elif inside and re.match(r'^=+$', l):
            inside = False
        elif inside and l.strip() and not SUMMARY_VOLATILE.search(l):
            out.append(l)

    return out


def table_counts(path):
    if not os.path.exists(path):
        return None

    db = sqlite3.connect('file:%s?mode=ro' % path, uri=True)
    tables = [r[0] for r in db.execute(
        "SELECT name FROM sqlite_master WHERE type = 'table' ORDER BY name")]
    counts = {t: db.execute('SELECT COUNT(*) FROM "%s"' % t).fetchone()[0]
              for t in tables}
    db.close()

    return counts


def find(d, suffix):
    # <prefix><suffix>: The prefix is the user's choice
    hits = sorted(f for f in os.listdir(d) if f.endswith(suffix))
    return os.path.join(d, hits[0]) if hits else None


def diff(name, old, new, report):
    if old is None and new is None:
        return True

    if old is None or new is None:
        report.append('%s: present in only one run' % name)
        return False

    if old == new:
        return True

    o, n = set(old), set(new)
    report.append('%s: differ' % name)
    for l in sorted(o - n)[:10]:
        report.append('  - %s' % l[:160])
    for l in sorted(n - o)[:10]:
        report.append('  + %s' % l[:160])
    if o == n:
        report.append('  (same lines, different order)')

    return False


def compare_dirs(old_dir, new_dir, out=sys.stdout):
    report = []
    same = True

    # The calls: They must be identical
    same &= diff('VCF', vcf_lines(find(old_dir, '.vcf') or ''),
                 vcf_lines(find(new_dir, '.vcf') or ''), report)

    old_rn = read_lines(find(old_dir, '.rnames.tsv') or '')
    new_rn = read_lines(find(new_dir, '.rnames.tsv') or '')
    same &= diff('Read names',
                 sorted(old_rn) if old_rn is not None else None,
                 sorted(new_rn) if new_rn is not None else None, report)

    same &= diff('Predicted insertions',
                 read_lines(find(old_dir, '.insertions.fa') or ''),
                 read_lines(find(new_dir, '.insertions.fa') or ''), report)

    # The summaries: A difference must be explained
    for log in LOGS:
        same &= diff('Summary of %s' % log,
                     summary_lines(os.path.join(old_dir, log)),
                     summary_lines(os.path.join(new_dir, log)), report)

    print('== Calls and summaries', file=out)
    print('\n'.join(report) if report else 'Identical', file=out)

    # The database: Informative
    old_db, new_db = find(old_dir, '.db'), find(new_dir, '.db')

    if old_db and new_db:
        oc, nc = table_counts(old_db), table_counts(new_db)
        print('\n== Database rows (old -> new)', file=out)
        for t in sorted(set(oc) | set(nc)):
            o, n = oc.get(t, '-'), nc.get(t, '-')
            print('  %-24s %12s %12s%s' % (t, o, n, '' if o == n else '  *'),
                  file=out)

        os_, ns = os.path.getsize(old_db), os.path.getsize(new_db)
        print('\n== Database size\n  %.1f MiB -> %.1f MiB (%.1f%%)'
              % (os_ / 2**20, ns / 2**20, 100.0 * (ns - os_) / os_), file=out)

    print('\n%s' % ('SAME CALLS' if same else 'DIFFERENT CALLS'), file=out)
    return same


def run_sim(sider, sim_dir, run_dir, ref):
    # The steps of tests/sim/run_sim_test.py, keeping the logs
    os.makedirs(run_dir)
    steps = [
        [sider, 'process-sample', '-q', '-s', '-l', os.path.join(run_dir, 'ps.log'),
         '-a', os.path.join(sim_dir, 'genes.gtf'), '-o', run_dir, '-p', 'sider',
         os.path.join(sim_dir, 'sim.sam')],
        # In place, as the user's runs: The same database
        [sider, 'merge-call', '-q', '-I', '-l', os.path.join(run_dir, 'mc.log'),
         '-p', 'sider'] + (['-r', ref] if ref else [])
        + [os.path.join(run_dir, 'sider.db')],
        [sider, 'make-vcf', '-q', '-l', os.path.join(run_dir, 'vcf.log'),
         '-r', os.path.join(sim_dir, 'ref.fa'), '-o', run_dir, '-p', 'sider',
         os.path.join(run_dir, 'sider.db')],
    ]

    for cmd in steps:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def compare_sim(old, new, sim):
    # Several simulations: Different seeds and shapes
    scenarios = [
        ['-s', '1'],
        ['-s', '2', '--truncation', '200'],
        ['-s', '3', '--tsd-len', '0', '--polya-len', '0'],
        ['-s', '4', '--n-noise-clips', '20'],
        ['-s', '5', '--n-junction-reads', '0'],
    ]
    same = True

    with tempfile.TemporaryDirectory(prefix='sider_cmp_') as d:
        for i, args in enumerate(scenarios):
            for ref in (True, False):
                sim_dir = os.path.join(d, 'sim%d' % i)
                if not os.path.exists(sim_dir):
                    subprocess.run([sys.executable, sim, '-o', sim_dir] + args,
                                   check=True, stdout=subprocess.DEVNULL)

                tag = '%d%s' % (i, 'r' if ref else '')
                old_dir = os.path.join(d, 'old' + tag)
                new_dir = os.path.join(d, 'new' + tag)
                r = os.path.join(sim_dir, 'ref.fa') if ref else None
                run_sim(old, sim_dir, old_dir, r)
                run_sim(new, sim_dir, new_dir, r)

                print('#### %s%s' % (' '.join(args), '' if ref else ' (no reference)'))
                same &= compare_dirs(old_dir, new_dir)
                print()

    print('ALL SAME CALLS' if same else 'SOME CALLS DIFFER')
    return same


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sim', action='store_true',
                   help='run two sider executables on simulated data')
    p.add_argument('args', nargs='+', metavar='ARG',
                   help='OLD_DIR NEW_DIR, or OLD_SIDER NEW_SIDER SIM with --sim')
    a = p.parse_args()

    if a.sim:
        if len(a.args) != 3:
            p.error('--sim needs OLD_SIDER NEW_SIDER SIM')
        return 0 if compare_sim(*a.args) else 1

    if len(a.args) != 2:
        p.error('needs OLD_DIR NEW_DIR')

    return 0 if compare_dirs(*a.args) else 1


if __name__ == '__main__':
    sys.exit(main())
