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
Evaluate a sideRETRO VCF against the truth of the SANDY benchmark.

A call is placed at its breakpoints (LBKP, RBKP) or, without them,
at POS + 1. It matches an event when they are at the same chromosome
and the call is within --window bases of the event (its breakpoints
or, for a CNV, its span):

  TP           A retrocopy event whose gene is at the call's PG.
               If several calls match an event, the best one (PASS
               first, then the nearest) is the TP and the others
               are FP 'duplicate'.
  FP <class>   The call matches a decoy (alu, l1, cnv, unprocessed)
               or a retrocopy with another parental gene
               ('wrong_gene').
  FP background  The call matches no event.

With --db (the database of merge-call), a retrocopy event without
any abnormal alignment within --flank bases of its site is left out
of the recall and of the FN: the reads give sideRETRO nothing to
call, so it is a failure of the simulation, not of the caller.

Precision, recall and F1 are given for the PASS calls and for all
the calls. Recall is broken down by the event attributes, and each
filter is counted at the TP and FP calls. The attributes of the TP
calls are compared to the expected ones.

Output files (with --prefix):

  <prefix>.events.tsv  One row per retrocopy event with its best call
  <prefix>.calls.tsv   One row per call with its classification

Example:

  $ tests/bench/bench_eval.py --prefix bench/eval bench/truth.tsv \\
      bench/out.vcf --db bench/out.db
"""

import argparse
import csv
import gzip
import sqlite3
import sys

DECOYS = ['alu', 'l1', 'cnv', 'unprocessed']
LONG_DEL = 1000

def read_truth(path):
	with open(path) as fh:
		rows = list(csv.DictReader(fh, delimiter='\t'))
	for r in rows:
		pts = [int(r['span_start']), int(r['span_end'])]
		for k in ('lbkp', 'rbkp'):
			if r[k] != '.':
				pts.append(int(r[k]))
		r['lo'], r['hi'] = min(pts), max(pts)
	return rows

def read_vcf(path):
	opener = gzip.open if path.endswith('.gz') else open
	calls = []
	with opener(path, 'rt') as fh:
		for line in fh:
			if line.startswith('#'):
				continue
			f = line.rstrip('\n').split('\t')
			info = {}
			for item in f[7].split(';'):
				k, _, v = item.partition('=')
				info[k] = v if _ else True
			gt = '.'
			if len(f) > 9:
				keys = f[8].split(':')
				vals = f[9].split(':')
				if 'GT' in keys:
					gt = vals[keys.index('GT')]
			pts = [int(f[1]) + 1]
			for k in ('LBKP', 'RBKP'):
				if k in info:
					pts.append(int(info[k]))
			calls.append({
				'chr': f[0], 'pos': int(f[1]), 'id': f[2], 'filter': f[6],
				'info': info, 'gt': gt, 'lo': min(pts), 'hi': max(pts),
				'pg': set(info.get('PG', '').split(','))
			})
	return calls

def count_reads(path, truth, flank):
	"""Abnormal alignments stored around each retrocopy event."""
	db = sqlite3.connect('file:%s?mode=ro' % path, uri=True)
	for ev in truth:
		if ev['class'] != 'retrocopy':
			continue
		ev['reads'] = db.execute(
			'SELECT COUNT(*) FROM alignment WHERE chr = ? AND pos BETWEEN ? AND ?',
			(ev['chr'], ev['lo'] - flank, ev['hi'] + flank)).fetchone()[0]
	db.close()

def distance(call, ev):
	if call['hi'] < ev['lo']:
		return ev['lo'] - call['hi']
	if call['lo'] > ev['hi']:
		return call['lo'] - ev['hi']
	return 0

def classify(calls, truth, window):
	by_chr = {}
	for ev in truth:
		by_chr.setdefault(ev['chr'], []).append(ev)
	best = {}
	for call in calls:
		near = [(distance(call, ev), ev) for ev in by_chr.get(call['chr'], [])]
		near = sorted((d, ev['id'], ev) for d, ev in near if d <= window)
		call['event'] = None
		call['class'] = 'background'
		for d, _, ev in near:
			if ev['class'] == 'retrocopy' and ev['gene'] in call['pg']:
				call['event'], call['dist'] = ev, d
				call['class'] = 'TP'
				key = (call['filter'] != 'PASS', d)
				if ev['id'] not in best or key < best[ev['id']][0]:
					best[ev['id']] = (key, call)
				break
		else:
			if near:
				ev = near[0][2]
				call['event'] = ev
				call['class'] = ev['class'] if ev['class'] in DECOYS \
					else 'wrong_gene'
	for call in calls:
		if call['class'] == 'TP' and best[call['event']['id']][1] is not call:
			call['class'] = 'duplicate'
	return {k: v[1] for k, v in best.items()}

def ratio(a, b):
	return a / b if b else float('nan')

def metrics(calls, events, best, only_pass):
	sel = [c for c in calls if not only_pass or c['filter'] == 'PASS']
	tp = sum(1 for ev in events if ev['id'] in best
			and (not only_pass or best[ev['id']]['filter'] == 'PASS'))
	fp = {}
	for c in sel:
		if c['class'] != 'TP':
			fp[c['class']] = fp.get(c['class'], 0) + 1
	nfp = sum(fp.values())
	p, r = ratio(tp, tp + nfp), ratio(tp, len(events))
	f1 = ratio(2 * p * r, p + r)
	return tp, nfp, len(events) - tp, p, r, f1, fp

def found(ev, best, only_pass):
	c = best.get(ev['id'])
	return c is not None and (not only_pass or c['filter'] == 'PASS')

def groups(ev):
	g = []
	tsite = ev['tsite']
	if tsite == 'DEL':
		tsite = 'DEL long' if -int(ev['tsdlen']) >= LONG_DEL else 'DEL short'
	g.append(('tsite', tsite))
	g.append(('zygosity', ev['zygosity']))
	g.append(('context', ev['context']))
	g.append(('exons', 'single' if ev['exons'] == '1' else 'multi'))
	ret, full = int(ev['retained_len']), int(ev['transcript_len'])
	g.append(('length', 'full' if ret == full
		else '< 500' if ret < 500 else '500-1000' if ret <= 1000 else '> 1000'))
	g.append(('poly-A', 'none' if ev['polya_len'] == '0' else 'present'))
	return g

def report(calls, truth, best, out):
	events = [ev for ev in truth if ev['class'] == 'retrocopy']
	no_reads = [ev for ev in events if ev.get('reads', 1) == 0]
	events = [ev for ev in events if ev.get('reads', 1) > 0]
	w = out.write

	if no_reads:
		w('== Left out: %d retrocopies without abnormal reads ==\n'
				% len(no_reads))
		for ev in no_reads:
			w('%s %s:%s %s%s\n' % (ev['id'], ev['chr'], ev['lbkp'], ev['gene'],
				' (called %s)' % best[ev['id']]['id'] if ev['id'] in best else ''))
		w('\n')

	w('== Precision and recall: %d retrocopies, %d calls ==\n'
			% (len(events), len(calls)))
	w('%-6s %4s %4s %4s %9s %9s %9s  %s\n'
			% ('calls', 'TP', 'FP', 'FN', 'precision', 'recall', 'F1', 'FP by class'))
	for label, only_pass in (('PASS', True), ('all', False)):
		tp, fp, fn, p, r, f1, by = metrics(calls, events, best, only_pass)
		w('%-6s %4d %4d %4d %9.3f %9.3f %9.3f  %s\n' % (label, tp, fp, fn,
			p, r, f1, ', '.join('%s=%d' % kv for kv in sorted(by.items()))))

	w('\n== Recall by attribute (PASS / any call) ==\n')
	table = {}
	for ev in events:
		for key in groups(ev):
			t = table.setdefault(key, [0, 0, 0])
			t[0] += 1
			t[1] += found(ev, best, True)
			t[2] += found(ev, best, False)
	for (attr, value), (n, fp_, fa) in sorted(table.items()):
		w('%-9s %-10s %3d  %3d %6.2f  %3d %6.2f\n'
				% (attr, value, n, fp_, ratio(fp_, n), fa, ratio(fa, n)))

	w('\n== Filters at the calls (a call may fail several) ==\n')
	filt = {}
	for c in calls:
		kind = 'TP' if c['class'] == 'TP' else 'FP'
		for name in c['filter'].split(';'):
			f = filt.setdefault(name, {'TP': 0, 'FP': 0})
			f[kind] += 1
	w('%-16s %5s %5s\n' % ('filter', 'TP', 'FP'))
	for name, f in sorted(filt.items(), key=lambda kv: (kv[0] != 'PASS', kv[0])):
		w('%-16s %5d %5d\n' % (name, f['TP'], f['FP']))

	w('\n== Attributes of the TP calls (any filter) ==\n')
	acc = {}
	def tally(name, ok):
		a = acc.setdefault(name, [0, 0])
		a[0] += 1
		a[1] += bool(ok)
	for ev in events:
		c = best.get(ev['id'])
		if c is None:
			continue
		info = c['info']
		for k in ('LBKP', 'RBKP'):
			tally(k + ' reported', k in info)
			if k in info:
				d = abs(int(info[k]) - int(ev[k.lower()]))
				tally(k + ' exact', d == 0)
				tally(k + ' within 5', d <= 5)
		if 'TSDLEN' in info:
			tally('TSDLEN exact', int(info['TSDLEN']) == int(ev['tsdlen']))
			tally('TSITE', info.get('TSITE') == ev['tsite'])
		if 'POLARITY' in info:
			tally('POLARITY', info['POLARITY'] == ev['orientation'])
		gt = c['gt'].replace('|', '/')
		tally('GT', gt == ('1/1' if ev['zygosity'] == 'HO' else '0/1'))
		if ev['polya_len'] != '0':
			tally('POLYA found', 'POLYA' in info)
		else:
			tally('POLYA absent', 'POLYA' not in info)
		if ev['tsite'] == 'DEL' and -int(ev['tsdlen']) >= LONG_DEL:
			tally('TSITEMERGED', 'TSITEMERGED' in info)
	for name, (n, ok) in sorted(acc.items()):
		w('%-18s %3d/%-3d %6.2f\n' % (name, ok, n, ratio(ok, n)))

def write_tables(calls, truth, best, prefix):
	cols = ['id', 'chr', 'tsite', 'tsdlen', 'zygosity', 'context', 'gene', 'exons',
			'retained_len', 'transcript_len', 'polya_len', 'orientation',
			'lbkp', 'rbkp']
	vcols = ['LBKP', 'RBKP', 'TSDLEN', 'TSITE', 'POLARITY', 'POLYA', 'EEJ',
			'TSITEMERGED']
	with open(prefix + '.events.tsv', 'w') as fh:
		fh.write('\t'.join(cols + ['reads', 'call', 'filter', 'gt']
			+ ['call_' + k.lower() for k in vcols]) + '\n')
		for ev in truth:
			if ev['class'] != 'retrocopy':
				continue
			c = best.get(ev['id'])
			row = [ev[k] for k in cols] + [str(ev.get('reads', '.'))]
			if c is None:
				row += ['.'] * (3 + len(vcols))
			else:
				row += [c['id'], c['filter'], c['gt']] \
					+ [str(c['info'].get(k, '.')) for k in vcols]
			fh.write('\t'.join(row) + '\n')
	with open(prefix + '.calls.tsv', 'w') as fh:
		fh.write('id\tchr\tpos\tfilter\tpg\tclass\tevent\tevent_class\n')
		for c in calls:
			ev = c['event']
			fh.write('\t'.join([c['id'], c['chr'], str(c['pos']), c['filter'],
				c['info'].get('PG', '.'), c['class'],
				ev['id'] if ev else '.', ev['class'] if ev else '.']) + '\n')

def main():
	p = argparse.ArgumentParser(
		description='Evaluate a sideRETRO VCF against the benchmark truth',
		formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog=__doc__)
	p.add_argument('truth', help='truth.tsv from sandy_events.py')
	p.add_argument('vcf', help='VCF from sider make-vcf')
	p.add_argument('--window', type=int, default=100,
		help='Maximum distance of a call to an event [default: %(default)s]')
	p.add_argument('--db',
		help='Database of merge-call: leave out the events without reads')
	p.add_argument('--flank', type=int, default=1000,
		help='Distance to look for the reads of an event [default: %(default)s]')
	p.add_argument('--prefix', help='Write <prefix>.events.tsv and .calls.tsv')
	args = p.parse_args()

	truth = read_truth(args.truth)
	calls = read_vcf(args.vcf)
	if args.db:
		count_reads(args.db, truth, args.flank)
	best = classify(calls, truth, args.window)
	report(calls, truth, best, sys.stdout)
	if args.prefix:
		write_tables(calls, truth, best, args.prefix)

if __name__ == '__main__':
	main()
