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
Design the structural variations of the SANDY benchmark.

The events are written as a SANDY genomic variation TSV: SANDY
replaces the 'reference' sequence at 'position' (1-based) by the
'alternate' one. Each insertion carries an anchor base or the
target site itself, so that the three target sites are built as:

    DUP of L   REF = S         ALT = S R S    S: L bases at the cut
    BLUNT      REF = b         ALT = b R      b: base before the cut
    DEL of d   REF = b D       ALT = b R      D: d bases at the cut

where R is the inserted sequence. The expected breakpoints, as
sideRETRO reports them, are written to the truth table.

Classes (the counts and distributions are at DEFAULT_CONFIG and
can be overridden by --config, a JSON file with the same keys):

  retrocopy    Spliced transcript of a protein coding gene, truncated
               at the 5' end, with a poly-A tail, inserted at an L1
               endonuclease motif (TTTT^AA or TT^AAAA). The true
               positives.
  alu, l1      Mobile element consensus (5' truncated for L1) with a
               poly-A tail and a TSD, at an L1 endonuclease motif.
  cnv          Tandem duplication or deletion spanning exons.
  unprocessed  Genomic segment of a gene, introns included, inserted
               blunt far from it.

The decoys must not be called as retrocopies. The output is
deterministic for a given seed and inputs.

Output files (into --output-dir):

  events.tsv    SANDY genomic variation TSV
  truth.tsv     One row per event with the expected values
  config.json   The configuration used

Example:

  $ tests/bench/sandy_events.py --reference hg38.fa \\
      --annotation gencode.gtf.gz --alu ALU.fa --l1 LINE1.fa \\
      --seed 17 --output-dir bench
"""

import argparse
import bisect
import copy
import gzip
import json
import os
import random
import re
import sys

PRIMARY = ['chr%d' % i for i in range(1, 23)] + ['chrX']

DEFAULT_CONFIG = {
	# Keep the events apart, so that their reads do not mix
	'min_event_distance': 10000,
	# Bases without N around an event
	'clean_flank': 2000,
	'retrocopy': {
		'n': 30,
		# Number of events of each target site
		'tsite': {'dup': 20, 'blunt': 3, 'del': 4, 'long_del': 3},
		'dup_len': [7, 20],
		'del_len': [1, 50],
		'long_del_len': [1000, 20000],
		# Fraction of full length insertions. The others retain
		# a uniform length from min_retained to the full length
		'full_length_frac': 0.33,
		'min_retained': 200,
		'min_transcript_len': 300,
		# Number of single exon parental genes
		'single_exon': 3,
		# Number of insertions without poly-A
		'no_polya': 3,
		'polya_len': [10, 60],
		'polya_interrupt_frac': 0.3,
		'hetero_frac': 0.67,
		# Number of insertion sites of each context
		'context': {'intergenic': 15, 'intronic': 12, 'exonic': 3},
		# Minimum distance to the parental gene at the same chromosome
		'min_parental_distance': 1000000
	},
	'alu': {
		'n': 30,
		'tsd_len': [7, 20],
		'polya_len': [10, 40],
		'divergence': 0.01,
		'hetero_frac': 0.67
	},
	'l1': {
		'n': 20,
		'tsd_len': [7, 20],
		'polya_len': [10, 60],
		'divergence': 0.005,
		# Most L1 insertions are 5' truncated: The retained length
		# is min_retained plus an exponential of mean mean_extra
		'full_length_frac': 0.1,
		'min_retained': 300,
		'mean_extra': 900,
		'hetero_frac': 0.67
	},
	'cnv': {
		'n_dup': 5,
		'n_del': 5,
		'len': [1000, 100000],
		'hetero_frac': 0.67
	},
	'unprocessed': {
		'n': 10,
		'len': [1000, 20000],
		'flank': 200,
		'hetero_frac': 0.67
	}
}

TRUTH_COLUMNS = [
	'id', 'class', 'chr', 'tsv_pos', 'ref_len', 'alt_len', 'lbkp',
	'rbkp', 'tsdlen', 'tsite', 'orientation', 'zygosity', 'gene',
	'transcript', 'exons', 'transcript_len', 'retained_len',
	'retained_exons', 'polya_len', 'context', 'span_start',
	'span_end'
]

COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

def revcomp(seq):
	return seq.translate(COMPLEMENT)[::-1]

def die(msg):
	sys.exit('%s: Err: %s' % (os.path.basename(sys.argv[0]), msg))

def merge_config(base, override):
	for key, value in override.items():
		if key not in base:
			die("Unknown config key '%s'" % key)
		if isinstance(base[key], dict) and key != 'context' and key != 'tsite':
			merge_config(base[key], value)
		else:
			base[key] = value
	return base

class Fasta:
	"""Random access to an uncompressed FASTA by its '.fai' index."""

	def __init__(self, path):
		fai = path + '.fai'
		if path.endswith('.gz'):
			die("'%s': Use an uncompressed FASTA" % path)
		if not os.path.exists(fai):
			die("'%s' not found: Index the FASTA with 'samtools faidx'" % fai)
		self.fh = open(path, 'rb')
		self.index = {}
		self.order = []
		with open(fai) as fh:
			for line in fh:
				name, length, offset, lbases, lwidth = line.split('\t')[:5]
				self.index[name] = (int(length), int(offset), int(lbases),
						int(lwidth))
				self.order.append(name)

	def length(self, chrom):
		return self.index[chrom][0]

	def fetch(self, chrom, start, end):
		"""Sequence at [start, end), 0-based, uppercase."""
		length, offset, lbases, lwidth = self.index[chrom]
		start = max(0, start)
		end = min(length, end)
		if end <= start:
			return ''
		first = offset + (start // lbases) * lwidth + start % lbases
		last = offset + ((end - 1) // lbases) * lwidth + (end - 1) % lbases
		self.fh.seek(first)
		data = self.fh.read(last - first + 1).decode('ascii')
		return data.replace('\n', '').replace('\r', '').upper()

def read_fasta_first(path):
	"""Name and sequence of the first record of a small FASTA."""
	opener = gzip.open if path.endswith('.gz') else open
	name, seq = None, []
	with opener(path, 'rt') as fh:
		for line in fh:
			line = line.strip()
			if line.startswith('>'):
				if name is not None:
					break
				name = line[1:].split()[0]
			elif name is not None:
				seq.append(line)
	if name is None:
		die("'%s': No sequence" % path)
	return name, ''.join(seq).upper()

class Intervals:
	"""Overlap queries over intervals [start, end), 0-based."""

	def __init__(self, items):
		items.sort(key=lambda x: (x[0], x[1]))
		self.starts = [x[0] for x in items]
		self.items = items
		self.max_end = []
		m = -1
		for x in items:
			m = max(m, x[1])
			self.max_end.append(m)

	def overlap(self, start, end):
		hits = []
		i = bisect.bisect_left(self.starts, end) - 1
		while i >= 0 and self.max_end[i] > start:
			if self.items[i][1] > start:
				hits.append(self.items[i])
			i -= 1
		return hits

class Annotation:
	"""Protein coding genes and their canonical transcripts."""

	ATTR = re.compile(r'(\S+) "([^"]*)"')

	def __init__(self, path, chroms):
		opener = gzip.open if path.endswith('.gz') else open
		tx = {}
		genes = []
		with opener(path, 'rt') as fh:
			for line in fh:
				if line.startswith('#'):
					continue
				f = line.rstrip('\n').split('\t')
				if len(f) < 9 or f[0] not in chroms:
					continue
				if f[2] != 'gene' and f[2] != 'exon':
					continue
				attr = dict(self.ATTR.findall(f[8]))
				if attr.get('gene_type') != 'protein_coding':
					continue
				name = attr.get('gene_name', attr.get('gene_id'))
				start, end = int(f[3]) - 1, int(f[4])
				if f[2] == 'gene':
					genes.append((f[0], start, end, name))
					continue
				if attr.get('transcript_type') != 'protein_coding':
					continue
				tid = attr['transcript_id']
				t = tx.get(tid)
				if t is None:
					t = tx[tid] = {
						'id': tid, 'gene': name, 'chr': f[0],
						'strand': f[6], 'exons': [],
						'canonical': 'Ensembl_canonical' in f[8]
					}
				t['exons'].append((start, end))

		# One transcript per gene: The canonical or the longest
		best = {}
		for t in tx.values():
			t['exons'].sort()
			t['len'] = sum(e - s for s, e in t['exons'])
			key = (t['canonical'], t['len'], t['id'])
			if t['gene'] not in best or key > best[t['gene']][0]:
				best[t['gene']] = (key, t)
		self.transcripts = sorted((v[1] for v in best.values()),
				key=lambda t: t['gene'])
		if not self.transcripts:
			die("No protein coding transcript at '%s'" % path)

		self.genes = {}
		self.exons = {}
		by_chr_genes = {}
		by_chr_exons = {}
		for c, s, e, n in genes:
			by_chr_genes.setdefault(c, []).append((s, e, n))
		for t in tx.values():
			for s, e in t['exons']:
				by_chr_exons.setdefault(t['chr'], []).append((s, e, t['gene']))
		for c in chroms:
			self.genes[c] = Intervals(by_chr_genes.get(c, []))
			self.exons[c] = Intervals(by_chr_exons.get(c, []))

	def context(self, chrom, pos):
		if self.exons[chrom].overlap(pos, pos + 1):
			return 'exonic'
		if self.genes[chrom].overlap(pos, pos + 1):
			return 'intronic'
		return 'intergenic'

class Designer:

	def __init__(self, args, config):
		self.rng = random.Random(args.seed)
		self.cfg = config
		self.fa = Fasta(args.reference)
		self.chroms = [c for c in PRIMARY if c in self.fa.index]
		if not self.chroms:
			die("No primary chromosome (chr1..chr22, chrX) at '%s'"
					% args.reference)
		self.weights = [self.fa.length(c) for c in self.chroms]
		self.ann = Annotation(args.annotation, set(self.chroms))
		self.alu = read_fasta_first(args.alu) if args.alu else None
		self.l1 = read_fasta_first(args.l1) if args.l1 else None
		# Occupied regions: (chr, start, end)
		self.occupied = []
		# Genes already used by an event
		self.used_genes = set()
		self.events = []

	# Helpers

	def zygosity(self, frac):
		return 'HE' if self.rng.random() < frac else 'HO'

	def uniform(self, lim):
		return self.rng.randint(lim[0], lim[1])

	def polya(self, lim, interrupt_frac=0.0):
		n = self.uniform(lim)
		tail = ['A'] * n
		if n > 8 and self.rng.random() < interrupt_frac:
			tail[self.rng.randint(3, n - 4)] = self.rng.choice('CGT')
		return ''.join(tail)

	def mutate(self, seq, divergence):
		seq = list(seq)
		for i in range(len(seq)):
			if self.rng.random() < divergence:
				seq[i] = self.rng.choice([b for b in 'ACGT' if b != seq[i]])
		return ''.join(seq)

	def is_free(self, chrom, start, end):
		d = self.cfg['min_event_distance']
		for c, s, e in self.occupied:
			if c == chrom and start - d < e and end + d > s:
				return False
		return True

	def is_clean(self, chrom, start, end):
		f = self.cfg['clean_flank']
		if start - f < 0 or end + f > self.fa.length(chrom):
			return False
		return 'N' not in self.fa.fetch(chrom, start - f, end + f)

	def find_cut(self, chrom, start, end):
		"""First L1 endonuclease cut at [start, end): TTTT^AA, TT^AAAA."""
		seq = self.fa.fetch(chrom, start, end + 6)
		cuts = []
		m = re.search('TTTTAA', seq)
		if m:
			cuts.append(start + m.start() + 4)
		m = re.search('TTAAAA', seq)
		if m:
			cuts.append(start + m.start() + 2)
		cuts = [c for c in cuts if c < end]
		return min(cuts) if cuts else None

	def random_window(self, context):
		"""A window of 500 bases at the asked context."""
		if context == 'exonic':
			t = self.rng.choice(self.ann.transcripts)
			s, e = self.rng.choice(t['exons'])
			return t['chr'], s, e
		if context == 'intronic':
			t = self.rng.choice(self.ann.transcripts)
			if len(t['exons']) < 2:
				return None
			i = self.rng.randrange(len(t['exons']) - 1)
			s, e = t['exons'][i][1], t['exons'][i + 1][0]
			if e - s < 1000:
				return None
			p = self.rng.randrange(s, e - 500)
			return t['chr'], p, p + 500
		chrom = self.rng.choices(self.chroms, self.weights)[0]
		p = self.rng.randrange(0, self.fa.length(chrom) - 500)
		return chrom, p, p + 500

	def pick_site(self, context, span, parent=None, motif=True):
		"""Cut position (0-based) for an event of 'span' reference bases."""
		for _ in range(200000):
			w = self.random_window(context)
			if w is None:
				continue
			chrom, s, e = w
			cut = self.find_cut(chrom, s, e) if motif \
				else self.rng.randrange(s, e)
			if cut is None:
				continue
			if self.ann.context(chrom, cut) != context:
				continue
			if parent is not None and parent['chr'] == chrom:
				ps, pe = parent['exons'][0][0], parent['exons'][-1][1]
				d = self.cfg['retrocopy']['min_parental_distance']
				if cut + span > ps - d and cut < pe + d:
					continue
			if not self.is_free(chrom, cut - 1, cut + span):
				continue
			if not self.is_clean(chrom, cut - 1, cut + span):
				continue
			return chrom, cut
		die("Failed to find a site at context '%s'" % context)

	def transcript_seq(self, t):
		seq = ''.join(self.fa.fetch(t['chr'], s, e) for s, e in t['exons'])
		return seq if t['strand'] == '+' else revcomp(seq)

	def retained_exons(self, t, retained):
		"""Number of exons at the last 'retained' transcript bases."""
		lens = [e - s for s, e in t['exons']]
		if t['strand'] == '-':
			lens.reverse()
		acm, n = 0, 0
		for length in reversed(lens):
			n += 1
			acm += length
			if acm >= retained:
				break
		return n

	# Insertions at a target site

	def insertion(self, cls, chrom, cut, tsite, size, seq, orientation,
			zygosity, extra):
		"""Add an insertion of 'seq' (sense) at 'cut' (0-based)."""
		r = seq if orientation == '+' else revcomp(seq)
		if tsite == 'DUP':
			ref = self.fa.fetch(chrom, cut, cut + size)
			alt = ref + r + ref
			pos, lbkp, rbkp, tsdlen = cut + 1, cut + size, cut + 1, size
			span = (cut, cut + size)
		else:
			d = size if tsite == 'DEL' else 0
			ref = self.fa.fetch(chrom, cut - 1, cut + d)
			alt = ref[0] + r
			pos, lbkp, rbkp, tsdlen = cut, cut, cut + d + 1, -d
			span = (cut - 1, cut + d)
		self.occupied.append((chrom, span[0], span[1]))
		ev = {
			'class': cls, 'chr': chrom, 'tsv_pos': pos, 'ref': ref,
			'alt': alt, 'lbkp': lbkp, 'rbkp': rbkp, 'tsdlen': tsdlen,
			'tsite': tsite, 'orientation': orientation,
			'zygosity': zygosity, 'context': self.ann.context(chrom, cut),
			'span_start': span[0] + 1, 'span_end': span[1]
		}
		ev.update(extra)
		self.events.append(ev)

	def design_retrocopies(self):
		c = self.cfg['retrocopy']
		tsites = []
		for key, name in (('dup', 'DUP'), ('blunt', 'BLUNT'),
				('del', 'DEL'), ('long_del', 'DEL')):
			lim = {'dup': c['dup_len'], 'blunt': [0, 0], 'del': c['del_len'],
					'long_del': c['long_del_len']}[key]
			tsites += [(name, lim)] * c['tsite'][key]
		if len(tsites) != c['n']:
			die("retrocopy: The target sites sum %d, not n=%d"
					% (len(tsites), c['n']))
		self.rng.shuffle(tsites)

		single = [t for t in self.ann.transcripts if len(t['exons']) == 1
				and t['len'] >= c['min_transcript_len']]
		multi = [t for t in self.ann.transcripts if len(t['exons']) > 1
				and t['len'] >= c['min_transcript_len']]
		parents = self.rng.sample(single, c['single_exon']) \
			+ self.rng.sample(multi, c['n'] - c['single_exon'])
		self.rng.shuffle(parents)

		no_polya = set(self.rng.sample(range(c['n']), c['no_polya']))
		contexts = [k for k in sorted(c['context'])
				for _ in range(c['context'][k])]
		if len(contexts) != c['n']:
			die("retrocopy: The contexts sum %d, not n=%d"
					% (len(contexts), c['n']))
		self.rng.shuffle(contexts)

		for i, (t, (tsite, lim)) in enumerate(zip(parents, tsites)):
			self.used_genes.add(t['gene'])
			size = self.uniform(lim)
			full = self.transcript_seq(t)
			if self.rng.random() < c['full_length_frac'] \
					or len(full) <= c['min_retained']:
				retained = len(full)
			else:
				retained = self.rng.randint(c['min_retained'], len(full))
			tail = '' if i in no_polya \
				else self.polya(c['polya_len'], c['polya_interrupt_frac'])
			context = contexts[i]
			span = size if tsite != 'BLUNT' else 1
			chrom, cut = self.pick_site(context, span, parent=t)
			self.insertion('retrocopy', chrom, cut, tsite, size,
					full[len(full) - retained:] + tail,
					self.rng.choice('+-'), self.zygosity(c['hetero_frac']), {
						'gene': t['gene'], 'transcript': t['id'],
						'exons': len(t['exons']), 'transcript_len': len(full),
						'retained_len': retained,
						'retained_exons': self.retained_exons(t, retained),
						'polya_len': len(tail)
					})

	def design_mobile(self, cls, element):
		c = self.cfg[cls]
		if c['n'] == 0:
			return
		if element is None:
			die("'%s' events need its consensus (--%s)" % (cls, cls))
		name, consensus = element
		consensus = consensus.rstrip('A')
		for _ in range(c['n']):
			if cls == 'l1' and self.rng.random() >= c['full_length_frac']:
				retained = min(len(consensus), c['min_retained']
						+ int(self.rng.expovariate(1.0 / c['mean_extra'])))
			else:
				retained = len(consensus)
			seq = self.mutate(consensus[len(consensus) - retained:],
					c['divergence'])
			tail = self.polya(c['polya_len'])
			size = self.uniform(c['tsd_len'])
			context = self.rng.choice(['intergenic', 'intergenic', 'intronic'])
			chrom, cut = self.pick_site(context, size)
			self.insertion(cls, chrom, cut, 'DUP', size, seq + tail,
					self.rng.choice('+-'), self.zygosity(c['hetero_frac']), {
						'gene': name, 'transcript_len': len(consensus),
						'retained_len': retained, 'polya_len': len(tail)
					})

	def free_transcripts(self, min_exons):
		return [t for t in self.ann.transcripts
				if len(t['exons']) >= min_exons and t['gene'] not in self.used_genes]

	def design_cnvs(self):
		c = self.cfg['cnv']
		kinds = ['DUP'] * c['n_dup'] + ['DEL'] * c['n_del']
		for kind in kinds:
			for _ in range(100000):
				t = self.rng.choice(self.free_transcripts(2))
				s, e = self.rng.choice(t['exons'])
				length = max(self.uniform(c['len']), e - s + 2)
				a = self.rng.randint(max(1, e - length), s)
				b = a + length
				chrom = t['chr']
				if b >= self.fa.length(chrom) or not self.is_free(chrom, a - 1, b) \
						or not self.is_clean(chrom, a - 1, b):
					continue
				break
			else:
				die('Failed to find a CNV')
			self.used_genes.add(t['gene'])
			if kind == 'DUP':
				ref = self.fa.fetch(chrom, b - 1, b)
				alt = ref + self.fa.fetch(chrom, a, b)
				pos = b
			else:
				ref = self.fa.fetch(chrom, a - 1, b)
				alt = ref[0]
				pos = a
			self.occupied.append((chrom, a - 1, b))
			exons = sum(1 for s, e in t['exons'] if s < b and e > a)
			self.events.append({
				'class': 'cnv', 'chr': chrom, 'tsv_pos': pos, 'ref': ref,
				'alt': alt, 'tsite': kind, 'zygosity': self.zygosity(c['hetero_frac']),
				'gene': t['gene'], 'transcript': t['id'], 'exons': exons,
				'context': 'exonic', 'span_start': a + 1, 'span_end': b
			})

	def design_unprocessed(self):
		c = self.cfg['unprocessed']
		for _ in range(c['n']):
			for _ in range(100000):
				t = self.rng.choice(self.free_transcripts(3))
				i = self.rng.randrange(len(t['exons']) - 1)
				j = self.rng.randrange(i + 1, len(t['exons']))
				a = t['exons'][i][0] - c['flank']
				b = t['exons'][j][1] + c['flank']
				if c['len'][0] <= b - a <= c['len'][1] \
						and 'N' not in self.fa.fetch(t['chr'], a, b):
					break
			else:
				die('Failed to find an unprocessed segment')
			self.used_genes.add(t['gene'])
			seq = self.fa.fetch(t['chr'], a, b)
			if t['strand'] == '-':
				seq = revcomp(seq)
			chrom, cut = self.pick_site(
					self.rng.choice(['intergenic', 'intronic']), 1,
					parent=t, motif=False)
			self.insertion('unprocessed', chrom, cut, 'BLUNT', 0, seq,
					self.rng.choice('+-'), self.zygosity(c['hetero_frac']), {
						'gene': t['gene'], 'transcript': t['id'],
						'exons': j - i + 1, 'transcript_len': b - a,
						'retained_len': b - a
					})

	def design(self):
		self.design_retrocopies()
		self.design_mobile('alu', self.alu)
		self.design_mobile('l1', self.l1)
		self.design_cnvs()
		self.design_unprocessed()

		rank = {c: i for i, c in enumerate(self.fa.order)}
		self.events.sort(key=lambda ev: (rank[ev['chr']], ev['tsv_pos']))
		count = {}
		for ev in self.events:
			count[ev['class']] = count.get(ev['class'], 0) + 1
			ev['id'] = '%s_%03d' % (ev['class'].upper(), count[ev['class']])
		return count

	def write(self, outdir):
		with open(os.path.join(outdir, 'events.tsv'), 'w') as fh:
			fh.write('#seqid\tposition\tid\treference\talternate\tgenotype\n')
			for ev in self.events:
				fh.write('\t'.join([ev['chr'], str(ev['tsv_pos']), ev['id'],
					ev['ref'], ev['alt'], ev['zygosity']]) + '\n')
		with open(os.path.join(outdir, 'truth.tsv'), 'w') as fh:
			fh.write('\t'.join(TRUTH_COLUMNS) + '\n')
			for ev in self.events:
				ev['ref_len'] = len(ev['ref'])
				ev['alt_len'] = len(ev['alt'])
				fh.write('\t'.join(str(ev.get(k, '.')) for k in TRUTH_COLUMNS)
						+ '\n')
		with open(os.path.join(outdir, 'config.json'), 'w') as fh:
			json.dump(self.cfg, fh, indent=2)
			fh.write('\n')

def main():
	p = argparse.ArgumentParser(
		description='Design the SANDY benchmark events',
		formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog=__doc__)
	p.add_argument('--reference', required=True,
		help='Uncompressed reference FASTA, indexed (.fai)')
	p.add_argument('--annotation', required=True,
		help='GENCODE GTF, the same given to sideRETRO')
	p.add_argument('--alu', help='Alu consensus FASTA (first record)')
	p.add_argument('--l1', help='L1 consensus FASTA (first record)')
	p.add_argument('--config', help='JSON overriding DEFAULT_CONFIG')
	p.add_argument('--seed', type=int, default=17)
	p.add_argument('--output-dir', default='.')
	args = p.parse_args()

	config = copy.deepcopy(DEFAULT_CONFIG)
	if args.config:
		with open(args.config) as fh:
			merge_config(config, json.load(fh))

	os.makedirs(args.output_dir, exist_ok=True)
	d = Designer(args, config)
	count = d.design()
	d.write(args.output_dir)
	print('Designed %d events: %s' % (len(d.events),
		', '.join('%s=%d' % kv for kv in sorted(count.items()))),
		file=sys.stderr)

if __name__ == '__main__':
	main()
