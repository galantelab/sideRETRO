#!/usr/bin/env python3
"""Generate the figures of the sideRETRO documentation (SVG).

One palette and one style for all of them. Usage:
    make_figures.py OUTPUT_DIR
"""

import os
import sys

# Palette
INK = '#1F2A37'      # text
MUTED = '#64748B'    # secondary text, links
GENE = '#2B7A9B'     # parental gene, strokes
EXON = '#CDE7F0'     # exons, read 1
DARK = '#1B5E7B'     # read 2, poly-A
CHR = '#F2C14E'      # chromosome of the insertion
CHR_D = '#B8860B'
ACCENT = '#D9480F'   # insertion point, breakpoints
TSD = '#F6A46B'      # target site duplication
INTRON = '#7E5BB5'   # retained intron
GREY = '#CBD5E1'     # neutral
BG = '#FFFFFF'
FONT = 'Helvetica, Arial, sans-serif'


class Fig:
    def __init__(self, w, h):
        self.w, self.h = w, h
        self.e = []

    def add(self, s):
        self.e.append(s)

    def rect(self, x, y, w, h, fill, stroke=None, rx=4, sw=1.5, dash=None, op=None):
        a = 'x="%g" y="%g" width="%g" height="%g" rx="%g" fill="%s"' % (x, y, w, h, rx, fill)
        if stroke:
            a += ' stroke="%s" stroke-width="%g"' % (stroke, sw)
        if dash:
            a += ' stroke-dasharray="%s"' % dash
        if op is not None:
            a += ' fill-opacity="%g"' % op
        self.add('<rect %s/>' % a)

    def text(self, x, y, s, size=15, fill=INK, anchor='start', weight='normal', style=None):
        a = 'x="%g" y="%g" font-size="%g" fill="%s" text-anchor="%s"' % (x, y, size, fill, anchor)
        if weight != 'normal':
            a += ' font-weight="%s"' % weight
        if style:
            a += ' font-style="%s"' % style
        s = s.replace('&', '&amp;').replace('<', '&lt;')
        self.add('<text %s>%s</text>' % (a, s))

    def line(self, x1, y1, x2, y2, stroke=MUTED, sw=1.5, dash=None, arrow=False, cap='butt'):
        a = 'x1="%g" y1="%g" x2="%g" y2="%g" stroke="%s" stroke-width="%g" stroke-linecap="%s"' % (
            x1, y1, x2, y2, stroke, sw, cap)
        if dash:
            a += ' stroke-dasharray="%s"' % dash
        if arrow:
            a += ' marker-end="url(#arrow-%s)"' % stroke.strip('#')
            self.need_arrow(stroke)
        self.add('<line %s/>' % a)

    def path(self, d, stroke=MUTED, sw=1.5, dash=None, arrow=False, fill='none'):
        a = 'd="%s" fill="%s" stroke="%s" stroke-width="%g" stroke-linecap="round"' % (d, fill, stroke, sw)
        if dash:
            a += ' stroke-dasharray="%s"' % dash
        if arrow:
            a += ' marker-end="url(#arrow-%s)"' % stroke.strip('#')
            self.need_arrow(stroke)
        self.add('<path %s/>' % a)

    def circle(self, x, y, r, fill, stroke=None, sw=1.5):
        a = 'cx="%g" cy="%g" r="%g" fill="%s"' % (x, y, r, fill)
        if stroke:
            a += ' stroke="%s" stroke-width="%g"' % (stroke, sw)
        self.add('<circle %s/>' % a)

    def need_arrow(self, color):
        if not hasattr(self, 'arrows'):
            self.arrows = set()
        self.arrows.add(color)

    # Shared pieces ---------------------------------------------------

    def read(self, x, y, label, w=40, h=26, dash=None):
        """A read: R1 light, R2 dark, S (supplementary) dark."""
        dark = label in ('R2', 'S')
        self.rect(x, y, w, h, DARK if dark else EXON, GENE, rx=5, dash=dash)
        self.text(x + w / 2, y + h / 2 + 5, label, size=13,
                  fill='#FFFFFF' if dark else INK, anchor='middle', weight='bold')

    def gene(self, y, x0=150, x1=730, label='Parental gene', sub='reference genome',
             exons=((200, 290), (360, 480), (560, 700)), line=GENE):
        """Parental gene: The line (introns) and its exons."""
        self.line(x0, y, x1, y, stroke=line, sw=8, cap='round')
        self.rect(exons[0][0] - 22, y - 12, 22, 24, EXON, GENE, rx=3)
        self.rect(exons[-1][1], y - 12, 22, 24, EXON, GENE, rx=3)
        for a, b in exons:
            self.rect(a, y - 20, b - a, 40, EXON, GENE)
        if label:
            self.text(x1 + 24, y - 1, label, weight='bold')
            self.text(x1 + 24, y + 17, sub, size=12, fill=MUTED)
        return exons

    def retro(self, y, x0=150, x1=730, flank=CHR, label='Retrocopy', sub='sample genome',
              start=280, widths=(90, 120, 140), polya=70):
        """The retrocopy at the sample: Flanks, spliced exons and poly-A."""
        self.line(x0, y, x1, y, stroke=flank, sw=8, cap='round')
        x = start
        self.rect(x - 22, y - 12, 22, 24, EXON, GENE, rx=3)
        bounds = []
        for w in widths:
            bounds.append((x, x + w))
            x += w
        self.rect(start, y - 20, x - start, 40, EXON, GENE)
        for a, b in bounds[1:]:
            self.line(a, y - 20, a, y + 20, stroke=GENE)
        self.rect(x, y - 12, 22, 24, EXON, GENE, rx=3)
        self.rect(x + 22, y - 12, polya, 24, DARK, GENE, rx=3)
        self.text(x + 22 + polya / 2, y + 5, 'AAAA', size=12, fill='#FFFFFF', anchor='middle',
                  weight='bold')
        if label:
            self.text(x1 + 24, y - 1, label, weight='bold')
            self.text(x1 + 24, y + 17, sub, size=12, fill=MUTED)
        return bounds, x + 22 + polya

    def chrom(self, y, x0=150, x1=730, label='Insertion site', sub='reference genome', color=CHR):
        self.line(x0, y, x1, y, stroke=color, sw=8, cap='round')
        if label:
            self.text(x1 + 24, y - 1, label, weight='bold')
            self.text(x1 + 24, y + 17, sub, size=12, fill=MUTED)

    def pair(self, x1, x2, y, w=40):
        """Mates of a fragment, linked."""
        self.line(x1 + w, y + 13, x2, y + 13, stroke=GENE, sw=1.5)

    def link(self, x1, y1, x2, y2, color=MUTED):
        """From where a read was sequenced to where it was aligned."""
        my = (y1 + y2) / 2
        self.path('M%g,%g C%g,%g %g,%g %g,%g' % (x1, y1, x1, my, x2, my, x2, y2),
                  stroke=color, dash='5,4', arrow=True)

    def caption(self, y, s, x=None):
        self.text(self.w / 2 if x is None else x, y, s, size=13, fill=MUTED,
                  anchor='middle' if x is None else 'start')

    def svg(self):
        defs = ''
        for c in sorted(getattr(self, 'arrows', [])):
            defs += ('<marker id="arrow-%s" markerWidth="9" markerHeight="9" refX="7" refY="4.5" '
                     'orient="auto"><path d="M1,1 L8,4.5 L1,8 z" fill="%s"/></marker>'
                     % (c.strip('#'), c))
        return ('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d" '
                'font-family="%s">\n<defs>%s</defs>\n<rect width="%d" height="%d" fill="%s"/>\n%s\n</svg>\n'
                % (self.w, self.h, self.w, self.h, FONT, defs, self.w, self.h, BG, '\n'.join(self.e)))


W = 900
REF_Y = 50      # parental gene at the reference
ALN_Y = 86      # reads aligned under the gene


def alignment_fig(kind):
    """The origin of the reads at the sample and their alignment."""
    two = kind in ('chr', 'dist', 'sr')
    sample_y = 340 if two else 250
    seq_y = sample_y - 56
    site_y = 170
    site_aln = site_y + 16
    f = Fig(W, sample_y + 60)

    f.gene(REF_Y)
    flank = CHR if kind != 'dist' else GENE

    if two:
        if kind == 'dist':
            # The same chromosome: The line goes on, far away
            f.path('M150,%d C96,%d 96,%d 150,%d' % (REF_Y, REF_Y, site_y, site_y),
                   stroke=GENE, sw=8)
            f.chrom(site_y, label='Same chromosome', sub='a distant region', color=GENE)
        else:
            f.chrom(site_y, label='Insertion site', sub='another chromosome')

    bounds, end = f.retro(sample_y, flank=flank)
    f.text(24, REF_Y + 5, 'Aligned', size=13, fill=MUTED, weight='bold')
    f.text(24, sample_y + 5, 'Sequenced', size=13, fill=MUTED, weight='bold')

    def sequenced(x1, x2, l1='R1', l2='R2'):
        f.read(x1, seq_y, l1)
        f.read(x2, seq_y, l2)
        f.pair(x1, x2, seq_y)

    if kind == 'indist':
        # Both mates inside the last exon
        sequenced(540, 620)
        f.read(585, ALN_Y, 'R1')
        f.read(645, ALN_Y, 'R2')
        f.pair(585, 645, ALN_Y)
        f.link(560, seq_y - 4, 605, ALN_Y + 32)
        f.link(640, seq_y - 4, 665, ALN_Y + 32)
    elif kind == 'exon':
        # The fragment spans the junction between two exons
        sequenced(420, 520)
        f.read(410, ALN_Y, 'R1')
        f.read(580, ALN_Y, 'R2')
        f.link(440, seq_y - 4, 430, ALN_Y + 32)
        f.link(540, seq_y - 4, 600, ALN_Y + 32)
    elif kind in ('chr', 'dist'):
        # One mate at the flank, the other inside the retrocopy
        sequenced(180, 300)
        f.read(330, site_aln, 'R1')
        f.read(225, ALN_Y, 'R2')
        f.link(200, seq_y - 4, 350, site_aln + 32)
        f.link(320, seq_y - 4, 245, ALN_Y + 32)
    elif kind == 'sr':
        # One mate crosses the junction: Flank + retrocopy
        f.read(150, seq_y, 'R1')
        f.rect(226, seq_y, 64, 26, DARK, GENE, rx=5)
        f.text(244, seq_y + 18, 'R2', size=13, fill='#FFFFFF', anchor='middle', weight='bold')
        f.line(258, seq_y - 8, 258, seq_y + 34, stroke=ACCENT, sw=2, dash='4,3')
        f.text(296, seq_y + 18, 'crosses the junction', size=12, fill=ACCENT)
        f.pair(150, 226, seq_y)
        # Aligned: Primary with its clip, and the supplementary
        f.read(330, site_aln, 'R1')
        f.rect(410, site_aln, 32, 26, DARK, GENE, rx=5)
        f.text(426, site_aln + 18, 'R2', size=13, fill='#FFFFFF', anchor='middle', weight='bold')
        f.rect(442, site_aln, 32, 26, BG, GENE, rx=5, dash='4,3')
        f.text(484, site_aln + 18, 'clip', size=12, fill=MUTED)
        f.read(206, ALN_Y, 'S', w=32)
        f.text(246, ALN_Y + 18, 'supplementary: the clipped part', size=12, fill=MUTED)
        f.link(170, seq_y - 4, 350, site_aln + 32)
        f.link(250, seq_y - 4, 426, site_aln + 32)
        f.path('M458,%d C458,140 222,150 222,%d' % (site_aln, ALN_Y + 32), stroke=MUTED,
               dash='5,4', arrow=True)

    return f


def clustering_fig():
    f = Fig(W, 470)
    f.gene(50)
    cy = 250
    f.chrom(cy, label='Insertion site', sub='reference genome')
    ip = 440
    # Mates at the exons and the reads piled at the insertion site
    pairs = [(203, 'R2', 300, 196, 'R1'), (385, 'R1', 345, 164, 'R2'), (249, 'R1', 345, 196, 'R2'),
             (571, 'R1', 480, 196, 'R2'), (617, 'R2', 525, 164, 'R1'), (663, 'R2', 525, 196, 'R1')]
    for gx, gl, sx, sy, sl in pairs:
        f.path('M%g,%g C%g,%g %g,%g %g,%g' % (gx + 20, 112, gx + 20, 150, sx + 20, 130, sx + 20, sy),
               stroke=MUTED, sw=1.2)
    for gx, gl, sx, sy, sl in pairs:
        f.read(gx, 86, gl)
        f.read(sx, sy, sl)
    f.rect(288, 154, 290, 78, 'none', GENE, rx=14, dash='6,4')
    f.text(590, 190, 'Cluster of', size=13, fill=GENE, weight='bold')
    f.text(590, 207, 'abnormal reads', size=13, fill=GENE, weight='bold')
    f.line(ip, 236, ip, 264, stroke=ACCENT, sw=3)
    f.text(ip, 288, 'Insertion point', size=13, fill=ACCENT, anchor='middle', weight='bold')
    # The putative retrocopy
    ry = 410
    bounds, end = f.retro(ry, label='Putative retrocopy', sub='sample genome')
    f.path('M%g,300 C%g,350 258,340 258,%g' % (ip, ip, ry - 26), stroke=ACCENT, sw=1.5, dash='5,4')
    f.path('M%g,300 C%g,350 %g,340 %g,%g' % (ip, ip, end, end, ry - 26), stroke=ACCENT, sw=1.5,
           dash='5,4')
    return f


def dbscan_fig():
    f = Fig(W, 330)
    y = 200
    f.line(60, y, 840, y, stroke=CHR, sw=8, cap='round')
    f.text(60, y + 30, 'Chromosome', size=13, fill=MUTED)

    def reads(xs, color, stroke, dash=None):
        for x, row in xs:
            f.rect(x, y - 24 - row * 16, 60, 10, color, stroke, rx=3, sw=1.2, dash=dash)

    c1 = [(110, 0), (128, 1), (140, 2), (152, 3), (166, 0), (176, 4), (186, 1), (204, 2), (222, 0)]
    c2 = [(560, 0), (574, 1), (588, 2), (602, 0), (612, 3), (626, 1), (648, 2), (668, 0)]
    noise = [(350, 0), (420, 0), (770, 0)]
    f.rect(98, y - 100, 196, 92, GENE, None, rx=12, op=0.10)
    f.rect(548, y - 86, 192, 78, GENE, None, rx=12, op=0.10)
    reads(c1, EXON, GENE)
    reads(c2, EXON, GENE)
    reads(noise, BG, MUTED, dash='4,3')
    f.text(196, y - 110, 'Cluster', size=14, fill=GENE, anchor='middle', weight='bold')
    f.text(644, y - 96, 'Cluster', size=14, fill=GENE, anchor='middle', weight='bold')
    f.text(415, y - 44, 'Noise', size=14, fill=MUTED, anchor='middle', weight='bold')
    f.text(800, y - 44, 'Noise', size=14, fill=MUTED, anchor='middle', weight='bold')
    # eps
    f.line(110, y + 52, 170, y + 52, stroke=INK, sw=1.5)
    f.line(110, y + 46, 110, y + 58, stroke=INK)
    f.line(170, y + 46, 170, y + 58, stroke=INK)
    f.text(182, y + 57, 'range of a read: its start and end are the point', size=13, fill=INK)
    f.rect(110, y + 76, 60, 10, EXON, GENE, rx=3, sw=1.2)
    f.text(182, y + 86, 'read in a dense region: at least min-points within eps', size=13)
    f.rect(110, y + 100, 60, 10, BG, MUTED, rx=3, sw=1.2, dash='4,3')
    f.text(182, y + 110, 'read in a sparse region: noise', size=13)
    return f


def genotype_fig():
    f = Fig(W, 380)
    y = 190
    ip = 450
    f.line(60, y, 840, y, stroke=CHR, sw=8, cap='round')
    # One decile of the read length around the insertion point
    f.rect(ip - 14, y - 150, 28, 300, ACCENT, None, rx=0, op=0.10)
    f.line(ip, y - 16, ip, y + 16, stroke=ACCENT, sw=3)
    f.text(ip + 22, y + 32, 'Insertion point', size=13, fill=ACCENT, weight='bold')

    # Abnormal reads: The alternate allele
    for x, row, lab in [(250, 0, 'R1'), (296, 1, 'R2'), (342, 0, 'R2'), (300, 2, 'R1'),
                        (500, 0, 'R1'), (546, 1, 'R2'), (592, 0, 'R1'), (540, 2, 'R2')]:
        f.read(x, y - 46 - row * 32, lab)
    f.rect(238, y - 122, 416, 108, 'none', GENE, rx=14, dash='6,4')
    f.text(64, y - 80, 'Abnormal reads', size=14, weight='bold')
    f.text(64, y - 62, 'evidence of the insertion', size=12, fill=MUTED)

    # Normal reads: Crossing the window or not
    def normal(x, row, w, cross):
        f.rect(x, y + 50 + row * 20, w, 12, DARK if cross else GREY, None, rx=4)
    for x, row, w, cross in [(390, 0, 120, True), (372, 1, 120, True), (414, 2, 120, True),
                             (250, 0, 120, False), (232, 1, 120, False), (330, 3, 112, False),
                             (540, 0, 120, False), (462, 3, 120, False), (566, 1, 120, False)]:
        normal(x, row, w, cross)
    f.text(64, y + 70, 'Normal reads', size=14, weight='bold')
    f.text(64, y + 88, 'proper pairs', size=12, fill=MUTED)
    f.rect(690, y + 52, 26, 12, DARK, None, rx=4)
    f.text(724, y + 63, 'crosses the window:', size=12)
    f.text(724, y + 79, 'reference allele', size=12)
    f.rect(690, y + 100, 26, 12, GREY, None, rx=4)
    f.text(724, y + 111, 'not counted', size=12)
    f.text(ip, y + 168, 'window: one decile of the read length', size=12, fill=ACCENT,
           anchor='middle')
    return f


def line1_fig():
    f = Fig(W, 250)
    y = 120
    f.line(40, y, 860, y, stroke=GENE, sw=8, cap='round')
    f.rect(40, y - 16, 56, 32, TSD, None)
    f.text(68, y + 5, 'TSD', size=13, anchor='middle', weight='bold')
    f.rect(804, y - 16, 56, 32, TSD, None)
    f.text(832, y + 5, 'TSD', size=13, anchor='middle', weight='bold')
    # 5' UTR
    f.rect(140, y - 22, 80, 44, EXON, GENE)
    f.text(180, y + 5, 'ORF0', size=14, anchor='middle', weight='bold')
    f.path('M112,%d L112,%d L150,%d' % (y - 6, y - 44, y - 44), stroke=INK, arrow=True)
    f.text(112, y - 54, 'TSS', size=13, weight='bold')
    f.path('M228,%d L228,%d L190,%d' % (y + 6, y + 44, y + 44), stroke=INK, arrow=True)
    f.text(236, y + 49, 'antisense promoter (ASP)', size=13)
    # ORF1
    f.rect(300, y - 22, 110, 44, EXON, GENE)
    f.text(355, y + 5, 'ORF1', size=14, anchor='middle', weight='bold')
    # ORF2: EN and RT
    f.rect(450, y - 22, 250, 44, EXON, GENE)
    f.rect(450, y - 22, 80, 44, CHR, GENE)
    f.text(490, y + 5, 'EN', size=14, anchor='middle', weight='bold')
    f.rect(570, y - 22, 130, 44, TSD, GENE)
    f.text(635, y + 5, 'RT', size=14, anchor='middle', weight='bold')
    f.text(575, y - 34, 'ORF2', size=14, anchor='middle', weight='bold')
    # Tail
    f.rect(730, y - 12, 62, 24, DARK, GENE, rx=3)
    f.text(761, y + 5, 'AAA(n)', size=12, fill='#FFFFFF', anchor='middle', weight='bold')
    # Regions
    f.line(100, y + 84, 290, y + 84, stroke=MUTED)
    f.text(195, y + 104, "5' UTR", size=13, fill=MUTED, anchor='middle')
    f.line(706, y + 84, 796, y + 84, stroke=MUTED)
    f.text(751, y + 104, "3' UTR", size=13, fill=MUTED, anchor='middle')
    f.text(490, y + 44, 'endonuclease', size=12, fill=MUTED, anchor='middle')
    f.text(635, y + 44, 'reverse transcriptase', size=12, fill=MUTED, anchor='middle')
    return f


def orientation_fig(same):
    f = Fig(W, 330)
    f.gene(50, x1=600, label=None, exons=((200, 270), (330, 420), (480, 580)))
    cy = 270
    f.chrom(cy, x1=600, label=None)
    gx = [205, 250, 355, 495, 540]
    sx = [230, 300, 370, 440, 510]
    if not same:
        sx = sx[::-1]
    for g, s in zip(gx, sx):
        f.line(g + 20, 112, s + 20, 218, stroke=MUTED, sw=1.2)
    for i, (g, s) in enumerate(zip(gx, sx)):
        f.read(g, 86, 'R2' if i % 2 else 'R1')
        f.read(s, 218, 'R1' if i % 2 else 'R2')
    f.text(24, 55, 'Parental gene', size=13, fill=MUTED, weight='bold')
    f.text(24, 104, 'Mates', size=13, fill=MUTED, weight='bold')
    f.text(24, 236, 'Cluster', size=13, fill=MUTED, weight='bold')
    f.text(24, 275, 'Insertion site', size=13, fill=MUTED, weight='bold')
    # Positions against positions
    ox, oy, s = 670, 220, 110
    f.line(ox, oy, ox + s + 10, oy, stroke=INK, arrow=True)
    f.line(ox, oy, ox, oy - s - 10, stroke=INK, arrow=True)
    for i in range(5):
        px = ox + 14 + i * 21
        py = oy - 14 - (i if same else 4 - i) * 21
        f.circle(px, py, 5, GENE)
    f.text(ox + s / 2, oy + 20, 'position at the gene', size=12, fill=MUTED, anchor='middle')
    f.text(ox - 8, oy - s - 20, 'position at the cluster', size=12, fill=MUTED)
    f.text(ox + s / 2, oy + 48, 'rho > 0: same strand' if same else 'rho < 0: opposite strands',
           size=14, fill=GENE, anchor='middle', weight='bold')
    return f


def breakpoints_fig():
    f = Fig(W, 400)
    # Reference
    y = 70
    f.text(24, y + 5, 'Reference', size=13, fill=MUTED, weight='bold')
    f.rect(130, y - 12, 300, 24, GREY, None)
    f.rect(430, y - 12, 60, 24, TSD, None, rx=0)
    f.rect(490, y - 12, 300, 24, GREY, None)
    f.text(280, y + 5, 'left flank', size=13, anchor='middle')
    f.text(460, y + 5, 'TSD', size=13, anchor='middle', weight='bold')
    f.text(640, y + 5, 'right flank', size=13, anchor='middle')
    f.line(431, y - 22, 431, y + 22, stroke=GENE, sw=2.5)
    f.text(431, y - 30, 'RBKP: first base of the right flank', size=12, fill=GENE, anchor='middle',
           weight='bold')
    f.line(489, y - 22, 489, y + 22, stroke=ACCENT, sw=2.5)
    f.text(489, y + 40, 'LBKP: last base of the left flank', size=12, fill=ACCENT, anchor='middle',
           weight='bold')
    # Sample
    y = 180
    f.text(24, y + 5, 'Sample', size=13, fill=MUTED, weight='bold')
    f.rect(110, y - 12, 180, 24, GREY, None)
    f.rect(290, y - 12, 60, 24, TSD, None, rx=0)
    f.rect(350, y - 12, 200, 24, EXON, GENE, rx=0)
    f.rect(550, y - 12, 50, 24, DARK, GENE, rx=0)
    f.rect(600, y - 12, 60, 24, TSD, None, rx=0)
    f.rect(660, y - 12, 150, 24, GREY, None)
    f.text(200, y + 5, 'left flank', size=13, anchor='middle')
    f.text(320, y + 5, 'TSD', size=13, anchor='middle', weight='bold')
    f.text(450, y + 5, 'retrocopy', size=13, anchor='middle')
    f.text(575, y + 5, 'AAAA', size=12, fill='#FFFFFF', anchor='middle', weight='bold')
    f.text(630, y + 5, 'TSD', size=13, anchor='middle', weight='bold')
    f.text(735, y + 5, 'right flank', size=13, anchor='middle')
    # Reads crossing the junctions
    for x, a, c in [(220, 130, 50), (250, 100, 80)]:
        row = 0 if a == 130 else 1
        ry = 224 + row * 20
        f.rect(x, ry, a, 12, MUTED, None, rx=3)
        f.rect(x + a, ry, c, 12, EXON, GENE, rx=3, sw=1.2, dash='4,3')
    for x, c, a in [(530, 70, 130), (560, 40, 100)]:
        row = 0 if c == 70 else 1
        ry = 224 + row * 20
        f.rect(x, ry, c, 12, EXON, GENE, rx=3, sw=1.2, dash='4,3')
        f.rect(x + c, ry, a, 12, MUTED, None, rx=3)
    f.text(220, 284, 'aligned part ends at LBKP; the clip', size=12, fill=INK)
    f.text(220, 300, 'is the beginning of the insertion', size=12, fill=INK)
    f.text(530, 284, 'aligned part starts at RBKP; the clip', size=12, fill=INK)
    f.text(530, 300, 'ends with the poly-A', size=12, fill=INK)
    # Legend
    ly = 340
    f.rect(110, ly, 26, 12, MUTED, None, rx=3)
    f.text(144, ly + 11, 'aligned to the reference', size=12)
    f.rect(320, ly, 26, 12, EXON, GENE, rx=3, sw=1.2, dash='4,3')
    f.text(354, ly + 11, 'soft clipped', size=12)
    f.rect(460, ly, 26, 12, TSD, None, rx=0)
    f.text(494, ly + 11, 'target site duplication: RBKP to LBKP', size=12)
    f.text(110, ly + 40, 'TSDLEN = LBKP - RBKP + 1. An insertion without duplication has '
           'LBKP = RBKP - 1 (TSDLEN = 0).', size=12, fill=MUTED)
    return f


def polya_fig():
    f = Fig(W, 330)

    def row(y, plus):
        f.text(24, y + 5, 'Retrocopy at +' if plus else 'Retrocopy at -', size=13, fill=MUTED,
               weight='bold')
        f.rect(150, y - 12, 150, 24, GREY, None)
        f.text(225, y + 5, 'left flank', size=13, anchor='middle')
        f.rect(650, y - 12, 150, 24, GREY, None)
        f.text(725, y + 5, 'right flank', size=13, anchor='middle')
        if plus:
            f.rect(300, y - 12, 280, 24, EXON, GENE, rx=0)
            f.rect(580, y - 12, 70, 24, DARK, GENE, rx=0)
            f.text(440, y + 5, "5'   mRNA   3'", size=13, anchor='middle')
            f.text(615, y + 5, 'AAAAA', size=12, fill='#FFFFFF', anchor='middle', weight='bold')
            f.rect(560, y + 24, 90, 12, EXON, GENE, rx=3, sw=1.2, dash='4,3')
            f.rect(650, y + 24, 110, 12, MUTED, None, rx=3)
            f.text(150, y + 62, 'The clip at RBKP ends with a run of A: The tail is at the right '
                   'junction', size=12)
        else:
            f.rect(300, y - 12, 70, 24, DARK, GENE, rx=0)
            f.rect(370, y - 12, 280, 24, EXON, GENE, rx=0)
            f.text(335, y + 5, 'TTTTT', size=12, fill='#FFFFFF', anchor='middle', weight='bold')
            f.text(510, y + 5, "3'   mRNA (reverse complement)   5'", size=13, anchor='middle')
            f.rect(190, y + 24, 110, 12, MUTED, None, rx=3)
            f.rect(300, y + 24, 90, 12, EXON, GENE, rx=3, sw=1.2, dash='4,3')
            f.text(150, y + 62, 'The clip at LBKP begins with a run of T: The tail is at the left '
                   'junction', size=12)

    row(60, True)
    row(190, False)
    f.text(150, 290, 'A clip with other bases before the run gives the exact tail length; a clip '
           'made only of A (or T)', size=12, fill=MUTED)
    f.text(150, 306, 'gives a lower bound (POLYAMIN). Sequences are read at the forward strand of '
           'the reference, as in SAM.', size=12, fill=MUTED)
    return f


def ends_fig():
    f = Fig(W, 460)
    exons = f.gene(60, exons=((200, 300), (380, 480), (560, 700)))
    f.text(340, 44, 'intron', size=12, fill=MUTED, anchor='middle')
    f.text(520, 44, 'intron', size=12, fill=MUTED, anchor='middle')
    for (a, b), n in zip(exons, ('exon 1', 'exon 2', 'exon 3')):
        f.text((a + b) / 2, 65, n, size=13, anchor='middle')
    p5, p3 = 240, 670
    f.line(p5, 30, p5, 96, stroke=ACCENT, sw=2.5)
    f.text(p5, 24, 'PG5P', size=12, fill=ACCENT, anchor='middle', weight='bold')
    f.line(p3, 30, p3, 96, stroke=ACCENT, sw=2.5)
    f.text(p3, 24, 'PG3P', size=12, fill=ACCENT, anchor='middle', weight='bold')

    # Insertion site and the clips
    y = 210
    f.text(24, y + 5, 'Insertion site', size=13, fill=MUTED, weight='bold')
    f.rect(150, y - 12, 150, 24, GREY, None)
    f.text(225, y + 5, 'left flank', size=13, anchor='middle')
    f.rect(300, y - 12, 100, 24, EXON, GENE, rx=3, sw=1.2, dash='4,3')
    f.text(372, y + 5, 'clip', size=13, anchor='middle')
    f.rect(300, y - 12, 44, 24, GENE, None, rx=3, op=0.35)
    f.rect(480, y - 12, 90, 24, EXON, GENE, rx=3, sw=1.2, dash='4,3')
    f.text(503, y + 5, 'clip', size=13, anchor='middle')
    f.rect(526, y - 12, 44, 24, GENE, None, rx=3, op=0.35)
    f.rect(570, y - 12, 60, 24, DARK, GENE, rx=0)
    f.text(600, y + 5, 'AAAA', size=12, fill='#FFFFFF', anchor='middle', weight='bold')
    f.rect(630, y - 12, 150, 24, GREY, None)
    f.text(705, y + 5, 'right flank', size=13, anchor='middle')
    f.link(322, y - 16, p5, 100, color=ACCENT)
    f.link(548, y - 16, p3, 100, color=ACCENT)
    f.text(150, y + 40, 'seed: the first 20 bases of the clips at LBKP', size=12)
    f.text(480, y + 40, 'seed: the last 20 bases before the poly-A', size=12)
    f.text(480, y + 56, 'of the clips at RBKP', size=12)

    # Predicted sequence
    y = 340
    f.text(24, y + 5, 'Predicted', size=13, fill=MUTED, weight='bold')
    x = 200
    for w, n in ((60, "exon 1'"), (100, 'exon 2'), (110, "exon 3'")):
        f.rect(x, y - 16, w, 32, EXON, GENE, rx=0)
        f.text(x + w / 2, y + 5, n, size=12, anchor='middle')
        x += w
    f.rect(x, y - 16, 70, 32, DARK, GENE, rx=0)
    f.text(x + 35, y + 5, 'poly-A', size=12, fill='#FFFFFF', anchor='middle', weight='bold')
    f.text(560, y - 2, 'reference bases of the exons from PG5P to PG3P,', size=12)
    f.text(560, y + 14, 'spliced, in the mRNA sense, plus the poly-A (SVLEN)', size=12)
    f.text(150, 410, 'The seeds are searched in both orientations. A forward match in a + gene '
           '(or a reverse one in a - gene)', size=12, fill=MUTED)
    f.text(150, 426, "means the insertion follows the mRNA sense at the reference. The junction "
           "with the poly-A is the 3' end.", size=12, fill=MUTED)
    return f


def intron_fig():
    sample_y, site_y = 340, 170
    seq_y, site_aln = sample_y - 56, site_y + 16
    f = Fig(W, sample_y + 60)
    exons = f.gene(REF_Y)
    # The retained intron, between the first and the second exons
    a, b = exons[0][1], exons[1][0]
    f.line(a + 5, REF_Y, b - 5, REF_Y, stroke=INTRON, sw=8)
    f.text((a + b) / 2, REF_Y - 14, 'intron', size=12, fill=INTRON, anchor='middle', weight='bold')
    f.chrom(site_y, label='Insertion site', sub='another chromosome')

    # The retrocopy kept the intron
    y = sample_y
    f.line(150, y, 730, y, stroke=CHR, sw=8, cap='round')
    x = 250
    f.rect(x - 22, y - 12, 22, 24, EXON, GENE, rx=3)
    f.rect(x, y - 20, 80, 40, EXON, GENE, rx=0)
    f.rect(x + 80, y - 20, 70, 40, INTRON, GENE, rx=0)
    f.text(x + 115, y + 5, 'intron', size=12, fill='#FFFFFF', anchor='middle', weight='bold')
    f.rect(x + 150, y - 20, 100, 40, EXON, GENE, rx=0)
    f.rect(x + 250, y - 20, 110, 40, EXON, GENE, rx=0)
    f.rect(x + 360, y - 12, 22, 24, EXON, GENE, rx=3)
    f.rect(x + 382, y - 12, 60, 24, DARK, GENE, rx=3)
    f.text(x + 412, y + 5, 'AAAA', size=12, fill='#FFFFFF', anchor='middle', weight='bold')
    f.text(754, y - 1, 'Retrocopy', weight='bold')
    f.text(754, y + 17, 'with a retained intron', size=12, fill=MUTED)
    f.text(24, REF_Y + 5, 'Aligned', size=13, fill=MUTED, weight='bold')
    f.text(24, sample_y + 5, 'Sequenced', size=13, fill=MUTED, weight='bold')

    # One mate at the flank, the other inside the intron
    f.read(165, seq_y, 'R1')
    f.read(345, seq_y, 'R2')
    f.pair(165, 345, seq_y)
    f.read(330, site_aln, 'R1')
    f.read(305, ALN_Y, 'R2')
    f.link(185, seq_y - 4, 350, site_aln + 32)
    f.link(365, seq_y - 4, 325, ALN_Y + 32)
    f.text(354, ALN_Y + 18, 'mate inside the intron', size=12, fill=INTRON, weight='bold')
    f.rect(300, site_y - 14, 100, 62, 'none', GENE, rx=12, dash='6,4')
    f.text(410, site_aln + 18, 'window of the retrocopy', size=12, fill=GENE)
    return f


FIGS = {
    'indistinguishable_alignment': lambda: alignment_fig('indist'),
    'abnormal_alignment_exon': lambda: alignment_fig('exon'),
    'abnormal_alignment_chr': lambda: alignment_fig('chr'),
    'abnormal_alignment_dist': lambda: alignment_fig('dist'),
    'abnormal_alignment_sr': lambda: alignment_fig('sr'),
    'abnormal_alignment_clustering': clustering_fig,
    'DBSCAN': dbscan_fig,
    'genotype': genotype_fig,
    'LINE1': line1_fig,
    'orientation_same_strand': lambda: orientation_fig(True),
    'orientation_opposite_strand': lambda: orientation_fig(False),
    'breakpoints': breakpoints_fig,
    'polya_polarity': polya_fig,
    'retrocopy_ends': ends_fig,
    'intron_retention': intron_fig,
}


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    for name, fn in FIGS.items():
        with open(os.path.join(out, name + '.svg'), 'w') as fp:
            fp.write(fn().svg())
    return 0


if __name__ == '__main__':
    sys.exit(main())
