.. _chap_output:

******
Output
******

The last step, ``make-vcf``, writes **three files** sharing the same
prefix (option ``-p``, *out* by default):

``<prefix>.vcf``
  One line per retrocopy insertion, in the `VCF
  <https://samtools.github.io/hts-specs/VCFv4.2.pdf>`_ format, with
  the annotation at the **INFO** column and the genotype of each
  sample at the **FORMAT** columns.

``<prefix>.rnames.tsv``
  The **reads** used to detect each retrocopy, by sample. So the user
  can go back to the alignment files and check the evidence by
  herself.

``<prefix>.insertions.fa``
  The **predicted sequence** of each retrocopy, in FASTA format.

The three files are tied together by the **retrocopy ID**
(``RTC000001``, ``RTC000002`` ...), which is the ``ID`` column of
the VCF, the first column of the read names file and the header of
each FASTA record. The IDs are assigned by ``merge-call``, so they
are stable for a given database, but **not** across different runs.

The VCF file
============

The header lists the version of **sideRETRO**, the reference genome
and its contigs (when ``-r`` is passed), the definition of every
INFO and FORMAT field, and finally the column names - with **one
column per sample**, named after the alignment file without its
extension.

Then, each line describes one retrocopy insertion:

+------------+-------------------------------------------------------------+
| Column     | Content                                                     |
+============+=============================================================+
| CHROM      | Chromosome of the insertion                                 |
+------------+-------------------------------------------------------------+
| POS        | The reference base **before** the insertion. When the two   |
|            | breakpoints are known, it is the base before the **TSD**,   |
|            | so the insertion is left aligned                            |
+------------+-------------------------------------------------------------+
| ID         | The retrocopy ID: ``RTC`` followed by six digits            |
+------------+-------------------------------------------------------------+
| REF        | The reference base at POS, or ``N`` without ``-r``          |
+------------+-------------------------------------------------------------+
| ALT        | ``<INS:ME:RTC>``: A **symbolic** allele for the insertion   |
|            | of a retrocopy. See the note below                          |
+------------+-------------------------------------------------------------+
| QUAL       | Phred scaled probability that the insertion is not real,    |
|            | from the genotype likelihoods, see :ref:`QUAL <qual>`       |
+------------+-------------------------------------------------------------+
| FILTER     | The checks of a canonical retrocopy that failed, separated  |
|            | by ``;``, or ``PASS``, see :ref:`FILTER <filter>`           |
+------------+-------------------------------------------------------------+
| INFO       | The annotation, see :ref:`INFO fields <info_fields>`        |
+------------+-------------------------------------------------------------+
| FORMAT     | ``GT:DP2:GL``, followed by one column per sample, see       |
|            | :ref:`FORMAT fields <format_fields>`                        |
+------------+-------------------------------------------------------------+

.. note::
   Why a **symbolic** ALT instead of the inserted sequence? Because
   the retrocopy sequence is **predicted** from the reference exons
   of the parental gene between the observed ends - it is not
   assembled from the reads. A sequence at the ALT column would be
   read by other tools as a precisely observed allele, which we
   cannot claim for the body of the insertion. The **observed** ends
   are at ``LEFT_SVINSSEQ`` and ``RIGHT_SVINSSEQ``, and the predicted
   sequence is at :ref:`the insertions file <insertions_file>`.

.. _qual:

QUAL
----

``QUAL`` is the Phred scaled probability that the site is **not**
variant: :math:`QUAL = -10 \log_{10} P(\text{all samples are } 0/0)`.
With a flat prior, each sample is homozygous reference with probability
:math:`L(0/0) / (L(0/0) + L(0/1) + L(1/1))`, from its genotype likelihoods
(``GL``). Samples without reads are not informative and are skipped; with
no informative sample, ``QUAL`` is ``.``. It is capped at 999.

.. note::
   ``QUAL`` inherits the :ref:`genotype model <method_genotype>`, which
   counts reads weighted by their mapping quality. It is a model based
   score, useful to **rank** the insertions, not a calibrated probability.

.. _filter:

FILTER
------

The ``FILTER`` column checks whether the insertion has the features of a
**canonical** retrocopy. It lists the checks that failed, separated by
``;``, or ``PASS`` when all of them passed. A check that **cannot be
evaluated** does not fail: Without both breakpoints, for instance, there is
no way to tell whether there is a TSD, so ``NoTSD`` is not raised. That is
what lets exome data or alignments without split reads (``bwa aln``) pass.

+-------------------+-------------------------------------------------------+
| LowQual           | ``QUAL`` below the threshold (``make-vcf -L``,        |
|                   | default 20)                                           |
+-------------------+-------------------------------------------------------+
| LowMAPQ           | Mean mapping quality of the supporting reads below    |
|                   | the threshold (``make-vcf -Q``, default 20). Typical  |
|                   | of reads in repeats or paralogs                       |
+-------------------+-------------------------------------------------------+
| AmbiguousParental | The parental gene is ambiguous: ``PGTYPE`` is not 1   |
+-------------------+-------------------------------------------------------+
| NoTSD             | Both breakpoints were found, but there is no target   |
|                   | site duplication (``TSDLEN <= 0``)                    |
+-------------------+-------------------------------------------------------+
| TSDOutlier        | ``TSDLEN`` out of the expected range: Greater than 50 |
|                   | or less than -20 - probably wrong breakpoints         |
+-------------------+-------------------------------------------------------+
| NoPolyA           | Both breakpoints were found, but there is no poly-A   |
|                   | tail at any of them                                   |
+-------------------+-------------------------------------------------------+
| NoJunction        | The parental gene has more than one exon and the      |
|                   | retrocopy should span a junction - its 5' end is not  |
|                   | known or not at the last exon - but no read spans it  |
+-------------------+-------------------------------------------------------+
| PolarityConflict  | The strand given by the poly-A disagrees with the one |
|                   | given by ``ORHO``                                     |
+-------------------+-------------------------------------------------------+

Any filter can be turned off with ``make-vcf -x`` (it may be passed more than
once), and the filters are applied by ``make-vcf``, so trying other thresholds
takes seconds - there is no need to run ``merge-call`` again. The high
confidence set is then::

   $ bcftools view -f PASS out.vcf

.. note::
   Failing a filter does not make an insertion false: It makes it **non
   canonical**. Some genuine retrocopies are inserted without a TSD, have lost
   the poly-A or come from a single exon parental gene. The filters point out
   what is missing, and ``HALLMARKS`` tells what was observed - the choice is
   up to you.

Precise and imprecise insertions
--------------------------------

The insertion point comes from the **clipped reads** crossing the
junctions between the genome and the retrocopy (see
:ref:`Breakpoints <method_breakpoints>`). When there is at least one
breakpoint, the retrocopy is **precise**: ``POS`` is the base before
the breakpoint and the ``SR`` field counts the supporting reads.

When there is no breakpoint at all - for instance, with low coverage
or when the aligner does not produce clipped alignments, as ``bwa
aln`` - the retrocopy is flagged as ``IMPRECISE``: ``POS`` is the
**mean** of the cluster window, whose limits relative to ``POS`` are
given by ``CIPOS``. All the annotation that depends on the junctions
(breakpoints, TSD, poly-A, retrocopy ends) is absent in that case.

.. _info_fields:

INFO fields
-----------

The fields are grouped by subject. Fields without a value are
**flags**: Their presence is the information.

**Site**

+------------+--------------------------------------------------------------+
| SVTYPE     | Always ``INS``                                               |
+------------+--------------------------------------------------------------+
| IMPRECISE  | Flag. No clipped read defines the insertion point, see above |
+------------+--------------------------------------------------------------+
| CIPOS      | Confidence interval around POS for imprecise insertions:     |
|            | The cluster window limits relative to POS                    |
+------------+--------------------------------------------------------------+
| DP         | Number of abnormal reads clustered at the insertion site,    |
|            | all samples together                                         |
+------------+--------------------------------------------------------------+
| SR         | Number of clipped reads supporting the breakpoints (left     |
|            | plus right). Only for precise insertions                     |
+------------+--------------------------------------------------------------+

**Breakpoints and target site duplication**

+----------------+----------------------------------------------------------+
| LBKP           | **Left breakpoint**: The last reference base of the left |
|                | flank, from the reads clipped at their right side        |
+----------------+----------------------------------------------------------+
| LSR            | Number of reads supporting LBKP                          |
+----------------+----------------------------------------------------------+
| RBKP           | **Right breakpoint**: The first reference base of the    |
|                | right flank, from the reads clipped at their left side   |
+----------------+----------------------------------------------------------+
| RSR            | Number of reads supporting RBKP                          |
+----------------+----------------------------------------------------------+
| TSDLEN         | Length of the **target site duplication**:               |
|                | ``LBKP - RBKP + 1``. Positive for a duplication, zero    |
|                | for a clean insertion and negative for a deletion at the |
|                | target site. Needs both breakpoints                      |
+----------------+----------------------------------------------------------+
| TSD            | The TSD sequence: The reference bases from RBKP to LBKP. |
|                | Needs ``-r`` and ``0 < TSDLEN <= 100``                   |
+----------------+----------------------------------------------------------+
| LEFT_SVINSSEQ  | **Observed** sequence of the insertion at the left       |
|                | junction: The consensus of the clipped reads at LBKP,    |
|                | read from the junction into the insertion. A base where  |
|                | the reads tie is ``N``                                   |
+----------------+----------------------------------------------------------+
| RIGHT_SVINSSEQ | Observed sequence at the right junction: The consensus   |
|                | of the clipped reads at RBKP, ending at the junction.    |
|                | A base where the reads tie is ``N``                      |
+----------------+----------------------------------------------------------+

A retrocopy with a TSD looks like this at the reference and at the
sequenced genome::

   reference:  ...left flank | TSD |                        right flank...
   genome:     ...left flank | TSD | retrocopy + poly-A | TSD | right flank...
                                  ^                          ^
                                LBKP                        RBKP

The reads from the left flank cross the junction after the **last**
base of the TSD (``LBKP``), while the reads from the right flank
cross it before the **first** base of the TSD (``RBKP``). That is
why ``RBKP`` comes **before** ``LBKP`` at the reference when there
is a duplication.

**Poly-A tail and polarity**

+-------------+--------------------------------------------------------------+
| POLYA       | Length of the poly-A tail observed in the reads crossing     |
|             | its junction. Genomic A bases adjacent to the tail cannot    |
|             | be told apart from it, so the run is reported as seen        |
+-------------+--------------------------------------------------------------+
| POLYAMIN    | Flag. POLYA is a **lower bound**: The tail extends beyond    |
|             | the longest clipped read                                     |
+-------------+--------------------------------------------------------------+
| POLYASR     | Number of reads supporting the poly-A                        |
+-------------+--------------------------------------------------------------+
| POLYASTRAND | Strand of the retrocopy given by the poly-A junction         |
|             | (``+``/``-``): The tail is at the right junction for ``+``   |
|             | and, reverse complemented, at the left junction for ``-``.   |
|             | Independent of POLARITY: If they disagree, the FILTER says   |
|             | ``PolarityConflict``                                         |
+-------------+--------------------------------------------------------------+
| POLARITY    | Orientation of the event at the reference (``+``/``-``),     |
|             | from ORHO                                                    |
+-------------+--------------------------------------------------------------+
| ORHO        | Spearman's rho between the positions of the reads at the     |
|             | parental gene and of their mates at the insertion site, see  |
|             | :ref:`Orientation <method_orientation>`. Only when the test  |
|             | is significant (``make-vcf -e``)                             |
+-------------+--------------------------------------------------------------+

**Hallmarks**

+------------+--------------------------------------------------------------+
| HALLMARKS  | The retrocopy hallmarks observed, separated by comma: ``TSD``|
|            | (target site duplication), ``POLYA`` (poly-A tail) and       |
|            | ``EEJ`` (reads spanning exon-exon junctions). With the       |
|            | filters, it tells a hallmark that was looked for and is      |
|            | missing (e.g. ``NoTSD``) from one that could not be looked   |
|            | for (neither ``TSD`` here nor ``NoTSD`` at ``FILTER``)       |
+------------+--------------------------------------------------------------+

**Parental gene**

+------------+--------------------------------------------------------------+
| PG         | Parental gene name(s), separated by ``/``                    |
+------------+--------------------------------------------------------------+
| PGTYPE     | How the parental gene was resolved: ``1`` single parental    |
|            | gene; ``2`` overlapped parental genes; ``4`` near parental   |
|            | genes; ``8`` hotspot - multiple parental genes with          |
|            | retrocopies at the same segment. Values may be summed        |
+------------+--------------------------------------------------------------+
| PGEXONS    | Number of exons of each parental gene, as ``gene:n``. A      |
|            | **single exon** gene has no exon-exon junction, so a missing |
|            | EEJ is expected for it                                       |
+------------+--------------------------------------------------------------+
| EEJ        | Number of reads spanning **exon-exon junctions** of the      |
|            | parental gene: Direct evidence of a processed transcript     |
+------------+--------------------------------------------------------------+
| EEJLIST    | The junctions spanned, as ``exon_id1-exon_id2`` separated by |
|            | ``/``. With alternative exons, it tells which isoform was    |
|            | retrocopied                                                  |
+------------+--------------------------------------------------------------+

.. note::
   The exon ids at ``EEJLIST``, ``PG5P`` and ``PG3P`` come from the annotation.
   When the annotation has no exon ids - a BED6 file, for instance - they are the
   exon coordinates, ``chr:start_end`` - e.g.
   ``EEJLIST=chr2:1001_1400-chr2:2001_2400`` - see :ref:`Annotation formats <annotation>`.

**Retrocopy ends and length** (need ``merge-call -r``)

+------------+--------------------------------------------------------------+
| PG5P       | The **5' end** of the retrocopy into the parental gene, as   |
|            | ``exon_id:position:reads``. The position is at the reference |
|            | genome. A 5' end inside the transcript means a **5'**        |
|            | **truncated** retrocopy                                      |
+------------+--------------------------------------------------------------+
| PG3P       | The **3' end**, before the poly-A, as                        |
|            | ``exon_id:position:reads``                                   |
+------------+--------------------------------------------------------------+
| SVLEN      | **Predicted** length of the insertion: The reference exon    |
|            | bases between the ends plus the poly-A. When PG3P is not     |
|            | observed, the 3' end is assumed at the end of the last exon  |
+------------+--------------------------------------------------------------+

**Genomic context**

+------------+--------------------------------------------------------------+
| EXONIC     | Host gene(s) when the insertion falls inside an exon         |
+------------+--------------------------------------------------------------+
| INTRONIC   | Host gene(s) when the insertion falls inside a gene, but not |
|            | at this gene exons                                           |
+------------+--------------------------------------------------------------+
| NEAR       | Gene(s) near the insertion, within ``make-vcf -n`` bases,    |
|            | for intergenic insertions                                    |
+------------+--------------------------------------------------------------+

.. _format_fields:

FORMAT fields
-------------

One column per sample, with:

+------------+--------------------------------------------------------------+
| GT         | Genotype: ``0/0`` homozygous reference, ``0/1``              |
|            | heterozygous, ``1/1`` homozygous alternate or ``./.`` when   |
|            | the likelihoods cannot decide                                |
+------------+--------------------------------------------------------------+
| DP2        | Read depth at this sample: Reference depth and alternate     |
|            | depth, separated by comma                                    |
+------------+--------------------------------------------------------------+
| GL         | Genotype likelihoods (log10) for ``0/0``, ``0/1`` and        |
|            | ``1/1``, see :ref:`Genotype <method_genotype>`               |
+------------+--------------------------------------------------------------+

An example
----------

The line below comes from the :ref:`simulated data <sim_test>`
shipped with **sideRETRO**: A retrocopy of the gene *PARENT*,
inserted at *chr1* with a 12 bases TSD and a 25 bases poly-A tail.
The INFO column was wrapped for reading::

   chr1  9999  RTC000001  A  <INS:ME:RTC>  999.0  PASS
     SVTYPE=INS;ORHO=0.514382;POLARITY=+;PG=PARENT;PGTYPE=1;PGEXONS=PARENT:3;
     HALLMARKS=TSD,POLYA,EEJ;EEJ=6;EEJLIST=ENSEP0-ENSEP1;PG5P=ENSEP0:1001:8;PG3P=ENSEP2:3398:3;SVLEN=1225;
     LEFT_SVINSSEQ=CACTGATACAGCAGACTCATTCCCCTATAAGCGAATCGTA;
     RIGHT_SVINSSEQ=ACTCAGCAGGACTCGGTTTAAAAAAAAAAAAAAAAAAAAAAAAAAA;
     DP=42;SR=16;LBKP=10011;LSR=8;RBKP=10000;RSR=8;TSDLEN=12;TSD=GAACGCCCCGTG;
     POLYA=27;POLYASR=8;POLYASTRAND=+
     GT:DP2:GL  1/1:0,42:-253,-12.6,-1.8e-05

Reading it: The insertion is **precise** (no ``IMPRECISE``), with 8
reads at each breakpoint (``LSR``, ``RSR``) and a TSD of 12 bases
(``TSDLEN``, ``TSD``). The poly-A has 27 bases - 25 simulated plus
two A at the end of the last exon - and puts the retrocopy at the
``+`` strand. Six reads span the junction between the first two exons
(``EEJ``, ``EEJLIST``), the retrocopy is **full length** (``PG5P`` at
the first base of exon *ENSEP0*) and its predicted length is 1225
bases: Three exons of 400 bases plus the poly-A. The single sample is
homozygous for the insertion, hence the maximum ``QUAL``, and the insertion
has all the hallmarks of a retrocopy (``HALLMARKS``), so it passes all the
filters (``PASS``).

.. _rnames_file:

The read names file
===================

A tab separated file with **one read per line**: The retrocopy ID, the
sample, the read name and its **role** in the evidence, after a header
line::

   #retrocopy_id	sample	read_name	role
   RTC000001	NA12878	A00123:45:HV7KJDSXY:1:1101:12345:6789	discordant
   RTC000001	NA12878	A00123:45:HV7KJDSXY:1:1102:2210:1003	split_left
   RTC000001	NA12878	A00123:45:HV7KJDSXY:1:1203:7781:2004	split_right,polya
   RTC000001	NA12778	A00456:12:HT2LMDSXY:2:2201:8812:4410	junction

Samples are named as in the VCF header. The file follows the VCF order and,
for each retrocopy and sample, the reads are sorted by name. These are the
**abnormal reads** clustered at the insertion site - exactly the evidence that
called the retrocopy. Each read name stands for the **fragment** (the read and
its mate), and its roles are:

+-------------+-------------------------------------------------------------+
| split_left  | Clipped at the accepted **left breakpoint** (counted at     |
|             | ``LSR``)                                                    |
+-------------+-------------------------------------------------------------+
| split_right | Clipped at the accepted **right breakpoint** (counted at    |
|             | ``RSR``)                                                    |
+-------------+-------------------------------------------------------------+
| polya       | Its clip has the **poly-A** tail (counted at ``POLYASR``)   |
+-------------+-------------------------------------------------------------+
| junction    | Spans an **exon-exon junction** of the parental gene        |
|             | (counted at ``EEJ``)                                        |
+-------------+-------------------------------------------------------------+
| discordant  | None of the above: A read at the insertion site whose mate  |
|             | falls into the parental gene                                |
+-------------+-------------------------------------------------------------+

A read may have more than one role, separated by comma - a read clipped at
the right breakpoint whose clip ends with the poly-A is ``split_right,polya``.
Only the reads supporting the **accepted** breakpoints are split reads: A read
clipped somewhere else, by an adapter for instance, is ``discordant``.

A read name never has a tab (see the `SAM specification
<https://samtools.github.io/hts-specs/SAMv1.pdf>`_), so the names are written
as they are - with the colons of Illumina, for instance - and the file can be
read by any tool, with no escaping. To inspect the reads crossing the left
breakpoint of a retrocopy in a sample::

   $ awk -F'\t' '$1 == "RTC000001" && $2 == "NA12878" && $4 ~ /split_left/ {print $3}' \
       out.rnames.tsv > reads.txt
   $ samtools view -N reads.txt NA12878.bam

In R or Python, read it as any table and group it by ``retrocopy_id`` and
``sample``::

   >>> import pandas as pd
   >>> rnames = pd.read_csv('out.rnames.tsv', sep='\t')
   >>> rnames.groupby(['#retrocopy_id', 'sample', 'role']).size()

.. _insertions_file:

The predicted insertions file
=============================

A FASTA file with one record per retrocopy whose **5' end** was found
(so it needs ``merge-call -r``)::

   >RTC000001 chr1:9999 PG=PARENT SVLEN=1225 predicted
   CACTGATACAGCAGACTCATTCCCCTATAAGCGAATCGTAAATTATAACCGTCATGAGGA
   ...

The header carries the ID, the insertion site, the parental gene and
the length. The sequence is written in the **mRNA sense** (5' to
3'), whatever the strand of the insertion, and it is made of:

- the reference bases of the parental exons from the 5' end (PG5P)
  to the 3' end (PG3P, or the end of the last exon when not
  observed), skipping the introns;
- the poly-A tail, with the observed length (POLYA).

The word *predicted* in the header is deliberate: The **junctions**
and the **poly-A** are observed in the reads, the **body** between
them is taken from the reference. It is right whenever the retrocopy
copied the annotated exons faithfully, and it says nothing about
variants carried from the parental allele, mutations acquired since
the insertion or unannotated isoforms. Use it as a template, not as
an assembly.
