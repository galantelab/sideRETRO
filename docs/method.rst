.. _chap_methodology:

***********
Methodology
***********

The sideRETRO **methodology** consists of detecting
abnormal (discordant) alignments in `SAM/BAM
<https://samtools.github.io/hts-specs/SAMv1.pdf>`_
file and, with an **unsupervised mashine learning**
algorithm, clustering these reads and genotype in order
to dicover somatic retrocopy insertions. Care is taken
to ensure the quality and consistency of the data,
taking into acount the features that characterize a
retrocopy mobilization, such as the absence of
**intronic** and **regulatory** regions.

.. note:: For more detail about the jargon, see `Retrocopy in a nutshell <retrocopy.rst>`_

What do you mean when you say "abnormal alignments"?
====================================================

When a structural variation, such as a retrotransposition,
occurs into an individual and her genome is sequenced with
a next-generation sequencing technology (e.g. `illumina
<https://www.illumina.com/>`_), we may **expect** that the
aligner (e.g. `BWA <http://bio-bwa.sourceforge.net/>`_,
`Bowtie <http://bowtie-bio.sourceforge.net/index.shtml>`_)
will be **confused** as to the origin of certain reads. As the
retrocopies come from a **mature mRNA**, reads from the
retrocopy may be **erroneously** aligned to an **exon** of the
**parental gene**:

.. image:: images/indistinguishable_alignment.png
   :scale: 25%
   :align: center

These kind of alignment may be called **indistinguishable**,
because they do not give any **clue** about the presence of
the retrocopy. However, for our luck, there are reads
with abnormal (discordant) alignments which could be
helpful according to their characteristics:

- Paired-end reads aligned at **different** exons
- Paired-end reads aligned at **different** chromosomes
- Paired-end reads aligned at **distant** regions
- Splitted reads (Reads with **supplementary** alignment)

We will talk about each one as best as we can in the
next lines.

Alignment at different exons
============================

When paired-end reads are mapped to contiguous exons and they
came from a genomic sequencing - which of course is not
expected.

.. image:: images/abnormal_alignment_exon.png
   :scale: 25%
   :align: center

This kind of alignment is useful for **assume** a
retrotransposition for the given parental gene, however it
is not possible to annotate the **genomic position** of the event.

Alignment at different chromosomes
==================================

When the retrotransposition does not occur into the **same** parental
gene chromosome, it may happen that one read come from a **near**
intergenic region and its pair from the somatic **retrocopy**. As the
retrocopy **does not exist** in the reference genome, the aligner will
**map** one read to the retrotransposition chromosome and its pair to
the parental gene **exon**.

.. image:: images/abnormal_alignment_chr.png
   :scale: 25%
   :align: center

This alignment is useful to **estimate** the genomic position of the
event, but not with so much **precision** concerning to the **insertion
point**.

Alignment at distant regions
============================

If a retrocopy is inserted into the **same** chromosome of its parental
gene, possibly it will occur at a **distant** location. As well as
*"alignment at different chromosomes"*, one read may come from a near
intergenic region and its mate from the somatic retrocopy. So when
the aligner try to map these reads, we will observe that one fall
inside the parental gene exon, while its pair is mapped to a **distant
region**.

.. image:: images/abnormal_alignment_dist.png
   :scale: 25%
   :align: center

Splitted reads
==============

The most **important** kind of alignment when detecting structural variations.
The splitted read may occur when **part** of the **same** read come from a near
intergenic region and part from the somatic retrocopy. When the aligner
**try** to map the read, it will need to **create** another one to represent
the splitted part, which is called **supplementary**.

.. image:: images/abnormal_alignment_sr.png
   :scale: 25%
   :align: center

This alignment is useful to detect the **insertion point** with a
**good precision**.

So far we can resume all abnormal alignments according to their power
to detect the retrotransposition coordinate and its exact insertion
point:

+--------------------------+------------+-----------------+
| Abnormal alignments      | Coordinate | Insertion point |
+==========================+============+=================+
| At different exons       |     NO     |      NO         |
+--------------------------+------------+-----------------+
| At different chromosomes |     YES    |      NO         |
+--------------------------+------------+-----------------+
| At distant regions       |     YES    |      NO         |
+--------------------------+------------+-----------------+
| Splitted read            |     YES    |      YES        |
+--------------------------+------------+-----------------+

sideRETRO uses **only** the abnormal alignments **capable** to
detect **at least** the coordinate, so those that fall into
*different exons* are dismissed.

