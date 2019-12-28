.. _chap_result:

*******
Results
*******

To assess **sideRETRO** performance and accuracy, we made a simulated
dataset whereupon we could have control of **true** positive and
negative values.

Dataset
========

Our dataset for testing is composed of 5 simulated human whole-genome
sequencing with 40x of depth and 10 randomly distributed retrocopies
each. In total, we had a list of 30 retrocopies consisting of the last
1000 bases of the largest transcript of the parental gene. All retrocopies
was stochastically designed for chromosome, position, strand and zygosity.

.. table:: 30 parental gene names chosen for the simulation
   :widths: auto

   +--------+--------+-------+---------+-------+
   | ACAP3  | AMY2B  | ATM   | BCL2    | BRAF  |
   +--------+--------+-------+---------+-------+
   | BRCA1  | BRIP1  | DVL1  | FANCD2  | FAT1  |
   +--------+--------+-------+---------+-------+
   | GNB1   | KLHL17 | MECP2 | MUTYH   | NPHP4 |
   +--------+--------+-------+---------+-------+
   | OR4F5  | PBX1   | PRCC  | PTEN    | RER1  |
   +--------+--------+-------+---------+-------+
   | SAMD11 | SET    | SIRT1 | SMARCE1 | SOX2  |
   +--------+--------+-------+---------+-------+
   | TMEM52 | TP53   | TPM3  | USP8    | XPO1  |
   +--------+--------+-------+---------+-------+

For each individual human whole-genome, we **sampled with replacement**
10 retrocopies from our list, so that some events were shared among
part of our subjects.

.. table:: All randomly designed retrocopies for simulation. \*HE heterozygous
   and \*HO homozygous alternate
   :widths: auto

   +---------------+------------+-----------+--------+----------+
   | Parental gene | Chromosome | Position  | Strand | Zygosity |
   +===============+============+===========+========+==========+
   | AMY2B         | chr5       | 122895832 |  \-    | HE       |
   +---------------+------------+-----------+--------+----------+
   | ACAP3         | chr11      | 133246346 |  \+    | H0       |
   +---------------+------------+-----------+--------+----------+
   | ATM           | chr19      | 16443740  |  \+    | HO       |
   +---------------+------------+-----------+--------+----------+
   | BCL2          | chr4       | 146453786 |  \+    | HO       |
   +---------------+------------+-----------+--------+----------+
   | BRAF          | chr18      | 22375812  |  \-    | HO       |
   +---------------+------------+-----------+--------+----------+
   | BRCA1         | chr7       | 102911369 |  \-    | HO       |
   +---------------+------------+-----------+--------+----------+
   | BRIP1         | chr12      | 113357216 |  \-    | HO       |
   +---------------+------------+-----------+--------+----------+
   | DVL1          | chr5       | 54105196  |  \-    | HE       |
   +---------------+------------+-----------+--------+----------+
   | FANCD2        | chr3       | 121339674 |  \+    | HE       |
   +---------------+------------+-----------+--------+----------+
   | FAT1          | chr16      | 65659952  |  \-    | HE       |
   +---------------+------------+-----------+--------+----------+
   | GNB1          | chr12      | 70734022  |  \+    | HE       |
   +---------------+------------+-----------+--------+----------+
   | MUTYH         | chr4       | 2716745   |  \+    | HO       |
   +---------------+------------+-----------+--------+----------+
   | PBX1          | chr21      | 22718521  |  \-    | HE       |
   +---------------+------------+-----------+--------+----------+
   | PRCC          | chr1       | 190777903 |  \-    | HO       |
   +---------------+------------+-----------+--------+----------+
   | PTEN          | chr2       | 140252560 |  \-    | HE       |
   +---------------+------------+-----------+--------+----------+
   | RER1          | chr13      | 55179109  |  \+    | HE       |
   +---------------+------------+-----------+--------+----------+
   | SAMD11        | chr4       | 14585135  |  \+    | HO       |
   +---------------+------------+-----------+--------+----------+
   | SET           | chr7       | 154178578 |  \-    | HO       |
   +---------------+------------+-----------+--------+----------+
   | SIRT1         | chrX       | 120716688 |  \-    | HO       |
   +---------------+------------+-----------+--------+----------+
   | SOX2          | chr10      | 88689163  |  \-    | HE       |
   +---------------+------------+-----------+--------+----------+
   | TMEM52        | chr1       | 82897536  |  \+    | HO       |
   +---------------+------------+-----------+--------+----------+
   | TP53          | chr5       | 42938944  |  \-    | HO       |
   +---------------+------------+-----------+--------+----------+
   | TPM3          | chr7       | 3488208   |  \-    | HE       |
   +---------------+------------+-----------+--------+----------+
   | USP8          | chr7       | 41333317  |  \+    | HO       |
   +---------------+------------+-----------+--------+----------+
   | XPO1          | chr10      | 33172062  |  \-    | HE       |
   +---------------+------------+-----------+--------+----------+

.. table:: Parental gene name for the retrocopies sampled by individual whole-genome
   :widths: auto

   +-------+---------+------+--------+---------+
   | IND1  | IND2    | IND3 | IND4   | IND5    |
   +=======+=========+======+========+=========+
   | ATM   | AMY2B   | ATM  | AMY2B  | AMY2B   |
   +-------+---------+------+--------+---------+
   | BCL2  | BRCA1   | GNB1 | BRAF   | BRIP1   |
   +-------+---------+------+--------+---------+
   | BRIP  | FANCD2  | MECP | BRIP1  | FANCD2  |
   +-------+---------+------+--------+---------+
   | DVL1  | FAT1    | PBX1 | DVL1   | NPHP4   |
   +-------+---------+------+--------+---------+
   | MECP  | KLHL17  | PRCC | FANCD2 | PTEN    |
   +-------+---------+------+--------+---------+
   | OR4F  | MUTYH   | PTEN | GNB1   | RER1    |
   +-------+---------+------+--------+---------+
   | PTEN  | RER1    | RER1 | KLHL17 | SMARCE1 |
   +-------+---------+------+--------+---------+
   | SET   | SET     | SAMD | MECP2  | TP53    |
   +-------+---------+------+--------+---------+
   | SOX2  | SMARCE1 | SIRT | SET    | TPM3    |
   +-------+---------+------+--------+---------+
   | XPO1  | SOX2    | XPO1 | TMEM52 | XPO1    |
   +-------+---------+------+--------+---------+

Simulation
==========

We used the **SANDY** tool (version v0.23), *A straightforward and complete next-generation
sequencing read simulator* [1]_, for simulate all 5 genomes according to the structural
variations that we desgined and according to the sampling. SANDY demands 2 steps for the
task: First we indexed all retrocopies for each individual and next we could simulate all
whole-genome sequencing.

.. code-block:: sh

   # Reference genome
   REF_FASTA=hg38.fa

   # Coverage depth
   DEPTH=40

   # Retrocopies by individual
   IND=(ind1.txt ind2.txt ind3.txt ind4.txt ind5.txt)

   # Index all retrocopies
   for ind in "${IND[@]}"; do
     sandy variation add --structural-variation=${ind%%.txt} $ind
   done

   # Simulate all genomes
   for ind in "${IND[@]}"; do
     sandy_index=${ind%%.txt}
     sandy genome \
       --id='%i.%U_%c:%S-%E_%v' \
       --structural-variation=$sandy_index \
       --output-dir=$sandy_index \
       --jobs=20 \
       --seed=42 \
       --quality-profile='hiseq_101' \
       --coverage=$DEPTH \
       --verbose \
       $REF_FASTA
   done

Running sideRETRO
=================

After our simulated dataset was ready, we could test the sideRETRO capabilities
to detect the designed somatic retrocopies.

.. code-block:: sh

   # Our simulated BAM files list
   LIST=(ind1/out.bam ind2/out.bam ind3/out.bam ind4/out.bam ind5/out.bam)

   # GENCODE annotation v31
   ANNOTATION=gencode.v31.annotation.gff3.gz

   # GENCODE reference genome
   REF_FASTA=hg38.fa

   # Run process-sample step
   sider process-sample \
     --cache-size=20000000 \
     --output-dir=result \
     --threads=5 \
     --max-distance=15000 \
     --alignment-frac=0.9 \
     --phred-quality=20 \
     --sorted \
     --log-file=ps.log \
     --annotation-file=$ANNOTATION \
     "${LIST[@]}"

   # Run merge-call step
   sider merge-call \
     --cache-size=20000000 \
     --epsilon=500 \
     --min-pts=20 \
     --genotype-support=5 \
     --near-gene-rank=3 \
     --log-file=mc.log \
     --threads=10 \
     --phred-quality=20 \
     --in-place \
     result/out.db

   # Finally run make-vcf
   sider make-vcf --reference-file=$REF_FASTA result/out.db

Analysis
========

References and Further Reading
==============================

.. [1] Miller, Thiago et al. (2019).
   galantelab/sandy: Release v0.23 (Version v0.23).
   Zenodo. http://doi.org/10.5281/zenodo.2589575.
