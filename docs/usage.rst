.. _chap_usage:

***************
Using sideRETRO
***************

General Syntax
==============

**sideRETRO** compiles to an executable called ``sider``, with three commands
run one after the other:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Command
     - What it does
   * - ``process-sample`` (``ps``)
     - Reads the **SAM/BAM/CRAM** files and stores their abnormal alignments
       at a SQLite3 database
   * - ``merge-call`` (``mc``)
     - Merges the databases, clusters the abnormal alignments, annotates the
       retrocopy insertions and genotypes them at each sample
   * - ``make-vcf`` (``vcf``)
     - Writes the retrocopies to a **VCF** and to its companion files: The
       supporting reads, the predicted sequences and the sequences for primer
       design

Called without any argument, ``sider`` shows its usage, and so does each
command. ``sider -V`` shows the installed version::

  $ sider -V
  sideRETRO 1.14.0

The full help is also at the `readthedocs page <https://sideretro.readthedocs.io>`_
and at ``sider --help`` or ``sider <command> --help``. See :ref:`A Practical
Workflow <pract_wf>` and :ref:`Running with Docker <run_dck>` for complete
examples.

.. _common_options:

Common options
--------------

The three commands share these options:

.. list-table::
   :header-rows: 1
   :widths: 25 60 15

   * - Option
     - Description
     - Default
   * - ``-h``, ``--help``
     - Show the help options
     -
   * - ``-q``, ``--quiet``, ``--silent``
     - Print only the errors - or nothing at the terminal, when a log file
       is given
     -
   * - ``-d``, ``--debug``
     - Increase the verbosity to the debug level
     -
   * - ``-l``, ``--log-file``
     - Print the log messages to a file as well
     -
   * - ``-o``, ``--output-dir``
     - Output directory, created if it does not exist
     - ``.``
   * - ``-p``, ``--prefix``
     - Prefix of the output files
     - ``out``

Each command ends with a **summary** of what it did - at the terminal and at the
log file (``-l``), unless ``-q`` is passed - and the command to run next. For
example, ``merge-call`` reports how many retrocopies were found and how many of
them have each hallmark::

  ======================== merge-call summary ========================

   Input
     Databases                1
     Samples                  1
     Abnormal alignments      98
     Clusters                 1 (after all filters)

   Retrocopies
     Total                    1
     Single parental gene     1 (100.0%)
     Precise                  1 (100.0%)
     Both breakpoints         1 (100.0%)
     Target site duplication  1 (100.0%)
     Poly-A tail              1 (100.0%)
     Exon-exon junctions      1 (100.0%)
     Predicted sequence       1 (100.0%)

   Retrocopies by sample (with alternate reads)
     sim                      1

   Output
     Database                 out/merged.db
     Next step                sider make-vcf -r ref.fa out/merged.db

     Elapsed time             1s
  ====================================================================

``process-sample`` reports the exons indexed and, for each alignment file, the
alignments read, the abnormal ones and how many of them overlap exons; and
``make-vcf`` reports the records written, how many of them ``PASS``, and the
output files.


Command ``process-sample``
==========================

``process-sample`` (``ps``) is the first step: It reads one or more alignment
files and keeps, at a SQLite3 database, the **abnormal** fragments with an
alignment at an **exon** of the annotation - the evidence that the next steps
work on. See the :ref:`methodology <chap_methodology>` for what makes a fragment
abnormal.

Its arguments are the alignment files, in SAM/BAM/CRAM format. Besides the
:ref:`common options <common_options>`:

.. list-table:: Mandatory options
   :header-rows: 1
   :widths: 25 60 15

   * - Option
     - Description
     - Default
   * - ``-a``, ``--annotation-file``
     - Exons of the protein coding genes, in GTF/GFF3 or BED6 format. See
       :ref:`Annotation formats <annotation>`
     -
   * - ``-i``, ``--input-file``
     - File with a newline separated list of alignment files. Not needed if
       the files are passed as arguments; if both are given, their union is
       used
     -

.. list-table:: Annotation options (GTF/GFF3 only)
   :header-rows: 1
   :widths: 25 60 15

   * - Option
     - Description
     - Default
   * - ``-T``, ``--gff-feature``
     - The feature (3rd column) of the exons
     - ``exon``
   * - ``-A``, ``--gff-attribute``
     - An attribute (9th column) to filter in, as ``key=value``.
       Alternatives separated by ``|`` must match any (OR); when passed
       multiple times, all of them must match (AND). Values match exactly
     - ``transcript_type=protein_coding``
   * - ``-N``, ``--gff-gene-name``
     - The attribute with the gene name. If missing, ``gene_id`` is used
     - ``gene_name``
   * - ``-E``, ``--gff-exon-id``
     - The attribute with the exon id. If missing, the exon is named
       ``chr:start_end``
     - ``exon_id``

.. list-table:: Read quality options
   :header-rows: 1
   :widths: 25 60 15

   * - Option
     - Description
     - Default
   * - ``-Q``, ``--phred-quality``
     - Minimum mapping quality (MAPQ) of the reads. A failing alignment drops
       its fragment, but a supplementary alignment below it is only left out
     - 20
   * - ``-M``, ``--max-base-freq``
     - Maximum frequency of a single base in a read: Above it, a read is of
       low complexity, such as poly-A
     - 0.75
   * - ``-D``, ``--deduplicate``
     - Remove the duplicated fragments: Reads and mates sharing the same 5'
       positions
     -

.. list-table:: Processing options
   :header-rows: 1
   :widths: 25 60 15

   * - Option
     - Description
     - Default
   * - ``-s``, ``--sorted``
     - Assume the reads are grouped by read name, even without the header tag
       ``SO:queryname``. Grouped files are read once, the others twice
     -
   * - ``-t``, ``--threads``
     - Number of threads. Each file uses one, and the spare ones decompress
       BAM/CRAM
     - 1
   * - ``-m``, ``--max-distance``
     - Maximum distance between the reads of a pair. Farther pairs are
       abnormal
     - 10000
   * - ``-F``, ``--alignment-frac``
     - Minimum overlap between an alignment and an exon, as a fraction of the
       **alignment**
     - 0.5
   * - ``-f``, ``--exon-frac``
     - Minimum overlap between an alignment and an exon, as a fraction of the
       **exon**. Beware: An exon longer than the read length divided by it
       loses all its reads, as the last exons with long 3' UTRs. A warning is
       given when it is set; prefer ``-F``
     - 1 base
   * - ``-e``, ``--either``
     - Either ``-f`` or ``-F`` must be satisfied. Without it, both must be
     -
   * - ``-r``, ``--reciprocal``
     - The overlap must be reciprocal: ``-F`` takes the value of ``-f``
     -
   * - ``-c``, ``--cache-size``
     - SQLite3 cache size in KiB
     - 200000

So, with three files - *f1.bam*, *f2.bam* and *f3.sam* - type::

  $ sider process-sample -a annotation_file.gtf f1.bam f2.bam f3.sam

or give them at a list (e.g. *my_files_list.txt*) with ``-i``::

  $ sider process-sample -a annotation_file.gtf -i my_files_list.txt

Both write a single database, *out.db* (the prefix changes with ``-p``): The
abnormal reads of all the files go to one table, **with their sequences**, which
are needed to annotate the junctions later. For one database per file, run one
``process-sample`` per file - ``merge-call`` merges them afterwards.

.. note::
   The alignment files do not need to be sorted. When they are **not**
   grouped by read name (``SO:queryname``), each file is read **twice**
   in order to catch all the alignments of each abnormal fragment. If
   there are more threads (``-t``) than files, the spare ones are used to
   decompress the BAM/CRAM.

A fuller example, with the short name of the command::

  $ sider ps \
      -i my_files_list.txt \
      -a annotation_file.gtf \
      -o output_dir \
      -p my_reads_database \
      -l my_log_file.log \
      -c 2000000 \
      -F 0.9 \
      -t 3

It writes *output_dir/my_reads_database.db* and the log to *my_log_file.log*,
with a cache of 2 GB, 3 threads to process the files in parallel and a stricter
overlap with the exons: 90% of each alignment (``-F``) instead of half.

To see ``process-sample`` chained in a real workflow, see :ref:`A Practical
Workflow <pract_wf>`.

.. _annotation:

Annotation formats
------------------

**sideRETRO** needs to know the **exons** of the protein coding genes: Their
coordinates, their **strand** and the **gene** they belong to. The exon ids are
optional, they only name the exons at the VCF (``EEJLIST``, ``PG5P``, ``PG3P``).
By default, it expects a `GENCODE <https://www.gencodegenes.org/>`_ annotation:
Lines with the feature ``exon``, the attribute ``transcript_type=protein_coding``,
the gene name at ``gene_name`` and the exon id at ``exon_id``. Each of them can
be changed:

* ``-T``: The feature (third column).
* ``-A``: The attributes to filter in. Alternatives inside one option are
  separated by ``|`` and **any** of them must match; when the option is passed
  more than once, **all** of them must match. Values match **exactly** - and,
  in GFF3, a value packed with commas (``tag=basic,CCDS``) matches each one of
  its parts.
* ``-N``: The attribute with the gene name. If it is missing at an exon,
  ``gene_id`` is used.
* ``-E``: The attribute with the exon id. If it is missing, the exon is named
  after its coordinates, ``chr:start_end``.

For example, an `Ensembl <https://www.ensembl.org/>`_ GTF names the biotype
``transcript_biotype``::

  $ sider ps -a Homo_sapiens.GRCh38.gtf.gz \
      -A transcript_biotype=protein_coding \
      -i my_files_list.txt

while protein coding *and* lncRNA transcripts tagged as ``basic`` would be::

  $ sider ps -a gencode.gtf \
      -A 'transcript_type=protein_coding|lncRNA' \
      -A tag=basic \
      -i my_files_list.txt

Any other annotation can be given as a **BED6** file of exons, with the gene name
at the 4th column and the strand at the 6th one - the file name must end in
``.bed`` or ``.bed.gz``. The exons are named after their coordinates::

  $ head -2 exons.bed
  chr1	65418	65433	OR4F5	0	+
  chr1	65519	65573	OR4F5	0	+

  $ sider ps -a exons.bed -i my_files_list.txt

If no exon passes the filter, ``process-sample`` stops at once and says so.

.. note::
   GTF/GFF3 is read as a flat file: Each exon line must carry its own gene
   name. Annotations that link the exons to their genes only through the
   ``Parent`` attribute (as Ensembl GFF3) must be flattened before - for
   example, converted to GTF with `gffread
   <https://github.com/gpertea/gffread>`_ - or given as BED6.

Command ``merge-call``
======================

``merge-call`` (``mc``) is the second step: It merges the databases of
``process-sample``, clusters the abnormal alignments into retrocopy insertions
and annotates them - parental gene, breakpoints, target site, poly-A tail,
exon-exon junctions, retrocopy ends and predicted sequence - and genotypes each
one at every sample, reading the alignment files again at the insertion site.

Its arguments are the SQLite3 databases from ``process-sample``. Besides the
:ref:`common options <common_options>`:

.. list-table:: Mandatory options
   :header-rows: 1
   :widths: 25 60 15

   * - Option
     - Description
     - Default
   * - ``-r``, ``--reference-file``
     - FASTA file of the reference genome. The clipped reads at the junctions
       are searched into the parental exons, for the retrocopy ends and its
       predicted sequence. Only the exons are read, so it costs little
     -
   * - ``-i``, ``--input-file``
     - File with a newline separated list of databases. Not needed if the
       databases are passed as arguments; if both are given, their union is
       used
     -

.. list-table:: Output options
   :header-rows: 1
   :widths: 25 60 15

   * - Option
     - Description
     - Default
   * - ``-I``, ``--in-place``
     - Merge all the databases into the first one, instead of creating a new
       database
     -
   * - ``-c``, ``--cache-size``
     - SQLite3 cache size in KiB
     - 200000

.. list-table:: Clustering options
   :header-rows: 1
   :widths: 25 60 15

   * - Option
     - Description
     - Default
   * - ``-e``, ``--epsilon``
     - DBSCAN: Maximum distance between two alignments inside a cluster
     - 300
   * - ``-m``, ``--min-pts``
     - DBSCAN: Minimum number of alignments to form a dense region
     - 10

.. list-table:: Blacklist options
   :header-rows: 1
   :widths: 25 60 15

   * - Option
     - Description
     - Default
   * - ``-b``, ``--blacklist-chr``
     - Do not cluster from or to this chromosome. It may be passed multiple
       times
     - ``chrM``
   * - ``-B``, ``--blacklist-region``
     - Regions where the clusters are discarded, in GTF/GFF3 or BED format. A
       GTF/GFF3 file is filtered by ``-T`` and ``-A``
     -
   * - ``-P``, ``--blacklist-padding``
     - Extend the blacklisted regions by N bases at both sides
     - 0
   * - ``-T``, ``--gff-feature``
     - The feature (3rd column) of the blacklisted regions
     - ``gene``
   * - ``-A``, ``--gff-attribute``
     - An attribute (9th column) to filter in, as ``key=value``, with the
       same rules of ``process-sample``
     - ``gene_type=processed_pseudogene|tag=retrogene``
   * - ``-N``, ``--gff-gene-name``
     - The attribute naming the blacklisted regions
     - ``gene_name``

.. list-table:: Annotation options
   :header-rows: 1
   :widths: 25 60 15

   * - Option
     - Description
     - Default
   * - ``-x``, ``--parental-distance``
     - Minimum distance between a cluster and its parental gene. Closer
       clusters are discarded
     - 1000000
   * - ``-g``, ``--genotype-support``
     - Minimum number of reads from a single sample to keep a cluster
     - 1
   * - ``-n``, ``--near-gene-distance``
     - Maximum distance, in number of genes, between the parental genes of
       overlapping clusters to merge them as near parental genes
     - 3
   * - ``-D``, ``--max-tsd-del``
     - Maximum length of a deletion at the target site to merge the
       retrocopies at its two flanks, each with a single breakpoint. 0 merges
       none
     - 50000

.. list-table:: Genotyping options
   :header-rows: 1
   :widths: 25 60 15

   * - Option
     - Description
     - Default
   * - ``-t``, ``--threads``
     - Number of threads
     - 1
   * - ``-Q``, ``--phred-quality``
     - Minimum mapping quality (MAPQ) of the reads supporting the reference
       allele
     - 20

The databases can be given directly or at a list::

  $ sider merge-call -r my_reference_genome.fa -I database1.db database2.db
  $ sider merge-call -r my_reference_genome.fa -I -i my_databases_list.txt

.. note::
   Without ``-I``, ``merge-call`` copies all the databases into a new one,
   *out.db* (the prefix changes with ``-p``). With ``-I``, it merges them into
   the **first** database, which is overwritten - faster and with no copy, but
   keep a backup if you need the original. Do not give the new database the
   name of an input one.

A fuller example, with the short name of the command::

  $ sider mc \
      -i my_databases_list.txt \
      -r my_reference_genome.fa \
      -o output_dir \
      -p my_database \
      -l my_log_file.log \
      -c 2000000 \
      -B my_black_list.bed \
      -g 5 \
      -t 3

Here a cluster needs at least 5 reads from a single sample (``-g``), and the
regions at *my_black_list.bed* are not clustered (``-B``) - for example, the
centromeres, or regions known for false insertions. The genotypes are computed
with 3 threads.

To see ``merge-call`` chained in a real workflow, see :ref:`A Practical
Workflow <pract_wf>`.

Command ``make-vcf``
====================

``make-vcf`` (``vcf``) is the last step: It writes the retrocopies annotated by
``merge-call`` to a **VCF**, with their genotypes and filters, and to the
companion files - a **TSV** with the reads supporting each retrocopy, a
**FASTA** with the predicted sequences and another **FASTA** with the sequences
for primer design. All of them are described at :ref:`Output <chap_output>`.

Its argument is the database processed by ``merge-call``. It reads only the
bases around the insertion sites from the reference genome. Besides the
:ref:`common options <common_options>`:

.. list-table:: Mandatory options
   :header-rows: 1
   :widths: 25 60 15

   * - Option
     - Description
     - Default
   * - ``-r``, ``--reference-file``
     - FASTA file of the reference genome, for the REF column, the contigs,
       the target site and the primer file
     -

.. list-table:: Annotation options
   :header-rows: 1
   :widths: 25 60 15

   * - Option
     - Description
     - Default
   * - ``-n``, ``--near-gene-dist``
     - Maximum distance between an intergenic insertion and a gene to report
       it as near (``NEAR``)
     - 10000
   * - ``-e``, ``--orientation-error``
     - Maximum p-value of the Spearman's rho to report the orientation
       (``ORHO``)
     - 0.05
   * - ``-i``, ``--min-ir``
     - Minimum number of reads of a sample, with the mate inside an intron of
       the parental gene, to report the intron as retained (``IR``,
       ``IRLIST``)
     - 2
   * - ``-t``, ``--min-tsd``
     - Minimum length of a target site duplication (``TSITE=DUP``). A shorter
       overlap of the breakpoints is a blunt insertion (``BLUNT``): The
       microhomology at the junctions makes it ambiguous
     - 5
   * - ``-f``, ``--primer-flank``
     - Number of bases of the reference at each side of a junction at the
       primer file
     - 200

The ``FILTER`` column lists the checks that failed, or ``PASS``: A retrocopy
with its insertion point well defined, a hallmark of retrotransposition, a copy
longer than a read and at a primary chromosome. The filters are described at :ref:`FILTER <filter>`;
these options set their thresholds:

.. list-table:: Filter options
   :header-rows: 1
   :widths: 25 60 15

   * - Option
     - Description
     - Default
   * - ``-L``, ``--min-qual``
     - ``LowQual``: Minimum ``QUAL``
     - 20
   * - ``-Q``, ``--min-mapq``
     - ``LowMAPQ``: Minimum mean mapping quality of the supporting reads
     - 20
   * - ``-S``, ``--min-sr``
     - ``LowSR``: Minimum number of split reads at the breakpoints, whatever
       their side
     - 3
   * - ``-T``, ``--max-tsd``
     - ``TSDOutlier``: Maximum length of the target site duplication. It must
       not be shorter than ``--min-tsd``
     - 25
   * - ``-P``, ``--primary-chr``
     - ``NonPrimary``: Extended regular expression of the primary chromosomes.
       Any other one is an alt, decoy, unplaced or HLA contig
     - ``^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$``
   * - ``-m``, ``--min-pgmates``
     - ``NoParentalMate``: Minimum number of fragments with the mate aligned
       end to end inside the parental exons (``PGMATES``). Is there a copy
       longer than a read?
     - 1
   * - ``-x``, ``--no-filter``
     - Disable a filter: ``LowQual``, ``LowMAPQ``, ``AmbiguousParental``,
       ``LowSR``, ``NoBreakpoint``, ``OneBreakpoint``, ``TSDOutlier``,
       ``NoHallmark``, ``NonPrimary`` or ``NoParentalMate``. It may be passed
       multiple times
     -

So, for a database like *my_database.db*, type::

  $ sider make-vcf -r my_reference_genome.fa my_database.db

It writes *out.vcf*, *out.rnames.tsv*, *out.insertions.fa* and *out.primer.fa*.
With more options, and the short name of the command::

  $ sider vcf my_database.db \
      -r my_reference_genome.fa \
      -o output_dir \
      -p my_retrocopies \
      -l my_log_file.log \
      -n 50000 \
      -x NonPrimary

Here an intergenic insertion is reported near the genes up to 50 kb away
(``-n``), and the retrocopies at alt, decoy or HLA contigs may ``PASS``
(``-x NonPrimary``). The filters can be changed without running ``merge-call``
again.

.. _sim_test:

A quick test with simulated data
================================

**sideRETRO** ships a **simulator** of retrocopy insertions, which is used by
its own test suite. It builds a small random genome, inserts a retrocopy of a
three exon gene - with a TSD, a poly-A tail and, if you like, a 5' truncation -
and writes the alignments an aligner would produce, with no aligner needed. So
the whole pipeline can be exercised in seconds::

  $ tests/sim/retrocopy_sim.py --seed 11 --output-dir sim
  $ sider process-sample -s -a sim/genes.gtf -o sim sim/sim.sam
  $ sider merge-call -r sim/ref.fa -o sim -p merged sim/out.db
  $ sider make-vcf -r sim/ref.fa -o sim -p merged sim/merged.db
  $ tests/sim/retrocopy_sim.py --check sim/merged.vcf sim/truth.tsv
  OK: sim/merged.vcf matches sim/truth.tsv

The file *sim/truth.tsv* has the expected values and ``--check`` compares the
VCF with them. Run ``tests/sim/retrocopy_sim.py --help`` for the simulation
parameters: TSD length, poly-A length, number of discordant and split reads
etc. The resulting VCF line is dissected at :ref:`Output <chap_output>`.

.. _sqlite:

The SQLite3 databases
=====================

The results of each step are kept at a `SQLite3 <https://www.sqlite.org>`_
database. A few practical notes:

Locks
  SQLite3 locks the database while writing. If another process holds the lock -
  another **sideRETRO** on the same database, a backup, a ``sqlite3`` session -
  **sideRETRO** waits and tries again for up to **60 seconds** before giving up.
  The waiting time, in milliseconds, can be changed with the environment variable
  ``SIDER_BUSY_TIMEOUT``::

    $ SIDER_BUSY_TIMEOUT=600000 sider merge-call -I -r hg38.fa out.db

Storage
  Prefer a **local disk** for the databases. Network filesystems (NFS, Lustre,
  GPFS) have unreliable locks and may fail to write. Input/output errors are
  **not** retried: SQLite3 has already rolled back the work in progress, and
  going on could leave an incomplete database. **sideRETRO** stops and says why,
  so the step can be run again.

Errors
  When something goes wrong, the message names the operation, the database, the
  SQLite3 error and what to do about it::

    FATAL db.c:77: Preparing failed at database 'out.db': database is locked (SQLite error 5)
    SELECT major_version,minor_version FROM schema LIMIT 1
    The database is in use by another process. Wait for it to finish, or raise
    the waiting time with the environment variable SIDER_BUSY_TIMEOUT (in milliseconds)

.. _cram:

Dealing with CRAM format
========================

Working with CRAM files may be a little **tricky**, mainly if you have downloaded
the data from a public repository. Let's take a look at two possible cases:

* Local alignment
* External alignment

Local alignment
---------------

In order to generate an alignment file in the CRAM format, first we need to
index the reference genome:

.. code-block:: sh

   # Index for BWA: .fa.amb, .fa.ann, .fa.bwt, .fa.pac, .fa.sa files
   bwa index hg38.fa

   # Index reference genome for CRAM: .fa.fai file
   samtools faidx hg38.fa

Then, we can align with :code:`bwa`:

.. code-block:: sh

   # Align with BWA and generate a CRAM
   bwa mem hg38.fa file_R1.fastq file_R2.fastq | \
      samtools view -T hg38.fa -C -o file.cram -

The alignment :file:`file.cram` can be processed with :code:`sider`, as long as
we don't change the reference genome and its index (:file:`.fa.fai`) path. If so,
we need to set the environment variables :file:`REF_PATH` and :file:`REF_CACHE`,
see :ref:`External alignment <extern_cram>`.

.. _extern_cram:

External alignment
------------------

When we download public data already aligned in the CRAM format, we may be
concerned about the reference genome index. Probably,  we won't have the
required genome index to read the :file:`.cram`, and the :file:`htslib`
library - used by :code:`sider` and :code:`samtools` - is able to download
the index from the `CRAM Reference Registry <http://www.ebi.ac.uk/ena/cram>`_.

However, in order to :file:`htslib` be able to accomplish this task, we need
to compile the library with the required flags and also we need to have the
required dependencies (as `libcurl <https://curl.haxx.se/libcurl/>`_).
Therefore to be able to read these files, without depending on these details,
we need to generate a new local index and set the environment variables -
:file:`REF_PATH` and :file:`REF_CACHE` - to the correct path:

.. code-block:: sh

   # Create cache dir
   mkdir -p /my/cache

   # Construct the index
   perl seq_cache_populate.pl -root /my/cache hg38.fa

   # Now before running samtools or sider, we need to
   # set the environment variables REF_PATH and REF_CACHE
   export REF_PATH=/my/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram
   export REF_CACHE=/my/cache/%2s/%2s/%s

   # So ...
   sider ps -a annot.gff3.gz -o result file.cram

The script :file:`seq_cache_populate.pl` can be found in the :file:`samtools`,
or at `seq_cache_populate.pl
<https://github.com/deweylab/RSEM/blob/master/samtools-1.3/misc/seq_cache_populate.pl>`_.

For more information, see `Samtools Workflow <https://www.htslib.org/workflow/>`_.

.. _pract_wf:

A Practical Workflow
====================

Now, let's do an interesting exercise, with real experimental data from the
`1000 Genomes Project <https://www.internationalgenome.org/>`_. (Warning: This example requires 16GB of RAM)

In order to run **sideRETRO** searching for retrocopies, we will download 2
whole-genome sequenced CRAM files, both aligned on the **gencode**'s
`hg38 <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz>`_
genome:
`NA12878 <ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram>`_
and
`NA12778 <ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239484/NA12778.final.cram>`_.

At the beginning of a run, the files listed below must be at the same directory where
the user is running **sideRETRO** or their correct paths must be supplied at the
correspondent option. Files are:

1. A GTF gene annotation file from gencode project
   (here :file:`gencode.v32.annotation.gtf`).

2. A FASTA file with the gencode's Human reference genome, version 38
   (here :file:`GRCh38_full_analysis_set_plus_decoy_hla.fa`).

3. A custom perl script, :code:`seq_cache_populate.pl`, to construct a new local index . 
   The :code:`seq_cache_populate.pl` script can be found in 
   `seq_cache_populate.pl <https://github.com/deweylab/RSEM/blob/master/samtools-1.3/misc/seq_cache_populate.pl>`_.

4. A custom perl script, :code:`analyser.pl`, to do the final analysis over the VCF file
   and produce the TSV file in a tabular format. The :code:`analyser.pl` script can be
   downloaded :download:`here <data/analyser.pl>`.

Also, we will set the environment variables :file:`REF_PATH` and :file:`REF_CACHE`, as
a requirement to work with CRAM files - more information at
:ref:`Dealing with CRAM format<extern_cram>`.

See the complete command sequence below for the whole analysis.

Tip: Copy and paste line by line in your terminal.

Tip 2: If you are running line by line in your terminal don't paste the "$" character. It is already in your terminal.

.. code-block:: sh

  # Do things inside a clean directory.
  # Average time: irrelevant
  $ mkdir -p sider_test
  $ cd sider_test

  # Download annotation from gencode
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz

  # Download the reference genome from 1000 genomes
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

  # Make the CRAM index
  # Create cache dir
  mkdir -p cache

  # create index
  perl seq_cache_populate.pl -root cache GRCh38_full_analysis_set_plus_decoy_hla.fa

  # Set environment variables
  export REF_PATH=$PWD/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram
  export REF_CACHE=$PWD/cache/%2s/%2s/%s

  # Create a download list (WGS.list) containing all files of interest.
  # Average time: irrelevant
  $ echo "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram" > WGS_download.list
  $ echo "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239484/NA12778.final.cram" >> WGS_download.list

  # Download all files: NA12878 and NA12778.
  # Average time: network dependent
  $ wget -c -i WGS_download.list

  # Create the list of BAM files.
  # Average time: irrelevant
  $ ls *.cram > WGS_genomes.list

  # First sideRETRO step: process-sample
  # Input file: WGS_genomes.list
  # Output file: 1000_genomes.db
  # Average time: 62m34.541
  $ sider process-sample \
      -i WGS_genomes.list \
      -a gencode.v32.annotation.gtf.gz \
      -p 1000_genomes \
      -c 2000000 \
      -F 0.9 \
      -t 2

  # Second sideRETRO step: merge-call
  # Input file: 1000_genomes.db
  # Output file: 1000_genomes.db (same file)
  # Average time: 62m34.541
  $ sider merge-call 1000_genomes.db \
      -c 2000000 \
      -x 1000000 \
      -g 5 \
      -r GRCh38_full_analysis_set_plus_decoy_hla.fa \
      -I \
      -t 2

  # Third sideRETRO step: make-vcf
  # Input file: 1000_genomes.db
  # Output files: 1000_genomes.vcf, 1000_genomes.rnames.tsv,
  #               1000_genomes.insertions.fa and 1000_genomes.primer.fa
  # Average time: 62m34.541
  $ sider make-vcf 1000_genomes.db \
      -p 1000_genomes \
      -r GRCh38_full_analysis_set_plus_decoy_hla.fa

  # Some analysis over the final VCF file.
  # Input file: 1000_genomes.vcf
  # Output file: 1000_genomes.tsv
  # Average time: 62m34.541
  $ perl analyser.pl 1000_genomes.vcf > 1000_genomes.tsv

This was a simple but complete pipeline to obtain a final TSV file with all
the relevant results in a tabular format ready to import in a R or Python script
and plot some graphics.

.. _run_dck:

Running with Docker
===================

Notwithstanding **sideRETRO**'s native run, user can happily run it from a
**Docker** image just prepending **Docker**'s directives to any example shown.
That is, with *Docker* installed and the image *galantelab/sider:latest* pulled
from `DockerHub <https://hub.docker.com/r/galantelab/sider>`_, replace the
ordinary ``sider`` command by
``docker run --rm -ti -v $(pwd):/home/sider -w /home/sider galantelab/sider``
- the image runs ``sider`` itself - like::

  $ docker run --rm -ti -v $(pwd):/home/sider -w /home/sider galantelab/sider \
    ps \
        -i my_files_list.txt \
        -a annotation_file.gtf \
        -o output_dir \
        -p my_reads_database \
        -l my_log_file.log \
        -c 2000000 \
        -F 0.9 \
        -t 3
