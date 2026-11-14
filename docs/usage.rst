.. _chap_usage:

***************
Using sideRETRO
***************

General Syntax
==============

**sideRETRO** has a very straightforward syntax. Basically, there are three main
commands, each one with a plethora of available options:

* ``process-sample``
* ``merge-call``
* ``make-vcf``

So, in order to test the installation process and run a first example, user can
call it without any argument from the command line, like this::

  $ sider
  Usage: sider [-hv]
         sider <command> [options]

  A pipeline for detecting
  Somatic Insertion of DE novo RETROcopies

  Options:
     -h, --help            Show help options
     -v, --version         Show current version
     -c, --cite            Show citation in BibTeX

  Commands:
     ps,  process-sample   Extract alignments related
                           an event of retrocopy
     mc,  merge-call       Discover and annotate
                           retrocopies
     vcf, make-vcf         Generate VCF file with all
                           annotate retrocopies

In the above situation, if **sideRETRO** was correctly installed, it will give
that default *usage* help.

Another classical example is to print **sideRETRO**'s installed version using
the ``-v`` option::

  $ sider --version
  sideRETRO 1.7.0

And, if the user need further help, he can find it both at the **sideRETRO**'s
`readthedocs page <https://sideretro.readthedocs.io>`_ or in the already
installed software documentation, from command line::

  $ sider --help

Please, see :ref:`A Practical Workflow <pract_wf>` and :ref:`Running with Docker
<run_dck>` sections for more examples and tips for using with **Docker**.

Now, to get more familiar with **sideRETRO** main commands and results, let's
try some basic examples for each command.

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

The first one is ``process-sample`` or ``ps`` for short, and was intended to act
as the *"evidence's grounding faith"* for **sideRETRO**. Here, we're saying
"first" because of an order in which the user must run the commands. The file
resultant from this command will become the input to the next one,
``merge-call``.

As explained in the `Introduction <intro.rst>`_ section, the command
``process-sample`` creates a database of abnormal reads from a SAM/BAM file set.
To do this, there are some mandatory options the user must supply to do a
correct search. Calling the command ``process-sample`` without any argument
will give a specific help where user can know all the mandatory options for
this command::

  $ sider process-sample

  Arguments:
     One or more alignment files in SAM/BAM/CRAM format

  Mandatory Options:
     -a, --annotation-file   Exons of the protein coding genes, in GTF/GFF3
                             or BED6 format. See 'Annotation Options'
     -i, --input-file        File with a newline separated list of alignment
                             files in SAM/BAM/CRAM format. Not mandatory if
                             the files are passed as arguments. If both are
                             given, their union is used

  Input/Output Options:
     -h, --help              Show help options
     -q, --quiet             Decrease verbosity to error messages only,
                             or suppress terminal output at all if
                             'log-file' is passed
         --silent            Same as '--quiet'
     -d, --debug             Increase verbosity to debug level
     -l, --log-file          Print log messages to a file
     -o, --output-dir        Output directory. It is created if it does not
                             exist [default:"."]
     -p, --prefix            Prefix of the output files [default:"out"]

  Annotation Options:
     A BED6 file (name ending in '.bed' or '.bed.gz') has one exon per
     line, the gene name at the 4th column and the strand at the 6th one.
     A GTF/GFF3 file is filtered by the options below
     -T, --gff-feature       The 'feature' (3rd column) of the exons
                             [default:"exon"]
     -A, --gff-attribute     An 'attribute' (9th column) to filter in, as
                             key=value. Alternatives separated by '|' must
                             match any (OR): 'k=v1|v2' or 'k1=v1|k2=v2'. If
                             passed multiple times, all of them must match
                             (AND). Values match exactly
                             [default:"transcript_type=protein_coding"]
     -N, --gff-gene-name     The attribute with the gene name. If missing,
                             'gene_id' is used [default:"gene_name"]
     -E, --gff-exon-id       The attribute with the exon id. If missing, the
                             exon is named 'chr:start_end' [default:"exon_id"]

  SQLite3 Options:
     -c, --cache-size        SQLite3 cache size in KiB [default:"200000"]

  Read Quality Options:
     -Q, --phred-quality     Minimum mapping quality (MAPQ) of the reads
                             [default:"8"]
     -M, --max-base-freq     Maximum frequency of a single base in a read.
                             Reads above it are discarded as low complexity
                             sequences, such as poly-A [default:"0.75"]
     -D, --deduplicate       Remove duplicated reads: Fragments whose reads
                             and mates share the same 5' positions

  Processing Options:
     -s, --sorted            Assume the reads are grouped by read name, even
                             without the header tag 'SO:queryname'. Grouped
                             files are read once, the others twice
     -t, --threads           Number of threads. Each file uses one thread and
                             the spare ones decompress BAM/CRAM [default:"1"]
     -m, --max-distance      Maximum distance between the reads of a pair.
                             Farther pairs are abnormal [default:"10000"]
     -f, --exon-frac         Minimum overlap between an alignment and an
                             exon, as a fraction of the exon
                             [default:"1e-09"; 1 base]
     -F, --alignment-frac    Minimum overlap between an alignment and an
                             exon, as a fraction of the alignment
                             [default:"1e-09"; 1 base]
     -e, --either            Either '-f' OR '-F' must be satisfied. Without
                             '-e', both must be
     -r, --reciprocal        The overlap must be reciprocal: '-F' is set to
                             the value of '-f'


So, supposing that the user has three files: *f1.bam*, *f2.bam*, *f3.sam*, he
can type::

  $ sider process-sample f2.bam f2.bam f3.sam \
      -a annotation_file.gtf

Note the mandatory ``-a`` option specifying the annotation file. And, in this
unique exception, we suppressed the ``-i`` mandatory option cause all the files
were explicitly called.

Let's see another example that shows the convenient use of the ``-i`` option to
call a list of input files (e.g. *my_files_list.txt*) instead of them directly::

  $ sider process-sample \
      -i my_files_list.txt \
      -a annotation_file.gtf

Both commands above will produce only one output database file *out.db*
containing all relevant reads for non-fixed retrocopies search, whose prefix
*out* can be easily changed with the ``-p`` option. The abnormal reads from
all input files will be merged in just one table, **with their sequences**,
which are needed to annotate the junctions later. To produce one database
for each input file separately, user must run one distinct instance of
**sideRETRO** per file.

.. note::
   The alignment files do not need to be sorted. When they are **not**
   grouped by read name (``SO:queryname``), each file is read **twice**
   in order to catch all the alignments of each abnormal fragment. If
   there are more threads (``-t``) than files, the spare ones are used to
   decompress the BAM/CRAM.

Some options' values can affect drastically the output. Let's play a little bit
with some of them while using the short version of the command ``ps``::

  $ sider ps \
      -i my_files_list.txt \
      -a annotation_file.gtf \
      -o output_dir \
      -p my_reads_database \
      -l my_log_file.log \
      -c 2000000 \
      -Q 20 \
      -F 0.9 \
      -t 3

Wow! The number of options can be overwhelming.

Here used ``-o`` option to specify the directory *output_dir* to write our
database as *my_reads_database.db* (``-p`` option). Also, we chose to save the
log messages in *my_log_file.log* file (``-l`` option), a cache size of 2Gb
(``-c`` option), a minimum phred score cutoff of 20 for alignments (``-Q``
option), a minimum overlap ratio of 0.9 for read alignments over exonic regions
(``-F`` option) and 3 threads to process those files in parallel (``-t`` option).

To see another example of the ``process-sample`` command chained in a real
workflow, please refer to the :ref:`A Practical Workflow <pract_wf>` section.

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

The second step in the **sideRETRO**'s *"journey for the truth of retrocopies"*
is the command ``merge-call`` or ``mc`` for short. The aim of this command is to
take the database created by ``process-sample`` step as input and populate more
tables in it, with information risen from a clustering process over the abnormal
reads regions.

Like ``process-sample``, ``merge-call`` has some mandatory options, which can be
known by calling it without any argument::

  $ sider merge-call

  Arguments:
     One or more SQLite3 databases generated by 'process-sample'

  Mandatory Options:
     -i, --input-file           File with a newline separated list of SQLite3
                                databases. Not mandatory if the databases are
                                passed as arguments. If both are given, their
                                union is used

  Input/Output Options:
     -h, --help                 Show help options
     -q, --quiet                Decrease verbosity to error messages only,
                                or suppress terminal output at all if
                                'log-file' is passed
         --silent               Same as '--quiet'
     -d, --debug                Increase verbosity to debug level
     -l, --log-file             Print log messages to a file
     -o, --output-dir           Output directory. It is created if it does not
                                exist [default:"."]
     -p, --prefix               Prefix of the output files [default:"out"]
     -I, --in-place             Merge all databases into the first one, instead
                                of creating a new database

  SQLite3 Options:
     -c, --cache-size           SQLite3 cache size in KiB [default:"200000"]

  Clustering Options:
     -e, --epsilon              DBSCAN: Maximum distance between two alignments
                                inside a cluster [default:"300"]
     -m, --min-pts              DBSCAN: Minimum number of alignments required to
                                form a dense region [default:"10"]

  Blacklist Options:
     -b, --blacklist-chr        Do not cluster from or to this chromosome. It may
                                be passed multiple times [default:"chrM"]
     -B, --blacklist-region     Regions where clusters are discarded, in GTF/GFF3
                                or BED format. A GTF/GFF3 file is filtered by
                                '-T' and '-A'
     -P, --blacklist-padding    Extend the blacklisted regions by N bases at both
                                sides [default:"0"]
     -T, --gff-feature          The 'feature' (3rd column) of the blacklisted
                                regions [default:"gene"]
     -A, --gff-attribute        An 'attribute' (9th column) to filter in, as
                                key=value. Alternatives separated by '|' must
                                match any (OR): 'k=v1|v2' or 'k1=v1|k2=v2'. If
                                passed multiple times, all of them must match
                                (AND). Values match exactly
                                [default:"gene_type=processed_pseudogene|tag=retrogene"]
     -N, --gff-gene-name        The attribute naming the blacklisted regions
                                [default:"gene_name"]

  Annotation Options:
     -x, --parental-distance    Minimum distance between a cluster and its
                                parental gene. Closer clusters are discarded
                                [default:"1000000"]
     -g, --genotype-support     Minimum number of reads from a single sample
                                (SAM/BAM/CRAM) to keep a cluster [default:"1"]
     -n, --near-gene-rank       Maximum distance, in number of genes, between
                                the parental genes of overlapping clusters to
                                merge them as near parental genes
                                [default:"3"]
     -r, --reference-file       FASTA file of the reference genome. With it,
                                the retrocopy ends are searched in the parental
                                exons and its sequence is predicted

  Genotyping Options:
     -t, --threads              Number of threads [default:"1"]
     -Q, --phred-quality        Minimum mapping quality (MAPQ) of the reads
                                supporting the reference allele [default:"8"]


And likewise, user can call a set of database files directly, or using a list of
files::

  $ sider merge-call database1.db database2.db -I

or ::

  $ sider merge-call -i my_databases_list.txt -I

.. note::
   Again, note the ``-I`` option that is not mandatory but would lead the creation
   of duplicated output databases if absent. This option do the clustering
   "in place" over the input files, overwriting them (so be careful). If user do
   not use the ``-p`` or ``-I`` options, the output files will be named *out.db*.

In a more sophisticated example, we will use the short version of the command
``mc``, with many other options::

  $ sider mc \
      -i my_databases_list.txt \
      -o output_dir \
      -p my_database \
      -l my_log_file.log \
      -I \
      -c 2000000 \
      -B my_black_list.bed \
      -x 1000000 \
      -g 5 \
      -Q 20 \
      -r my_reference_genome.fa \
      -t 3

Here, options ``-i``, ``-o``, ``-p``, ``-l``, ``-I``, ``-c``, ``-Q`` and ``-t``
keep the same meaning as they have in the ``process-sample`` command.
The others need some explanation. All we've done here was to ask for a minimum
number of 5 reads of contribution from each input SAM/BAM file to consider a
clustering region as a retrocopy candidate (with ``-g`` option); a minimum
distance of 1000000 bp from the parental gene to resolve some doubtful overlaps
(``-x`` option) and, importantly, a BED file with a list of regions to be
ignored at the clustering process called *my_black_list.bed* (``-B`` option).
This last option's file can describe entire chromosomes (like chrM) and many
chromosomal regions with poor insertion evidences taken from literature, like
centromeres. All specified regions won't be targets for clustering.

The option ``-r`` deserves a special note. ``merge-call`` finds the
**breakpoints**, the **TSD** length, the **poly-A** tail and the **exon-exon
junctions** from the reads alone. With the reference genome, it goes one step
further: It searches the clipped sequences at the junctions into the exons of
the parental gene, in order to find **where** the retrocopy begins and ends
(the 5' truncation) and to **predict** its sequence. See the
:ref:`methodology <chap_methodology>` for the details. The option is not
mandatory, but we recommend it: Only the exons of the parental genes are
read from the FASTA, so it costs little.

To see another example of the ``merge-call`` command chained in a real workflow,
please refer to the :ref:`A Practical Workflow <pract_wf>` section.

Command ``make-vcf``
====================

The third and last step to the **sideRETRO**'s *"crusade to retrocopies"* is the
``make-vcf`` command or ``vcf`` for short. This command takes the already
clustered tables in the database files populated at the ``merge-call`` step and
writes the results: A **VCF** file with all retrocopy insertions annotated, a
**TSV** file with the reads supporting each one, and a **FASTA** file with the
predicted sequences. All of them are described at :ref:`Output <chap_output>`.

This command has no mandatory options, but it is worth try to discover the
others::

  $ sider make-vcf

  Arguments:
     SQLite3 database processed by 'merge-call'

  Input/Output Options:
     -h, --help                 Show help options
     -q, --quiet                Decrease verbosity to error messages only,
                                or suppress terminal output at all if
                                'log-file' is passed
         --silent               Same as '--quiet'
     -d, --debug                Increase verbosity to debug level
     -l, --log-file             Print log messages to a file
     -o, --output-dir           Output directory. It is created if it does not
                                exist [default:"."]
     -p, --prefix               Prefix of the output files [default:"out"]

  Annotation Options:
     -n, --near-gene-dist       Maximum distance between an intergenic insertion
                                and a gene to report it as near (NEAR)
                                [default:"10000"]
     -e, --orientation-error    Maximum p-value of the Spearman's rho to report
                                the orientation (ORHO) [default:"0.05"]
     -r, --reference-file       FASTA file of the reference genome, for the REF
                                column, the contigs and the TSD sequence

  Filter Options:
     The FILTER column lists the checks of a canonical retrocopy that
     failed, or PASS. A check that cannot be evaluated does not fail
     -L, --min-qual             LowQual: Minimum QUAL [default:"20"]
     -Q, --min-mapq             LowMAPQ: Minimum mean mapping quality of the
                                supporting reads [default:"20"]
     -x, --no-filter            Disable a filter. It may be passed multiple
                                times: LowQual, LowMAPQ, AmbiguousParental,
                                NoTSD, TSDOutlier, NoPolyA, NoJunction and
                                PolarityConflict


So, in order to produce the results from a database input file like
*my_database.db*, just type::

  $ sider make-vcf my_database.db

This will produce the *out.vcf*, *out.rnames.tsv* and *out.insertions.fa*
files.

Let's add more options to customize it to our needs (with the short version of
the command only for symmetry)::

  $ sider vcf my_database.db \
      -o output_dir \
      -p my_retrocopies \
      -l my_log_file.log \
      -r my_reference_genome.fa \
      -n 50000

Command ``make-vcf`` is very simple and doesn't allow the user to use threads.
The only new options are ``-r``, which must specify the reference genome in
FASTA format (like **gencode**'s *Hg38.fa*) and ``-n``, where user can establish
a distance threshold for genes surrounding insertion points for additional
information in the output VCF file. Without ``-r``, the VCF ``REF`` column is
``N`` and the TSD sequence is not reported. Only the bases at the insertion
sites are read from the FASTA, not the whole genome.

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

    $ SIDER_BUSY_TIMEOUT=600000 sider merge-call -I out.db

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
      -Q 20 \
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
  # Output files: 1000_genomes.vcf, 1000_genomes.rnames.tsv
  #               and 1000_genomes.insertions.fa
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
That is, supposing the user has *Docker* installed and has pulled the image
*galantelab/sider:latest* from `DockerHub
<https://hub.docker.com/r/galantelab/sider>`_, he can just prepend
``docker --rm -ti -v $(pwd):/home/sider -w /home/sider galantelab/sider``
to the ordinary ``sider`` command, like::

  $ docker --rm -ti -v $(pwd):/home/sider -w /home/sider galantelab/sider \
    sider ps \
        -i my_files_list.txt \
        -a annotation_file.gtf \
        -o output_dir \
        -p my_reads_database \
        -l my_log_file.log \
        -c 2000000 \
        -Q 20 \
        -F 0.9 \
        -t 3
