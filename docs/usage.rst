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
  sideRETRO 1.0.0

And, if the user need further help, he can find it both at the **sideRETRO**'s
`readthedocs page <https://sideretro.readthedocs.io>`_ or in the already
installed software documentation, from command line::

  $ sider --help

Please, see :ref:`A Practical Workflow <pract_wf>` and :ref:`Running with Docker
<run_dck>` sections for more examples and tips for using with **Docker**.

Now, to get more familiar with **sideRETRO** main commands and results, let's
try some basic examples for each command.


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
   One or more alignment file in SAM/BAM format

Mandatory Options:
  -a, --annotation-file   Gene annotation on the reference genome
                          in GTF/GFF3 format. sider will look for 'exon'
                          with the attribute 'transcript_type=protein_coding'.
                          The attributes 'gene_name', 'gene_id' and 'exon_id'
                          are also required
  -i, --input-file        File containing a newline separated list of
                          alignment files in SAM/BAM/CRAM format.
                          This option is not mandatory if one or more
                          SAM/BAM/CRAM files are passed as argument.
                          If 'input-file' and arguments are set
                          concomitantly, then the union of all alignment
                          files is used

Input/Output Options:
  -h, --help              Show help options
  -q, --quiet             Decrease verbosity to error messages only
                          or suppress terminal outputs at all if
                          'log-file' is passed
  --silent                Same as '--quiet'
  -d, --debug             Increase verbosity to debug level
  -l, --log-file          Print log messages to a file
  -o, --output-dir        Output directory. Create the directory if it does
                          not exist [default:"."]
  -p, --prefix            Prefix output files [default:"out"]

SQLite3 Options:
  -c, --cache-size        Set SQLite3 cache size in KiB [default:"200000"]

Read Quality Options:
  -Q, --phred-quality     Minimum mapping quality of the reads required
                          [default:"8"]
  -M, --max-base-freq     Maximum base frequency fraction allowed
                          [default:"0.75"]
  -D, --deduplicate       Remove duplicated reads. Reads are considered
                          duplicates when they share the 5 prime positions
                          of both reads and read-pairs

Processing Options:
  -s, --sorted            Assume all reads are grouped by queryname, even if
                          there is no SAM/BAM/CRAM header tag 'SO:queryname'
  -t, --threads           Number of threads [default:"1"]
  -m, --max-distance      Maximum distance allowed between paired-end reads
                          [default:"10000"]
  -f, --exon-frac         Minimum overlap required as a fraction of exon
                          [default:"1e-09"; 1 base]
  -F, --alignment-frac    Minimum overlap required as a fraction of
                          alignment [default:"1e-09"; 1 base]
  -e, --either            The minimum fraction must be satisfied for at least
                          exon OR alignment. Without '-e', both fractions would
                          have to be satisfied
  -r, --reciprocal        The fraction overlap must be reciprocal for exon and
                          alignment. If '-f' is 0.5, then '-F' will be set to
                          0.5 as well

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
all input files will be merged in just one table. To produce one database for
each input file separately, user must run one distinct instance of
**sideRETRO** per file.

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
   One or more SQLite3 databases generated in the `process-sample
   <#command-process-sample>`_ step

Mandatory Options:
   -i, --input-file           File containing a newline separated list of
                              SQLite3 databases to be processed. This
                              option is not mandatory if one or more
                              SQLite3 databases are passed as argument.
                              If 'input-file' and arguments are set
                              concomitantly, then the union of all files
                              is used

Input/Output Options:
   -h, --help                 Show help options
   -q, --quiet                Decrease verbosity to error messages only
                              or suppress terminal outputs at all if
                              'log-file' is passed
   --silent                   Same as '--quiet'
   -d, --debug                Increase verbosity to debug level
   -l, --log-file             Print log messages to a file
   -o, --output-dir           Output directory. Create the directory if it does
                              not exist [default:"."]
   -p, --prefix               Prefix output files [default:"out"]
   -I, --in-place             Merge all databases with the first one of the list,
                              instead of creating a new file

SQLite3 Options:
   -c, --cache-size           Set SQLite3 cache size in KiB [default:"200000"]

Clustering Options:
   -e, --epsilon              DBSCAN: Maximum distance between two alignments
                              inside a cluster [default:"300"]
   -m, --min-pts              DBSCAN: Minimum number of points required to form a
                              dense region [default:"10"]

Filter & Annotation Options:
   -b, --blacklist-chr        Avoid clustering from and to this chromosome. This
                              option can be passed multiple times [default:"chrM"]
   -B, --blacklist-region     GTF/GFF3/BED blacklisted regions. If the file is in
                              GTF/GFF3 format, the user may indicate the 'feature'
                              (third column), the 'attribute' (ninth column) and
                              its values
   -P, --blacklist-padding    Increase the blacklisted regions ranges (left and right)
                              by N bases [default:"0"]
   -T, --gff-feature          The value of 'feature' (third column) for GTF/GFF3
                              file [default:"gene"]
   -H, --gff-hard-attribute   The 'attribute' (ninth column) for GTF/GFF3
                              file. It may be passed in the format key=value
                              (e.g. gene_type=pseudogene). Each value will match
                              as regex, so 'pseudogene' can capture IG_C_pseudogene,
                              IG_V_pseudogene etc. This option can be passed multiple
                              times and must be true in all of them
   -S, --gff-soft-attribute   Works as 'gff-hard-attribute'. The difference is
                              if this option is passed multiple times, it needs
                              to be true only once
                              [default:"gene_type=processed_pseudogene tag=retrogene"]
   -x, --parental-distance    Minimum distance allowed between a cluster and
                              its putative parental gene [default:"1000000"]
   -g, --genotype-support     Minimum number of reads coming from a given source
                              (SAM/BAM/CRAM) within a cluster [default:"3"]
   -n, --near-gene-rank       Minimum ranked distance between genes in order to
                              consider them close [default:"3"]

Genotyping Options:
   -t, --threads              Number of threads [default:"1"]
   -Q, --phred-quality        Minimum mapping quality used to define reference
                              allele reads [default:"8"]


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
      -C 15 \
      -t 3

Here, options ``-i``, ``-o``, ``-p``, ``-l``, ``-I``, ``-c``, ``-Q`` and ``-t``
keeps the same meaning as they have in the ``process-sample`` command.
The others need some explanation. All we've done here was to ask for a minimum
number of 5 reads of contribution from each input SAM/BAM file to consider a
clustering region as a retrocopy candidate (with ``-g`` option); a minimum
distance of 1000000 bp from the parental gene to resolve some doubtful overlaps
(``-x`` option), a minimum number of 15 crossing reads over the putative
insertion point to consider heterozygosis evidence (``-C``) and, importantly,
a BED file with a list of regions to be ignored at the clustering process called
*my_black_list.txt* (``-B`` option). This last option's file can describe
entire chromosomes (like chrM) and many chromosomal regions with poor insertion
evidences taken literature, like centromers. All specified regions won't be
targets for clustering.

To see another example of the ``merge-call`` command chained in a real workflow,
please refer to the :ref:`A Practical Workflow <pract_wf>` section.

Command ``make-vcf``
====================

The third and last step to the **sideRETRO**'s *"crusade to retrocopies"* is the
``make-vcf`` command or ``vcf`` for short. This command takes the already
clustered tables in the database files populated at the ``merge-call`` step and
creates one VCF file with all statistically significant retrocopy insertions
annotated in a convenient format.

This command has no mandatory options, but it is worth try to discover the
others::

  $ sider make-vcf

Arguments:
   SQLite3 database generated at `process-sample <#command-process-sample>`_
   and `merge-call <#command-merge-call>`_ steps

Input/Output Options:
   -h, --help                 Show help options
   -q, --quiet                Decrease verbosity to error messages only
                              or suppress terminal outputs at all if
                              'log-file' is passed
   --silent                   Same as '--quiet'
   -d, --debug                Increase verbosity to debug level
   -l, --log-file             Print log messages to a file
   -o, --output-dir           Output directory. Create the directory if it does
                              not exist [default:"."]
   -p, --prefix               Prefix output files [default:"out"]

Filter & Annotation Options:
   -n, --near-gene-dist       Minimum distance between genes in order to
                              consider them close [default:"10000"]
   -e, --orientation-error    Maximum error allowed for orientation rho
                              [default:"0.05"]
   -r, --reference-file       FASTA file for the reference genome

So, in order to produce a VCF file from a database input file like
*my_database.db*, just type::

  $ sider make-vcf my_database.db

This will produce a *out.vcf* output file.

Let's add more options to customize it to our needs (with the short version of
the command only for symmetry)::

  $ sider vcf my_database.db \
      -o output_dir \
      -p my_retrocopies \
      -l my_log_file.log \
      -r my_reference_genome.fa \
      -n 50000

Command ``make-vcf`` is very simple and don't allow the user to use threads.
The only new options are ``-r``, which must specify the reference genome in
FASTA format (like **gencode**'s *Hg38.fa*) and ``-n``, where user can establish
a distance threshold for genes surrounding insertion points for additional
information in the output VCF file.

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

   # Inde for BWA: .fa.amb, .fa.ann, .fa.bwt, .fa.pac, .fa.sa files
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
reqeuired dependencies (as `libcurl <https://curl.haxx.se/libcurl/>`_).
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

For more information, see `Samtools Worflow <https://www.htslib.org/workflow/>`_.

.. _pract_wf:

A Practical Workflow
====================

Now, let's do an interesting exercise, with real experimental data from the
`1000 Genomes Project <https://www.internationalgenome.org/>`_. (Warning: This example requires 16GB of RAM)

In order to run **siderRETRO** searching for retrocopies, we will download 2
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
      -I \
      -t 2

  # Second sideRETRO step: merge-call
  # Input file: 1000_genomes.db
  # Output file: 1000_genomes.vcf
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
