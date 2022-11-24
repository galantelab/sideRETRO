<p align="center"><a href="https://sideretro.readthedocs.io/en/latest/?badge=latest"><img src="docs/images/logo_sideRETRO.png" alt="sideRETRO" width="200"></a></p>
<h2 align="center">A pipeline for detecting Somatic Insertion of DE novo RETROcopies</h2>

<p align="center">
  <a href="https://github.com/galantelab/sideRETRO/actions/workflows/ci.yml"><img alt="" src="https://github.com/galantelab/sideRETRO/actions/workflows/ci.yml/badge.svg?branch=master" align="center"></a>
  <a href="https://sideretro.readthedocs.io/en/latest/?badge=latest"><img alt="" src="https://readthedocs.org/projects/sideretro/badge/?version=latest" align="center"></a>
  <a href="https://coveralls.io/github/galantelab/sideRETRO?branch=master"><img alt="" src="https://coveralls.io/repos/github/galantelab/sideRETRO/badge.svg?branch=master" align="center"></a>
</p>

**sideRETRO** is a bioinformatic tool devoted for the detection of somatic **retrocopy** insertion, also known as
**retroCNV**, in whole genome and whole exome sequencing data (WGS, WES). The program has been written from scratch
in C, and uses [HTSlib](http://www.htslib.org/) and [SQLite3](https://www.sqlite.org) libraries, in order to manage
**SAM/BAM/CRAM** reading and data analysis.

For full documentation, please visit <https://sideretro.readthedocs.io>.

### Features

When detecting retrocopies, **sideRETRO** can annotate several other features related to each event:

* **Parental gene**

   The gene which underwent retrotransposition process.

* **Genomic position**

   The genome coordinate where occurred the retrocopy integration event (chromosome:start-end).
   It includes the insertion point (the expected exact point of each retrocopy insertion).

* **Strandness**

   Detects the orientation of the insertion (+/-). It takes into account the orientation of insertion,
   whether in the leading (+) or lagging (-) DNA strand.

* **Genomic context**

   The retrocopy integration site context: If the retrotransposition event occurred at an intergenic or
   intragenic region - the latter can be splitted into exonic and intronic according to the host gene.

* **Genotype**

   When multiple individuals (genomes) are analyzed, **sideRETRO** discriminates events found in each one.
   That way, it is possible to distinguish whether an event is exclusive or shared among the cohort analyzed.

* **Haplotype**

   Our tool provides information about the ploidy of the event, i.e., whether it occurs in one or both homologous
   chromosomes (homozygous or heterozygous).

## Getting Started

### Installation

The project depends on [Meson build system](https://mesonbuild.com) and [Ninja](https://github.com/ninja-build/ninja)
to manage configuration and compilation process. They can be obtained using package manager or from source. For example,
using [Ubuntu](https://ubuntu.com) distribution:

```
$ sudo apt-get install python3 \
                       python3-pip \
                       python3-setuptools \
                       python3-wheel \
                       ninja-build
```

and then:

`$ pip3 install --user meson`

(or: `$ sudo apt install meson`)


Finally, clone this repository:

`$ git clone https://github.com/galantelab/sideRETRO.git`

Inside sideRETRO directory, run:

`$ meson build && ninja -C build`

You can find `sider` executable inside `build/src`. Optionally, install to system directories with:

`$ sudo ninja -C build install`

### Usage

**sideRETRO** compiles to an executable called `sider`, which has three subcommands: `process-sample`, `merge-call`
and `make-vcf`. The `process-sample` subcommand processes a list of **SAM/BAM/CRAM** files, and captures abnormal reads
that must be related to an event of retrocopy. All those data is saved to a **SQLite3 database** and then we come
to the second step `merge-call`, which processes the database and annotates all the retrocopies found. Finally we
can run the subcommand `make-vcf` and generate a file (in **VCF** format) with retrocopies and further information
about them.

```sh
# List of BAM files
$ cat 'my-bam-list.txt'
/path/to/file1.bam
/path/to/file2.bam
/path/to/file3.bam

# Run process-sample step
$ sider process-sample \
    --annotation-file='my-annotation.gtf' \
    --input-file='my-bam-list.txt'

$ ls -1
my-genome.fa
my-annotation.gtf
my-bam-list.txt
out.db

# Run merge-call step
$ sider merge-call --in-place out.db

# Run make-vcf step
$ sider make-vcf \
    --reference-file='my-genome.fa' out.db
```

## Citation

If sideRETRO was somehow useful in your research, please cite it:

```bib
@article{10.1093/bioinformatics/btaa689,
  author = {Miller, Thiago L A and Orpinelli, Fernanda and Buzzo, José Leonel L and Galante, Pedro A F},
  title = "{sideRETRO: a pipeline for identifying somatic and polymorphic insertions of processed pseudogenes or retrocopies}",
  journal = {Bioinformatics},
  year = {2020},
  month = {07},
  issn = {1367-4803},
  doi = {10.1093/bioinformatics/btaa689},
  url = {https://doi.org/10.1093/bioinformatics/btaa689},
  note = {btaa689},
}
```

## License

This is free software, licensed under:

`The GNU General Public License, Version 3, June 2007`

