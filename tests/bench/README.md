# SANDY benchmark

A whole genome simulation with known events, to measure the precision,
recall and F1 of sideRETRO and to tune its options and filters.

| file              | role                                                   |
|-------------------|--------------------------------------------------------|
| `sandy_events.py` | Design the events: SANDY TSV, truth table, config      |
| `run_bench.sh`    | Design, simulate (events and control), align, call, evaluate |
| `bench_eval.py`   | Compare a sideRETRO VCF to the truth table             |

## Events

The default mix (`DEFAULT_CONFIG` at `sandy_events.py`, overridden by a
JSON file with `--config`):

- 30 retrocopies, the true positives: spliced transcripts of protein
  coding genes, 5' truncated, with a poly-A tail, at L1 endonuclease
  motifs. Target sites: 20 DUP (7-20 bp), 3 BLUNT, 4 short DEL
  (1-50 bp) and 3 long DEL (1-20 kb). 3 single exon parental genes,
  3 without poly-A, 15 intergenic, 12 intronic and 3 exonic.
- Decoys, which must not be called: 30 Alu, 20 L1 (5' truncated),
  10 CNV (5 tandem duplications and 5 deletions spanning exons) and
  10 unprocessed copies (genomic gene segments with introns).

SANDY replaces the reference sequence by the alternate one, so each
insertion carries its target site:

| target site | REF       | ALT       | LBKP    | RBKP        | TSDLEN |
|-------------|-----------|-----------|---------|-------------|--------|
| DUP of L    | `S`       | `S R S`   | cut + L | cut + 1     | L      |
| BLUNT       | `b`       | `b R`     | cut     | cut + 1     | 0      |
| DEL of d    | `b D`     | `b R`     | cut     | cut + d + 1 | -d     |

where `cut` is the 0-based position of the cut, `S` the L bases from
it, `b` the base before it and `R` the insertion.

## Run

```
$ REF=GRCh38_full_analysis_set_plus_decoy_hla.fa \
  GTF=gencode.annotation.gtf ALU=ALU.fa L1=LINE1.fa OUT=builddir/bench \
  tests/bench/run_bench.sh
```

The reads are simulated from chr1-22 and chrX and aligned to the full
reference. The control run has the very same parameters without the
events: All its calls are false positives.
