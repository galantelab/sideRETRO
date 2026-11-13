# Simulation

`retrocopy_sim.py` generates a reproducible retrocopy insertion for
testing the whole pipeline (`process-sample`, `merge-call` and
`make-vcf`) without running an aligner. It writes the reference
genome, the parental gene annotation, the alignments in SAM format
and a `truth.tsv` with the expected values. It also compares the
resulting VCF with the truth.

```sh
$ tests/sim/retrocopy_sim.py --seed 11 --output-dir sim
$ sider process-sample -s -a sim/genes.gtf -o sim sim/sim.sam
$ sider merge-call -r sim/ref.fa -o sim -p merged sim/out.db
$ sider make-vcf -r sim/ref.fa -o sim -p merged sim/merged.db
$ tests/sim/retrocopy_sim.py --check sim/merged.vcf sim/truth.tsv
OK: sim/merged.vcf matches sim/truth.tsv
```

Run `retrocopy_sim.py --help` for the simulation parameters: TSD
length, poly-A length, number of discordant and split reads, etc.

`run_sim_test.py` runs a set of scenarios through the pipeline, with
and without the reference genome at `merge-call`, and is registered
at `meson test` as the `simulation` suite:

```sh
$ meson test -C builddir --suite simulation
```
