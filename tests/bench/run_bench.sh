#!/usr/bin/env bash
#
# sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
# Copyright (C) 2019-2026 Thiago L. A. Miller <tmiller@mochsl.org.br>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Run the SANDY benchmark: Design the events, simulate a whole genome
# with them ('events') and without them ('control') using the very
# same parameters, align with bwa-mem, run sideRETRO and evaluate.
#
# The reads are simulated from the primary chromosomes (chr1-22, chrX:
# a diploid female genome) and aligned to the full reference (alt,
# decoy and HLA contigs included), as a real sample.
#
# Every step is skipped if its output exists, so a broken run can be
# resumed. Settings come from the environment:
#
#   REF        Full reference FASTA for bwa and sideRETRO (bwa indexed)
#   GTF        GENCODE annotation
#   ALU, L1    Mobile element consensus FASTA
#   OUT        Output directory              [bench]
#   SEED       Seed for the design and SANDY [17]
#   COVERAGE   SANDY coverage                [30]
#   THREADS    Threads                       [8]
#   PHRED      process-sample -Q             [20]
#   ALIGN_FRAC process-sample -F             [0.5]
#   RUNS       Runs to do                    ["events control"]
#
# Example:
#
#   $ REF=hg38.fa GTF=gencode.gtf ALU=ALU.fa L1=LINE1.fa OUT=bench \
#       tests/bench/run_bench.sh

set -euo pipefail

: "${REF:?Set REF to the reference FASTA}"
: "${GTF:?Set GTF to the GENCODE annotation}"
: "${ALU:?Set ALU to the Alu consensus FASTA}"
: "${L1:?Set L1 to the L1 consensus FASTA}"
OUT=${OUT:-bench}
SEED=${SEED:-17}
COVERAGE=${COVERAGE:-30}
THREADS=${THREADS:-8}
PHRED=${PHRED:-20}
ALIGN_FRAC=${ALIGN_FRAC:-0.5}
RUNS=${RUNS:-"events control"}

BENCH=$(cd "$(dirname "$0")" && pwd)
VARIATION=bench_events_${SEED}

log() { echo "[$(date '+%F %T')] $*" >&2; }

mkdir -p "$OUT"

# Primary chromosomes for SANDY
PRIMARY_FA=$OUT/primary.fa
if [ ! -s "$PRIMARY_FA" ]; then
	log "Extract the primary chromosomes"
	samtools faidx "$REF" chr{1..22} chrX > "$PRIMARY_FA"
	samtools faidx "$PRIMARY_FA"
fi

# Events and truth
if [ ! -s "$OUT/events.tsv" ]; then
	log "Design the events"
	python3 "$BENCH/sandy_events.py" --reference "$PRIMARY_FA" \
		--annotation "$GTF" --alu "$ALU" --l1 "$L1" --seed "$SEED" \
		--output-dir "$OUT"
	sandy variation add -a "$VARIATION" "$OUT/events.tsv"
fi

for run in $RUNS; do
	dir=$OUT/$run
	mkdir -p "$dir"

	# Simulate: The control differs only by the variations
	if ! ls "$dir"/reads*R1*.fastq.gz > /dev/null 2>&1; then
		log "[$run] Simulate the reads"
		var_opt=()
		[ "$run" = events ] && var_opt=(-a "$VARIATION")
		sandy genome --verbose "${var_opt[@]}" \
			--seed="$SEED" --jobs="$THREADS" --coverage="$COVERAGE" \
			--sequencing-type=paired-end --read-mean=150 \
			--fragment-mean=500 --fragment-stdd=50 \
			--output-format=fastq.gz --output-dir="$dir" --prefix=reads \
			"$PRIMARY_FA"
	fi

	# Align
	if [ ! -s "$dir/reads.bam.bai" ]; then
		log "[$run] Align"
		r1=$(ls "$dir"/reads*R1*.fastq.gz)
		r2=$(ls "$dir"/reads*R2*.fastq.gz)
		bwa mem -t "$THREADS" -R "@RG\tID:$run\tSM:$run" "$REF" "$r1" "$r2" \
			| samtools sort -@ "$THREADS" -o "$dir/reads.bam" -
		samtools index "$dir/reads.bam"
	fi

	# sideRETRO
	if [ ! -s "$dir/out.vcf" ]; then
		log "[$run] Run sideRETRO"
		sider process-sample --annotation-file="$GTF" --output-dir="$dir" \
			--prefix=out --gff-feature=exon \
			--gff-attribute='transcript_type=protein_coding' \
			--gff-gene-name=gene_name --gff-exon-id=exon_id \
			--phred-quality="$PHRED" --alignment-frac="$ALIGN_FRAC" \
			--threads="$THREADS" --log-file="$dir/ps.log" --silent \
			"$dir/reads.bam"
		sider merge-call --reference-file="$REF" --output-dir="$dir" \
			--prefix=out --in-place --threads="$THREADS" \
			--log-file="$dir/mc.log" --silent "$dir/out.db"
		sider make-vcf --reference-file="$REF" --output-dir="$dir" \
			--prefix=out --log-file="$dir/vcf.log" --silent "$dir/out.db"
	fi

	# Evaluate: For the control, every call is a false positive
	log "[$run] Evaluate"
	python3 "$BENCH/bench_eval.py" --prefix "$dir/eval" "$OUT/truth.tsv" \
		"$dir/out.vcf" | tee "$dir/eval.txt"
done
