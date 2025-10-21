#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# High-Diversity FP Libraries — One-shot NGS Pipeline
# -----------------------------
# This script executes the full analysis for C1P and C2P as used in the manuscript.
# It assumes the following input files exist in the working directory:
#   C1P.6291.fastq.gz
#   C2P.6291.fastq.gz
#
# It also assumes the following scripts and references are present:
#   barcode_extract.py
#   collision.py
#   translate.toUpper.v2.aatrim.py
#   FP.C1.v2.genes
#   FP.C2.v2.genes
#   parse_SAM_alignments.py
#
# Outputs are written to ./analysis_output/
#
# Tunables
THREADS_STARCODE=2
THREADS_BBMAP=4
PY_THREADS_PARSE=10

# ---- helper for checking dependencies
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: required tool '$1' not found in PATH" >&2; exit 127; }; }

echo "[check] dependencies"
need starcode
need bbmap.sh
need python

echo "[check] scripts and references"
for f in barcode_extract.py collision.R translate.toUpper.v2.aatrim.py FP.C1.v2.genes FP.C2.v2.genes parse_SAM_alignments.py; do
  [[ -f "$f" ]] || { echo "ERROR: missing required file: $f" >&2; exit 2; }
done

mkdir -p analysis_output

echo "[0/7] Check for input files"
[[ -f C1P.6291.fastq.gz ]] || { echo "ERROR: missing C1P.6291.fastq.gz" >&2; exit 3; }
[[ -f C2P.6291.fastq.gz ]] || { echo "ERROR: missing C2P.6291.fastq.gz" >&2; exit 3; }

echo "[1/7] Extract barcode + gene (Python)"
python barcode_extract.py

echo "[2/7] Collapse barcodes with STARcode (LD=1, spheres)"
starcode -d 1 -t "${THREADS_STARCODE}" --print-clusters \
  -i parents.C1bc_stats_for_starcode.tsv \
  -o parents.C1_clustered_LD1.txt >> C1.sc.log 2>&1
starcode -d 1 -t "${THREADS_STARCODE}" --print-clusters \
  -i parents.C2bc_stats_for_starcode.tsv \
  -o parents.C2_clustered_LD1.txt >> C2.sc.log 2>&1

echo "[3/7] Build consensus per barcode (Python)"
python collision.py

echo "[4/7] Convert csv to fasta"
python all_info_csv_to_fasta.py parents.C1_all_info_v2.csv parents.C1_all_info.fasta --wrap 0
python all_info_csv_to_fasta.py parents.C2_all_info_v2.csv parents.C2_all_info.fasta --wrap 0
#awk -F, 'NR>1{print ">"$1"\n"$2}' parents.C1_all_info_v2.csv > parents.C1_all_info.fasta
#awk -F, 'NR>1{print ">"$1"\n"$2}' parents.C2_all_info_v2.csv > parents.C2_all_info.fasta

echo "[5/7] BBMap — DNA alignment (C1 then C2)"
bbmap.sh ref=FP.C1.v2.genes
bbmap.sh t="${THREADS_BBMAP}" in=parents.C1_all_info.fasta \
  outm=parents.C1.bc-nt.map.sam outu=parents.C1.bc-nt.map.unaligned.sam
bbmap.sh ref=FP.C2.v2.genes
bbmap.sh t="${THREADS_BBMAP}" in=parents.C2_all_info.fasta \
  outm=parents.C2.bc-nt.map.sam outu=parents.C2.bc-nt.map.unaligned.sam

echo "[6/7] Parse SAM + assign mutIDs (Python)"
python parse_SAM_alignments.py --lib codon1 --threads "${PY_THREADS_PARSE}"
python parse_SAM_alignments.py --lib codon2 --threads "${PY_THREADS_PARSE}"

echo "[7/7] Translate consensus sequences"
python translate.toUpper.v2.aatrim.py parents.C1_all_info.fasta > parents.C1_all_info.trans.fasta
python translate.toUpper.v2.aatrim.py parents.C2_all_info.fasta > parents.C2_all_info.trans.fasta

echo "[done] All steps completed. Key outputs in ./analysis_output/"
