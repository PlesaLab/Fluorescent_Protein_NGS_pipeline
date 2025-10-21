# High-Diversity FP Libraries — NGS Analysis Pipeline

**Manuscript:** Benabbas *et al.*, *High diversity gene libraries facilitate machine learning–guided exploration of fluorescent protein sequence space* (2025 preprint)

This repository contains scripts and a one-shot bash runner to reproduce the NGS processing for the **Parents Codon 1 (C1P)** and **Parents Codon 2 (C2P)** libraries.

This pipeline assumes you created your amplicons using primers with the following homology to pEVBC1:

* FWD:

 *TGTGAGCGGATAACAATTTCACACAGGAAACAGCTCATATG*
* REV:

 *CGAAAAGTGCCACCTGACGTCGTGC*


---

## Contents
- `run_pipeline.sh` — end-to-end pipeline runner
- `parse_SAM_alignments_MP.py` — parses SAM alignments and assigns `mutID`s (CLI w/ `--lib {codon1,codon2}`)
- `barcode_extract.py` — extracts 20-nt barcodes (fuzzy motif) and inserts between NdeI/KpnI; writes CSV/TSV files
- `collision.py` — Python consensus caller (replaces `collision.R`) implementing staged tie resolution
- `translate.toUpper.v2.aatrim.py` — translates FASTA (uppercase, trims at first `*`)
- Reference files:
  - `FP.C1.v2.genes` — DNA reference for C1
  - `FP.C2.v2.genes` — DNA reference for C2
  - `FP_12.proteins` — protein reference for `parse_SAM_alignments.py`

---

## Data
Raw PacBio reads (CCS) used in the manuscript:

- **C1P (Parents Codon 1):** SRA accession [SRX29434776](https://www.ncbi.nlm.nih.gov/sra/SRX29434776) → save as `C1P.6291.fastq.gz`

- **C2P (Parents Codon 2):** SRA accession [SRX29434777](https://www.ncbi.nlm.nih.gov/sra/SRX29434777) → save as `C2P.6291.fastq.gz`

---

## Requirements

### Binaries
- [STARcode](https://github.com/gui11aume/starcode) ≥ 1.3
- [BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) ≥ 39
- GNU coreutils, bash, gzip

### Python
- Python ≥ 3.9
- Packages: `biopython`, `hashlib `,`multiprocessing`,`pandas `
- use pip or conda to install

---

## Quickstart

1. Place `C1P.6291.fastq.gz` and `C2P.6291.fastq.gz` in the repo root.
2. Ensure `barcode_extract.py`, `collision.py`, `all_info_csv_to_fasta.py`, `parse_SAM_alignments.py`, `translate.toUpper.v2.aatrim.py`, and references are present.
3. Run:
   ```bash
   ./run_pipeline.sh
   ```

All outputs are written under `analysis_output/` unless noted.

---

## Pipeline Outline

### 1) Extract barcode and gene (Python)
```bash
python barcode_extract.py  # writes parents.C1.csv, parents.C1_noN_aa.csv, parents.C1_bc_stats*.{csv,tsv}, parents.C1_bc_list.csv (and C2 analogs)
```
**barcode_extract.py details**

- Fuzzy motif: `((TGGCTGCGGAAC)(....................)(GCACGACGTCAG)){e<4}` (allows up to 3 mismatches) to locate the 20-nt barcode.
- Trims the insert between **NdeI** (`CATATG`) and **KpnI** (`GGTACC` in this context searched with stop codon too: `TAAGGTACC`).
- Outputs:
  - `parents.C1.csv` / `parents.C2.csv` (BC→insert)
  - `parents.C1_noN_aa.csv` / `parents.C2_noN_aa.csv` (no 'N', plus AA translation)
  - `parents.C1_bc_stats.csv` / `parents.C2_bc_stats.csv` and `*_bc_stats_for_starcode.tsv`
  - `parents.C1_bc_list.csv` / `parents.C2_bc_list.csv`

### 2) Collapse barcodes (STARcode, LD=1, spheres)
```bash
starcode -d 1 -t 2 --print-clusters -i parents.C1bc_stats_for_starcode.tsv -o parents.C1_clustered_LD1.txt >> C1.sc.log 2>&1
starcode -d 1 -t 2 --print-clusters -i parents.C2bc_stats_for_starcode.tsv -o parents.C2_clustered_LD1.txt >> C2.sc.log 2>&1
```

### 3) Build barcode consensus (Python)

```bash
python collision.py
```

**collision.py details**

* Extracts inserts between **NdeI** (`CATATG`) and **KpnI** (`GGTACC`); rescues reads with multiple NdeI but one KpnI (takes last NdeI/KpnI).
* Determines per-barcode consensus using staged tie resolution:

  1. Majority reads
  2. Tie-resolved by identical translation
  3. Tie-resolved by longest trimmed translation
  4. Tie-resolved by most common trimmed translation
* Outputs:

  * `parents.C?_consensus_summary_v2.csv` — BC, consensus\_ntseq, total\_reads, consensus_call
  * `parents.C?_all_info_v2.csv` — final BC → ntseq + AA fields

### 4) Convert to fasta (Python)

convert the csv to fasta:

```bash
python all_info_csv_to_fasta.py parents.C1_all_info_v2.csv parents.C1_all_info.fasta --wrap 0 

python all_info_csv_to_fasta.py parents.C2_all_info_v2.csv parents.C2_all_info.fasta --wrap 0 
```

alternatively you could also use awk

```bash
awk -F, 'NR>1{print ">"$1"\n"$2}' parents.C1_all_info_v2.csv > parents.C1_all_info.fasta
awk -F, 'NR>1{print ">"$1"\n"$2}' parents.C2_all_info_v2.csv > parents.C2_all_info.fasta
```

### 5) Align at DNA level (BBMap)

Run:

```bash
bbmap.sh ref=FP.C1.v2.genes
bbmap.sh t=4 in=parents.C1_all_info.fasta outm=parents.C1.bc-nt.map.sam outu=parents.C1.bc-nt.map.unaligned.sam 

bbmap.sh ref=FP.C2.v2.genes
bbmap.sh t=4 in=parents.C2_all_info.fasta outm=parents.C2.bc-nt.map.sam outu=parents.C2.bc-nt.map.unaligned.sam 
```

### 6) Parse SAM + assign mutant IDs (Python)
```bash
python parse_SAM_alignments_MP_cli.py --lib codon1
python parse_SAM_alignments_MP_cli.py --lib codon2
```

### 7) Translate consensus sequences (Python)
```bash
python translate.toUpper.v2.aatrim.py parents.C1_all_info.fasta > parents.C1_all_info.trans.fasta
python translate.toUpper.v2.aatrim.py parents.C2_all_info.fasta > parents.C2_all_info.trans.fasta
```

---

## Expected Outputs (key files)

- `parents.C1.bc-nt.map.sam`, `parents.C2.bc-nt.map.sam` — nucleotide alignments (SAM)
- `analysis_output/C1seqs_mutID_all.csv`, `analysis_output/C1seqs_mutID_info_all.csv`
- `analysis_output/C2seqs_mutID_all.csv`, `analysis_output/C2seqs_mutID_info_all.csv`
- `parents.C1_all_info.trans.fasta`, `parents.C2_all_info.trans.fasta`
- `parents.C1_consensus_summary_v2.csv`, `parents.C2_consensus_summary_v2.csv`
- `parents.C1_all_info_v2.csv`, `parents.C2_all_info_v2.csv`

---

## Reproducibility Tips

- Use `--threads` in `parse_SAM_alignments_MP_cli.py` to increase throughput.
