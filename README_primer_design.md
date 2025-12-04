# Dial-out PCR Primer Design

This tool designs forward and reverse primers to specifically retrieve genes from barcode-tagged libraries using dial-out PCR.

## Inputs

One or more CSVs (one per library), each with four columns:

- `BC` — barcode sequence on the top strand (usually 20 nt).
- `consensus_ntseq` — variable gene region (~675 bp).
- `total_reads`
- `consensus_call`

Examples:

- `parents.C1_consensus_summary_v2py.csv`
- `parents.C2_consensus_summary_v2py.csv`

## Plasmid Context (Top Strand)

```
... UPSTREAM ... agctCATATG + [consensus_ntseq] + GGTACCtaaGTGTGGCTGCGGAAC + [BC] + gcacGACGTcaggtggcacttttcggggaaatgtgcgcggaacccctatt
```

Where:
- UPSTREAM (full in code) ends with `agctCATATG` (includes NdeI site).
- `consensus_ntseq` ≈ 675 bp, varies per gene.
- `BC` is random, typically 20 bp.
- The conserved suffix after BC is `gcacGACGTcaggtggcacttttcggggaaatgtgcgcggaacccctatt`.

## Primer Design Rules Implemented

### Forward primer
- Must **begin with** the conserved motif `agctCATATG` at its 5' end.
- 3' end lies **inside the gene** (consensus_ntseq).
- Includes **>= 6 nt** of the gene region.
- **Length >= 20 nt**. The script extends into the gene to hit length/Tm targets.

### Reverse primer
- **3' end fixed**: anneals to the **start of BC** (top strand).
- **Length >= 20 nt**. If needed, it extends by adding conserved sequence **after BC**.
- Built as reverse complement of the first L bases of (`BC + conserved_suffix_after_BC`).

## Filtering/Checks

For each gene, the script proposes up to **3 primer pairs** that pass all filters:

- **Tm >= 58 degC** (primer3; salt/dna conc configurable via CLI).
- **Hairpins**: reports dG (kcal/mol) from **seqfold** if available; otherwise falls back to primer3 hairpin dG (and notes method).
- **Homodimers**: no more than **4** contiguous complementary bases at the **3' end** (3'-anchored) for self-dimers.
- **Heterodimers (pair)**: same 3'-end <= 4 rule across the pair.

The 3'-end dimer rule is implemented via an anti-parallel alignment that finds the longest 3'-anchored contiguous complement in each direction.

## Output CSV

Original columns plus:

- `source_file` (which input CSV the row came from)
- `pair_index` (1..N per gene)
- `forward_primer`
- `reverse_primer`
- `forward_tm` (degC, primer3)
- `reverse_tm` (degC, primer3)
- `forward_hairpin_dG` (kcal/mol; seqfold if available, else primer3)
- `reverse_hairpin_dG` (kcal/mol)
- `forward_hairpin_method` (`seqfold` | `primer3` | `none`)
- `reverse_hairpin_method`
- `forward_homodimer_3p_run` (max contiguous 3' complement)
- `reverse_homodimer_3p_run`
- `heterodimer_3p_run`
- `note` (e.g., “No primer pair met filters”)

If no pair passes, a single row with empty primer fields is emitted with `note` explaining why.

## Install

```bash
python3 -m venv venv
source venv/bin/activate
pip install pandas biopython primer3-py seqfold
```

## Run

```bash
python primer_design.py --threads 4 \
   --pairs-per-gene 3 \
   --input parents.C1_consensus_summary_v2py.csv \
   parents.C2_consensus_summary_v2py.csv
```

Outputs:

- `parents.C1_consensus_summary_v2py_dialout.csv`
- `parents.C2_consensus_summary_v2py_dialout.csv`


### Optional flags

- `--tm-min` (default 58.0)
- `--mv-conc` (mM, default 50.0)
- `--dv-conc` (mM, default 1.5)
- `--dntp-conc` (mM, default 0.6)
- `--dna-conc` (uM, default 0.25)
- `--temp-c` (degC for hairpin eval, default 60.0)
- `--max-pairs` (default 3)
- `--f-len-min` `--f-len-max` (default 20, 34)
- `--r-len-min` `--r-len-max` (default 20, 34)
- `--max-3p-overlap` (default 4)
- `--pairs-per-gene` (default 1) controls how many passing primer pairs to keep per gene.
- `--threads` (default **4**) parallelizes per-row designs for faster runs.

## Assumptions

- Forward primer strictly **starts with** `agctCATATG` and only extends **into the gene** (not further upstream).

- Reverse primer always has its 3' terminal base mapped to the **first base of BC**; length is modulated by including more of the downstream conserved region on the template side.

- If `seqfold` API is unavailable, the script falls back to primer3 hairpin dG and records the method used per primer.
