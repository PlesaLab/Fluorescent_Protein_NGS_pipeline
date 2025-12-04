#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Dial-out PCR primer design with progress and multithreading.
- Progress: uses tqdm if available; otherwise prints periodic updates.
- Multithreading: --threads (default 4) parallelizes per-row design work.
- Per-file outputs: writes <basename>_dialout.csv for each input.
- Pairs per gene: --pairs-per-gene (default 1).

Dependencies:
  pip install pandas biopython primer3-py seqfold tqdm
(tqdm is optional but recommended; without it you'll get periodic prints.)
"""

import argparse
import sys
from typing import List, Dict, Optional, Any, Tuple
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

from Bio.Seq import Seq

# Optional progress bar
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except Exception:
    TQDM_AVAILABLE = False

try:
    import primer3
except ImportError as e:
    sys.stderr.write("ERROR: primer3-py is required. Install with: pip install primer3-py\n")
    raise

SEQFOLD_AVAILABLE = True
try:
    import seqfold
except Exception:
    SEQFOLD_AVAILABLE = False

UPSTREAM_FULL = "gtgtggaattgtgagcggataacaatttcacacaggaaacagct"
UPSTREAM_MUST_PREFIX = "agctCATATG"
DOWNSTREAM_AFTER_GENE = "GGTACCtaaGTGTGGCTGCGGAAC"
SUFFIX_AFTER_BC = "gcacGACGTcaggtggcacttttcggggaaatgtgcgcggaacccctatt"

UPSTREAM_FULL = UPSTREAM_FULL.upper()
UPSTREAM_MUST_PREFIX = UPSTREAM_MUST_PREFIX.upper()
DOWNSTREAM_AFTER_GENE = DOWNSTREAM_AFTER_GENE.upper()
SUFFIX_AFTER_BC = SUFFIX_AFTER_BC.upper()

def revcomp(s: str) -> str:
    return str(Seq(s).reverse_complement())

def calc_tm(seq: str, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=0.25) -> float:
    return primer3.calc_tm(
        seq,
        mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc,
    )

def hairpin_energy_seqfold(seq: str, temp_c: float = 60.0) -> Optional[float]:
    if not SEQFOLD_AVAILABLE:
        return None
    try:
        stems = seqfold.fold(seq, temp_c)
        if not stems:
            return 0.0
        return min(stem.dg for stem in stems)
    except Exception:
        return None

def hairpin_energy_primer3(seq: str, temp_c: float = 60.0) -> Optional[float]:
    try:
        hp = primer3.calc_hairpin(seq, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=0.25, temp_c=temp_c)
        if hp.structure_found:
            return hp.dg / 1000.0
        else:
            return 0.0
    except Exception:
        return None

def max_3prime_contig_complementarity(seq5: str, seq3: str) -> int:
    def longest_3prime_anchor(a, b):
        rb = revcomp(b)
        max_run = 0
        for k in range(1, min(len(a), len(rb)) + 1):
            seg_a = a[-k:]
            seg_rb = rb[:k]
            run = 0
            for i in range(1, k + 1):
                if seg_a[-i] == seg_rb[-i]:
                    run += 1
                else:
                    break
            if run > max_run:
                max_run = run
        return max_run
    return max(longest_3prime_anchor(seq5, seq3), longest_3prime_anchor(seq3, seq5))

def passes_3prime_rule_homodimer(seq: str, max_allowed: int = 4) -> bool:
    return max_3prime_contig_complementarity(seq, seq) <= max_allowed

def passes_3prime_rule_heterodimer(seq_f: str, seq_r: str, max_allowed: int = 4) -> bool:
    return max_3prime_contig_complementarity(seq_f, seq_r) <= max_allowed

def build_forward_candidates(consensus_ntseq: str, min_len: int = 20, max_len: int = 34) -> List[str]:
    prefix = UPSTREAM_MUST_PREFIX
    seq = consensus_ntseq.upper()
    min_L = max(min_len, len(prefix) + 6)
    cands = []
    for L in range(min_L, max_len + 1):
        k = L - len(prefix)
        if k <= 0:
            continue
        if k > len(seq):
            break
        primer = (prefix + seq[:k]).upper()
        cands.append(primer)
    return cands

def build_reverse_candidates(bc: str, min_len: int = 20, max_len: int = 34) -> List[str]:
    template = (bc.upper() + SUFFIX_AFTER_BC).upper()
    cands = []
    for L in range(min_len, max_len + 1):
        if L > len(template):
            break
        primer = revcomp(template[:L])
        cands.append(primer)
    return cands

def score_and_filter_pairs(forward_list: List[str], reverse_list: List[str],
                           tm_min: float = 58.0,
                           mv_conc: float = 50.0, dv_conc: float = 1.5, dntp_conc: float = 0.6, dna_conc: float = 0.25,
                           temp_c: float = 60.0,
                           max_three_prime_overlap: int = 4,
                           max_pairs: int = 1) -> List[Dict]:
    results = []
    def primer_metrics(pseq: str) -> Dict:
        tm = calc_tm(pseq, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc)
        hp_seqfold = hairpin_energy_seqfold(pseq, temp_c=temp_c) if SEQFOLD_AVAILABLE else None
        hp_fallback = hairpin_energy_primer3(pseq, temp_c=temp_c) if hp_seqfold is None else None
        hp_used = hp_seqfold if hp_seqfold is not None else hp_fallback
        hom_ok = passes_3prime_rule_homodimer(pseq, max_three_prime_overlap)
        hom_score = max_3prime_contig_complementarity(pseq, pseq)
        return {
            "tm": round(tm, 2),
            "hairpin_dG_kcalmol": None if hp_used is None else round(hp_used, 2),
            "hairpin_method": "seqfold" if (hp_seqfold is not None) else ("primer3" if hp_fallback is not None else "none"),
            "homodimer_3p_run": hom_score,
            "homodimer_pass": hom_ok
        }

    f_cache, r_cache = {}, {}
    for f in forward_list:
        f_cache[f] = primer_metrics(f)
    for r in reverse_list:
        r_cache[r] = primer_metrics(r)

    for f in forward_list:
        fm = f_cache[f]
        if fm["tm"] < tm_min or not fm["homodimer_pass"]:
            continue
        for r in reverse_list:
            rm = r_cache[r]
            if rm["tm"] < tm_min or not rm["homodimer_pass"]:
                continue
            het_ok = passes_3prime_rule_heterodimer(f, r, max_three_prime_overlap)
            het_score = max_3prime_contig_complementarity(f, r)
            if not het_ok:
                continue
            results.append({
                "forward_primer": f,
                "reverse_primer": r,
                "forward_tm": fm["tm"],
                "reverse_tm": rm["tm"],
                "forward_hairpin_dG": fm["hairpin_dG_kcalmol"],
                "reverse_hairpin_dG": rm["hairpin_dG_kcalmol"],
                "forward_hairpin_method": fm["hairpin_method"],
                "reverse_hairpin_method": rm["hairpin_method"],
                "forward_homodimer_3p_run": fm["homodimer_3p_run"],
                "reverse_homodimer_3p_run": rm["homodimer_3p_run"],
                "heterodimer_3p_run": het_score
            })
            if len(results) >= max_pairs:
                return results
    return results[:max_pairs]

def process_row(row_dict, tm_min, mv_conc, dv_conc, dntp_conc, dna_conc, temp_c,
                max_pairs, f_len_min, f_len_max, r_len_min, r_len_max, max_three_prime_overlap):
    bc = str(row_dict["BC"]).strip().upper()
    gene = str(row_dict["consensus_ntseq"]).strip().upper()

    fwd_cands = build_forward_candidates(gene, min_len=f_len_min, max_len=f_len_max)
    rev_cands = build_reverse_candidates(bc, min_len=r_len_min, max_len=r_len_max)

    pairs = score_and_filter_pairs(
        fwd_cands, rev_cands,
        tm_min=tm_min,
        mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc,
        temp_c=temp_c,
        max_three_prime_overlap=max_three_prime_overlap,
        max_pairs=max_pairs
    )

    if not pairs:
        return [{
            **row_dict,
            "pair_index": 1,
            "forward_primer": "",
            "reverse_primer": "",
            "forward_tm": "",
            "reverse_tm": "",
            "forward_hairpin_dG": "",
            "reverse_hairpin_dG": "",
            "forward_hairpin_method": "none",
            "reverse_hairpin_method": "none",
            "forward_homodimer_3p_run": "",
            "reverse_homodimer_3p_run": "",
            "heterodimer_3p_run": "",
            "note": "No primer pair met filters"
        }]

    out = []
    for i, p in enumerate(pairs, start=1):
        out.append({
            **row_dict,
            "pair_index": i,
            **p
        })
    return out

def process_file(path: str, tm_min: float, mv_conc: float, dv_conc: float, dntp_conc: float, dna_conc: float,
                 temp_c: float, max_pairs: int, f_len_min: int, f_len_max: int, r_len_min: int, r_len_max: int,
                 max_three_prime_overlap: int, threads: int, show_progress: bool) -> pd.DataFrame:
    df = pd.read_csv(path, sep=",")
    from pathlib import Path


    required = {"BC", "consensus_ntseq"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")

    inputs = [(i, df.iloc[i].to_dict()) for i in range(len(df))]
    results_rows = []

    if threads <= 1:
        iterator = inputs
        if show_progress and TQDM_AVAILABLE:
            iterator = tqdm(iterator, total=len(inputs), desc=f"Designing ({Path(path).name})")
        for i, rowd in iterator:
            rows = process_row(rowd, tm_min, mv_conc, dv_conc, dntp_conc, dna_conc, temp_c, max_pairs,
                               f_len_min, f_len_max, r_len_min, r_len_max, max_three_prime_overlap)
            results_rows.append((i, rows))
            if show_progress and not TQDM_AVAILABLE and (i+1) % 500 == 0:
                print(f"{Path(path).name}: processed {i+1}/{len(inputs)} rows...")
    else:
        with ThreadPoolExecutor(max_workers=threads) as ex:
            futures = {ex.submit(process_row, rowd, tm_min, mv_conc, dv_conc, dntp_conc, dna_conc, temp_c, max_pairs,
                                 f_len_min, f_len_max, r_len_min, r_len_max, max_three_prime_overlap): i
                       for i, rowd in inputs}
            if show_progress and TQDM_AVAILABLE:
                for fut in tqdm(as_completed(futures), total=len(futures), desc=f"Designing ({Path(path).name})"):
                    i = futures[fut]
                    rows = fut.result()
                    results_rows.append((i, rows))
            else:
                done = 0
                for fut in as_completed(futures):
                    i = futures[fut]
                    rows = fut.result()
                    results_rows.append((i, rows))
                    done += 1
                    if show_progress and not TQDM_AVAILABLE and done % 500 == 0:
                        print(f"{Path(path).name}: processed {done}/{len(inputs)} rows...")

    results_rows.sort(key=lambda x: x[0])
    flat = []
    for _, group in results_rows:
        flat.extend(group)
    return pd.DataFrame(flat)

def main():
    import argparse
    from pathlib import Path
    ap = argparse.ArgumentParser(description="Design dial-out PCR primers per gene using fixed context and BC.")
    ap.add_argument("--input", "-i", nargs="+", required=True,
                    help="Input CSV files (e.g., parents.C1_consensus_summary_v2py.csv parents.C2_consensus_summary_v2py.csv)")
    ap.add_argument("--pairs-per-gene", type=int, default=1, help="How many primer pairs to keep per gene (default 1)")
    ap.add_argument("--tm-min", type=float, default=58.0, help="Minimum primer Tm (degC)")
    ap.add_argument("--mv-conc", type=float, default=50.0, help="Monovalent salt concentration, mM")
    ap.add_argument("--dv-conc", type=float, default=1.5, help="Divalent cation concentration, mM")
    ap.add_argument("--dntp-conc", type=float, default=0.6, help="dNTP concentration, mM")
    ap.add_argument("--dna-conc", type=float, default=0.25, help="Primer strand concentration, uM")
    ap.add_argument("--temp-c", type=float, default=60.0, help="Temperature (degC) for hairpin estimates")
    ap.add_argument("--f-len-min", type=int, default=20, help="Forward primer minimum length")
    ap.add_argument("--f-len-max", type=int, default=34, help="Forward primer maximum length")
    ap.add_argument("--r-len-min", type=int, default=20, help="Reverse primer minimum length")
    ap.add_argument("--r-len-max", type=int, default=34, help="Reverse primer maximum length")
    ap.add_argument("--max-3p-overlap", type=int, default=4, help="Max allowed contiguous 3'-complementarity for (homo/hetero)dimers")
    ap.add_argument("--threads", type=int, default=4, help="Number of threads to use (default 4)")
    ap.add_argument("--no-progress", action="store_true", help="Disable progress indicators")

    args = ap.parse_args()
    show_progress = not args.no_progress

    for path in args.input:
        df = process_file(
            path,
            tm_min=args.tm_min,
            mv_conc=args.mv_conc, dv_conc=args.dv_conc, dntp_conc=args.dntp_conc, dna_conc=args.dna_conc,
            temp_c=args.temp_c,
            max_pairs=args.pairs_per_gene,
            f_len_min=args.f_len_min, f_len_max=args.f_len_max,
            r_len_min=args.r_len_min, r_len_max=args.r_len_max,
            max_three_prime_overlap=args.max_3p_overlap,
            threads=args.threads,
            show_progress=show_progress
        )
        inpath = Path(path)
        df.insert(0, "source_file", inpath.name)
        out_path = inpath.with_name(inpath.stem + "_dialout.csv")
        df.to_csv(out_path, index=False)
        print(f"Wrote {out_path} with {len(df)} rows. SEQFOLD_AVAILABLE={SEQFOLD_AVAILABLE}, TQDM_AVAILABLE={TQDM_AVAILABLE}")

if __name__ == "__main__":
    main()
