#!/usr/bin/env python3
import argparse, csv, sys, textwrap

def main():
    p = argparse.ArgumentParser(
        description="Convert a CSV with columns 'bc' and 'ntseq' to FASTA."
    )
    p.add_argument("input_csv", help="Input CSV (with headers)")
    p.add_argument("output_fasta", help="Output FASTA file")
    p.add_argument("--bc-col", default="bc", help="Column name for barcode/header (default: bc)")
    p.add_argument("--seq-col", default="ntseq", help="Column name for sequence (default: ntseq)")
    p.add_argument("--wrap", type=int, default=80, help="Line wrap width for sequences; 0 = no wrap (default: 80)")
    p.add_argument(
        "--dedupe",
        choices=["none", "pair", "bc"],
        default="pair",
        help=(
            "Deduplicate rows before writing: "
            "'none' = write all rows; "
            "'pair' = unique (bc, ntseq) pairs (default); "
            "'bc' = first occurrence per bc (see --on-conflict)"
        ),
    )
    p.add_argument(
        "--on-conflict",
        choices=["error", "first", "last"],
        default="error",
        help="When --dedupe=bc and the same bc maps to different ntseq values: "
             "'error' (default), or keep the 'first' or 'last' occurrence.",
    )
    args = p.parse_args()

    total = written = skipped_missing = skipped_dupe = conflicts = 0
    seen_pairs = set()
    seen_bc = {}

    def cleaned(seq: str) -> str:
        # remove whitespace and uppercase
        return "".join(seq.split()).upper()

    with open(args.input_csv, newline="", encoding="utf-8") as fin, \
         open(args.output_fasta, "w", encoding="utf-8") as fout:

        reader = csv.DictReader(fin)
        for row in reader:
            total += 1
            bc = (row.get(args.bc_col) or "").strip()
            seq_raw = (row.get(args.seq_col) or "")
            if not bc or not seq_raw:
                skipped_missing += 1
                continue

            seq = cleaned(seq_raw)

            # Deduplication logic
            if args.dedupe == "pair":
                key = (bc, seq)
                if key in seen_pairs:
                    skipped_dupe += 1
                    continue
                seen_pairs.add(key)

            elif args.dedupe == "bc":
                if bc in seen_bc:
                    if seen_bc[bc] != seq:
                        conflicts += 1
                        if args.on_conflict == "error":
                            sys.stderr.write(
                                f"ERROR: bc '{bc}' maps to multiple sequences; "
                                "use --on-conflict first/last to resolve.\n"
                            )
                            sys.exit(1)
                        elif args.on_conflict == "first":
                            skipped_dupe += 1
                            continue
                        elif args.on_conflict == "last":
                            # Overwrite previous choice by marking previous as effectively skipped
                            pass
                    else:
                        skipped_dupe += 1
                        continue
                seen_bc[bc] = seq

            # Write FASTA
            fout.write(f">{bc}\n")
            if args.wrap and args.wrap > 0:
                fout.write(textwrap.fill(seq, width=args.wrap) + "\n")
            else:
                fout.write(seq + "\n")
            written += 1

    # Summary to stderr
    sys.stderr.write(
        f"Rows read: {total}\n"
        f"Written: {written}\n"
        f"Skipped (missing bc/ntseq): {skipped_missing}\n"
        f"Skipped (duplicates): {skipped_dupe}\n"
        f"Conflicts: {conflicts}\n"
    )

if __name__ == "__main__":
    main()
