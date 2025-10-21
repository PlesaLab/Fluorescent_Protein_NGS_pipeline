#!/usr/bin/env python3
# barcode_extract_fastq.py
# ------------------
# Extracts 20-nt barcodes from PacBio reads in FASTQ.GZ format using a fuzzy-matched flanking motif.
# Ignores quality scores.
#
# Inputs:
#   - C1P.6291.fastq.gz
#   - C2P.6291.fastq.gz
#
# Outputs (per library):
#   - parents.C?.csv                        : barcode → insert (NdeI..KpnI)
#   - parents.C?_noN_aa.csv                : barcode → insert (no 'N') → AA translation
#   - parents.C?_bc_stats.csv              : barcode read counts
#   - parents.C?_bc_stats_for_starcode.tsv : barcode read counts (tab-sep; Starcode input)
#   - parents.C?_bc_list.csv               : barcode, full upstream sequence (for auditing)
#
# Notes:
# - Uses fuzzy motif ((TGGCTGCGGAAC)(....................)(GCACGACGTCAG)){e<4}
#   to locate 20-nt barcode; searches forward then reverse-complement.
# - Extracts inserts between NdeI (CATATG) and KpnI (TAAGGTACC).
# - Translation-based QC: frame check and premature stop detection.

import regex
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import collections
import csv

# Input FASTQ.GZ files
inputz = ["C1P.6291.fastq.gz", "C2P.6291.fastq.gz"]
outputz = ["parents.C1.fasta", "parents.C2.fasta"]

for i in range(len(inputz)):

    input_file = inputz[i]
    output_file = outputz[i]

    d = collections.defaultdict(list)
    dnum = dict()
    dnum3 = dict()
    outputseqs = []
    outputseqs_noN = []
    outputseqs_noN_aa = []
    doutputseqs_noN = dict()
    doutputseqs_noN_aa = dict()

    motif = r'((TGGCTGCGGAAC)(....................)(GCACGACGTCAG)){e<4}'  # fuzzy 20nt barcode flanked by 12bp motifs

    # --- Read FASTQ.GZ directly ---
    with gzip.open(input_file, "rt") as handle:
        records = list(SeqIO.parse(handle, "fastq"))

    print(f"{len(records)} total records in {input_file}")

    count_no_ndei = 0
    count_no_kpni = 0
    bafframecount = 0
    rc_count = 0
    goodaa = 0

    for idx, record in enumerate(records, 1):
        if idx % 50000 == 0:
            print("On record:", idx)

        seq = str(record.seq).upper()
        rc_flag = 0

        match = regex.search(motif, seq, regex.BESTMATCH)
        if match is None:
            seq = str(record.seq.reverse_complement()).upper()
            match = regex.search(motif, seq, regex.BESTMATCH)
            rc_flag = 1
            rc_count += 1

        if match is not None:
            barcode = match.group(3)
            bc_index = seq.find(barcode)
            sequence = seq[0:bc_index]
            d[barcode].append(sequence)

            if barcode not in dnum:
                dnum[barcode] = 1
                ndei_index = sequence.find("CATATG")
                if ndei_index != -1:
                    kpni_index = sequence.find("TAAGGTACC")
                    if kpni_index != -1:
                        seq_between_ndei_kpni = sequence[ndei_index + 6:kpni_index]
                        dnum3[barcode] = seq_between_ndei_kpni
                        outputseqs.append(SeqRecord(id=barcode, seq=Seq(seq_between_ndei_kpni), description=""))
                        if "N" not in seq_between_ndei_kpni:
                            aaseq = Seq(seq_between_ndei_kpni).translate()
                            outputseqs_noN.append(SeqRecord(id=barcode, seq=Seq(seq_between_ndei_kpni), description=""))
                            outputseqs_noN_aa.append(SeqRecord(id=barcode, seq=aaseq))
                            doutputseqs_noN[barcode] = seq_between_ndei_kpni
                            doutputseqs_noN_aa[barcode] = str(aaseq)
                            if len(seq_between_ndei_kpni) % 3 != 0:
                                bafframecount += 1
                            if "*" not in str(aaseq):
                                goodaa += 1
                    else:
                        count_no_kpni += 1
                else:
                    count_no_ndei += 1
            else:
                dnum[barcode] += 1

    print(f"{rc_count} reads required reverse complementing")
    print(f"{len(d)} barcodes found")
    print(f"{count_no_ndei} without NdeI site")
    print(f"{count_no_kpni} without KpnI site (but have NdeI)")
    print(f"{len(outputseqs_noN)} non-ambiguous inserts")
    print(f"{goodaa} valid ORFs (no premature stop codon)")

    # --- Write outputs (same as before) ---
    with open(output_file.replace(".fasta", ".csv"), "w", newline="") as outfile:
        writer = csv.writer(outfile)
        for k, v in dnum3.items():
            writer.writerow([k, str(v)])

    with open(output_file.replace(".fasta", "_noN_aa.csv"), "w", newline="") as outfile:
        writer = csv.writer(outfile)
        for k, v in doutputseqs_noN.items():
            writer.writerow([k, str(v), str(doutputseqs_noN_aa[k])])

    with open(output_file.replace(".fasta", "bc_stats.csv"), "w", newline="") as outfile:
        writer = csv.writer(outfile)
        for k, v in dnum.items():
            writer.writerow([k, str(v)])

    with open(output_file.replace(".fasta", "bc_stats_for_starcode.tsv"), "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        for k, v in dnum.items():
            writer.writerow([k, str(v)])

    with open(output_file.replace(".fasta", "bc_list.csv"), "w", newline="") as outfile:
        writer = csv.writer(outfile)
        for k, v in d.items():
            for seq_str in v:
                writer.writerow([k, seq_str])

    print(f"Finished processing {input_file}")
