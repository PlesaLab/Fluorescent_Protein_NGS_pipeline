#!/usr/bin/env python3
"""
Parse SAM alignments and assign mutant IDs from translated sequences.

This script loads variant translations, their best nucleotide-level alignments
(SAM), and a reference protein FASTA. For each translated variant, it performs
a global pairwise amino-acid alignment versus the closest reference protein and
generates a compact mutant identifier (mutID). It outputs:
  1) A per-barcode table (BC→mutID, mutations, cigar, reads).
  2) A per-mutant table (mutID→IDalign, #BCs, mutations, seq, pct_ident, BC list, BC code).

Usage
-----
Analyze the "codon1" library:
    python parse_SAM_alignments_MP_cli.py --lib codon1

Analyze the "codon2" library:
    python parse_SAM_alignments_MP_cli.py --lib codon2

Optional paths:
    --base-path PATH           Prefix for all input relative paths (default: "")
    --analysis-dir PATH        Output directory (default: <base-path>/analysis_output/)

Notes
-----
- Library-specific file stems are selected by --lib to avoid editing the script.
- Core computation and output formats are unchanged from the original script.
"""

import argparse
import time
import csv
import hashlib
import multiprocessing
from Bio import SeqIO, pairwise2

# ------------------------------
# Utilities
# ------------------------------

def getFastaSeqs(filename):
    """Return dict {ID_without_suffix: protein_seq_wo_initial_M} from FASTA."""
    fastaseqs = {}
    with open(filename) as handle:
        for seqrec in SeqIO.parse(handle, "fasta"):
            fastaseqs[seqrec.id.split(".")[0]] = str(seqrec.seq[1:])  # drop leading 'M'
    return fastaseqs

def genMutantName(datain):
    """Pairwise-align (globalxx) and derive mutID and metrics.

    Returns (fail_flag, pct_ident, mutations, mutID, seq, ID).
    """
    ID, refseq, seq = datain
    fail_flag = 1
    pct_ident = 0.0
    mutations = -10000
    mutant_name = ""

    aln = pairwise2.align.globalxx(refseq, seq, one_alignment_only=True)
    if aln:
        fail_flag = 0
        matches = aln[0][2]
        pct_ident = matches / len(refseq)

        sizeread = len(seq.strip())
        sizeref = len(refseq.strip())

        if sizeread <= sizeref:
            mutations = int(sizeref - matches)
        else:
            # Negative sentinel for reads longer than reference
            mutations = -int(sizeref - matches + (sizeread - sizeref)) - 10000

        if (mutations < 6) and (mutations > -1):
            # Enumerate explicit mutations
            read_seq_align = aln[0][1]
            ref_seq_align = aln[0][0]
            ref_ins = 0
            name = [ID]
            for i in range(len(read_seq_align)):
                if read_seq_align[i] != ref_seq_align[i]:
                    pos = i + 1 - ref_ins + 1  # original indexing retained
                    if read_seq_align[i] == '-':         # deletion
                        name.append(f"_{ref_seq_align[i]}{pos}X")
                    elif ref_seq_align[i] == '-':         # insertion
                        name.append(f"_X{pos}{read_seq_align[i]}")
                        ref_ins += 1
                    else:                                  # substitution
                        name.append(f"_{ref_seq_align[i]}{pos}{read_seq_align[i]}")
            mutant_name = "".join(name)
        else:
            # Hash sequence for compact ID when many edits
            hex_dig = hashlib.sha256(seq.encode('utf-8')).hexdigest()
            mutant_name = f"{ID}_{hex_dig}"

    return fail_flag, pct_ident, mutations, mutant_name, seq, ID

# ------------------------------
# Main
# ------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Parse SAM alignments and generate mutant IDs with library selection."
    )
    parser.add_argument("--lib", choices=["codon1", "codon2"], required=True,
                        help="Library selection: codon1 → C1 files, codon2 → C2 files.")
    parser.add_argument("--base-path", default="", help="Base path prefix for inputs.")
    parser.add_argument("--analysis-dir", default=None,
                        help="Output directory (default: <base-path>/analysis_output/)." )
    parser.add_argument("--threads", type=int, default=10, help="Worker processes (default: 10).")
    parser.add_argument("--protein-fasta", default="FP_12.proteins",
                        help="Reference protein FASTA (default: FP_12.proteins)." )

    args = parser.parse_args()

    start_time = time.time()

    base_path = args.base_path
    analysis_dir = args.analysis_dir if args.analysis_dir is not None else base_path + 'analysis_output/'
    mapping_path_prefix = base_path

    # Library-specific stems
    stems = {
        "codon2": {
            "BC_info_file": 'parents.C2_all_info_v2',
            "trans_file":   'parents.C2_all_info_v2.csv',
            "sam_file":     'parents.C2.bc-nt',
            "mutID_BC_out": analysis_dir + 'C2seqs_mutID_all.csv',
            "mutID_info_out": analysis_dir + 'C2seqs_mutID_info_all.csv',
        },
        "codon1": {
            "BC_info_file": 'parents.C1_all_info_v2',
            "trans_file":   'parents.C1_all_info_v2.csv',
            "sam_file":     'parents.C1.bc-nt',
            "mutID_BC_out": analysis_dir + 'C1seqs_mutID_all.csv',
            "mutID_info_out": analysis_dir + 'C1seqs_mutID_info_all.csv',
        },
    }[args.lib]

    # Inputs
    BC_info_file = [stems["BC_info_file"]]
    trans_files  = [stems["trans_file"]]
    sam_files    = [stems["sam_file"]]
    protein_seq_file = [args.protein_fasta]

    # Outputs
    mutID_BC_out   = [stems["mutID_BC_out"]]
    mutID_info_out = [stems["mutID_info_out"]]

    total_num_libs = 1
    counter_lib = 0

    for _ in range(total_num_libs):
        # Accumulators
        BClist = []
        BCreads_frommapcsv = {}
        seq_to_BC_dict = {}
        seq_to_name_dict = {}
        seq_to_BCcode_dict = {}
        seq_to_info_dict = {}
        seq_to_info2_dict = [dict() for _ in range(2)]  # [mutations, pct_ident]
        seq_to_mutID_dict = {}
        mutIDcount = {}

        NUM_PROCS = args.threads

        print("load reference protein seqs and translate")
        proseqs = getFastaSeqs(protein_seq_file[counter_lib])

        # Load good barcodes after pruning
        print("Opening " + BC_info_file[counter_lib] + '.csv')
        barcode_count = 0
        for line in open(mapping_path_prefix + BC_info_file[counter_lib] + '.csv'):
            listWords = line.split(',')
            BC = listWords[0]
            if not BC or BC == '"':
                print("ERROR blank BC")
            BClist.append(BC)
            BCreads_frommapcsv[BC] = 0
            barcode_count += 1
        print(str(barcode_count) + " good barcodes loaded")

        # Load translations, truncate at first stop (*)
        print("loop over translated mapping files and load all protein sequences in-frame until the first stop codon")
        zero_length_trans_seq = 0
        maptransdict = {}
        prematurestop = []
        running_count = 0
        current_transmap = trans_files[counter_lib]

        print("Opening " + current_transmap)

        count_GT_noStop = 0
        count_GT_Stop = 0
        count_noGT_noStop = 0
        count_noGT_Stop = 0

        for line in open(mapping_path_prefix + current_transmap):
            # Columns: BC, NA, SEQ
            listWords = line.split(',')
            BC = listWords[0]
            seqt = listWords[4]
            stop_index = seqt.find('*')
            if stop_index == -1:
                seqtcut = seqt
                count_noGT_noStop += 1
                BC_code = "NoGT_noStop"
                prematurestop.append(BC)
            else:
                seqtcut = seqt[0:stop_index]
                count_noGT_Stop += 1
                BC_code = "NoGT_Stop"

            # Map seq to BCs (space-delimited) and annotate
            if seqtcut in seq_to_BC_dict:
                prior = seq_to_BC_dict[seqtcut].split(" ")
                if BC not in prior:
                    seq_to_BC_dict[seqtcut] = seq_to_BC_dict[seqtcut] + " " + BC
            else:
                seq_to_BC_dict[seqtcut] = BC
            maptransdict[BC] = seqtcut
            seq_to_BCcode_dict[seqtcut] = BC_code

            if len(seqtcut) == 0:
                zero_length_trans_seq += 1

            running_count += 1
            if running_count % 1000000 == 0:
                print(str(running_count), end='\r')

        total_BCs = len(maptransdict)
        print(str(total_BCs), " translated sequences loaded, of which ", str(len(prematurestop)), "has a premature stop codon")
        print("No GT, no stop ", str(count_noGT_noStop), " ", str(round(count_noGT_noStop / total_BCs * 100)), "%")
        print("No GT, has stop ", str(count_noGT_Stop), " ", str(round(count_noGT_Stop / total_BCs * 100)), "%")
        print("GT, no stop ", str(count_GT_noStop), " ", str(round(count_GT_noStop / total_BCs * 100)), "%")
        print("GT, has stop ", str(count_GT_Stop), " ", str(round(count_GT_Stop / total_BCs * 100)), "%")
        print("zero_length_trans_seq ", str(zero_length_trans_seq), " ", str(round(zero_length_trans_seq / total_BCs * 100)), "%")
        print(str(len(maptransdict)), " BCs loaded")
        print(str(len(seq_to_BCcode_dict)), " unique seqs loaded")

        # Map BC → best ID from SAM
        print("for each BC find the corresponding ID from alignment SAM files")
        IDalign_count = 0
        running_count = 0
        sum_nonzeropos = 0

        current_sam = sam_files[counter_lib]
        print("opening " + current_sam)
        for line in open(mapping_path_prefix + current_sam + '.map.sam'):
            if line[0] == '@':
                continue
            listWords = line.split("\t")
            BC = listWords[0].strip()
            stop_index = BC.find('_part_2')
            if stop_index != -1:
                BC = BC[0:stop_index]

            IDalign = listWords[2].split(".")[0]
            pos = int(listWords[3])
            cigar = listWords[5]
            if pos > 1:
                sum_nonzeropos += 1

            if BC in maptransdict and IDalign != '*':
                seq = maptransdict[BC]
                seq_to_info_dict[seq] = IDalign + " " + cigar
                if IDalign:
                    IDalign_count += 1

            running_count += 1
            if running_count % 1000000 == 0:
                print(str(running_count), end='\r')

        print(str(len(seq_to_info_dict)), " aligned IDs matched to seqs")
        print('Processed ', str(running_count), ' BCs')
        print('', str(IDalign_count), ' with an ID')
        print(str(sum_nonzeropos), ' BCs with pos > 1')

        # Collect good seqs
        print("for each good BC generate a list of good seq")
        no_translation_counter = 0
        no_align_counter = 0
        for BC in BClist:
            if BC in maptransdict:
                seq = maptransdict[BC]
                if seq in seq_to_info_dict:
                    seq_to_name_dict[seq] = ''
                else:
                    no_align_counter += 1
            else:
                no_translation_counter += 1
        print(str(no_translation_counter) + " good BCs with no translation.")
        print(str(no_align_counter) + " good BCs with no alignment.")
        print(str(len(seq_to_name_dict)) + " good seqs added to dictionary list.")

        # Build inputs for alignment workers
        print("determine new mutant (mutID) name for each translated seq")
        func_args_pack = []
        for seq in seq_to_name_dict:
            ID = seq_to_info_dict[seq].split()[0]
            refseq = proseqs[ID]
            func_args_pack.append([ID, refseq, seq])
        print("finished packing data, starting worker pool")

        # Parallel alignment → mutID
        pool = multiprocessing.Pool(NUM_PROCS)
        results = pool.imap_unordered(genMutantName, func_args_pack, chunksize=1000)

        pairwise_align_fail_count = 0
        running_count = 0
        for fail_flag, pct_ident, mutations, mutant_name, seq, ID in results:
            if fail_flag == 0:
                seq_to_info2_dict[1][seq] = pct_ident
                seq_to_info2_dict[0][seq] = mutations
                seq_to_mutID_dict[seq] = mutant_name
            else:
                pairwise_align_fail_count += 1
                print('ID: ' + ID)
                print('refe: ' + proseqs[ID])
                print('read: ' + seq)

            running_count += 1
            if running_count % 1000 == 0:
                print(f"{running_count} of {len(seq_to_name_dict)} {100*running_count/len(seq_to_name_dict):.1f}", end='\r')

        pool.close()
        pool.join()
        print('workers done, pool closed')
        print(str(pairwise_align_fail_count) + ' pairwise alignments failed')
        print("processed ", str(running_count), " total seqs")
        print("", str(len(seq_to_mutID_dict)), " total unique mutIDs")

        # Output 1: BC→mutID table
        csvfile = open(mutID_BC_out[counter_lib], 'w')
        fieldnames = ['BC', 'mutID', 'mutations', 'cigar', 'reads']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        BC_out_count = 0
        for BC in BClist:
            if BC in maptransdict:
                seq = maptransdict[BC]
                if seq and (seq in seq_to_info_dict) and (seq in seq_to_mutID_dict):
                    if seq not in seq_to_BC_dict:
                        # Rare: sequence not in seq_to_BC_dict; continue but note in stdout.
                        print(seq + " not in seq_to_BC_dict")
                    writer.writerow({
                        'BC': BC,
                        'mutID': seq_to_mutID_dict[seq],
                        'mutations': str(seq_to_info2_dict[0][seq]),
                        'cigar': seq_to_info_dict[seq].split()[1],
                        'reads': BCreads_frommapcsv[BC],
                    })
                    BC_out_count += 1
                    mutID = seq_to_mutID_dict[seq]
                    mutIDcount[mutID] = mutIDcount.get(mutID, 0) + 1

        print(str(BC_out_count) + " BCs saved in " + mutID_BC_out[counter_lib])
        csvfile.close()

        # Output 2: per-mutant info
        csvfile = open(mutID_info_out[counter_lib], 'w')
        fieldnames = ['mutID', 'IDalign', 'numBCs', 'mutations', 'seq', 'pct_ident', 'BCs', 'BCcode']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        mutID_count = 0
        for seq in seq_to_info_dict:
            if seq and (seq in seq_to_mutID_dict) and (seq in seq_to_BC_dict):
                mutID = seq_to_mutID_dict[seq]
                writer.writerow({
                    'mutID': mutID,
                    'IDalign': seq_to_info_dict[seq].split()[0],
                    'numBCs': mutIDcount[mutID],
                    'mutations': str(seq_to_info2_dict[0][seq]),
                    'seq': seq,
                    'pct_ident': str(seq_to_info2_dict[1][seq]),
                    'BCs': seq_to_BC_dict[seq],
                    'BCcode': seq_to_BCcode_dict[seq],
                })
                mutID_count += 1
        csvfile.close()
        print(str(mutID_count) + " mutIDs saved in " + mutID_info_out[counter_lib])

        counter_lib += 1

    print("--- %s seconds ---" % (time.time() - start_time))
