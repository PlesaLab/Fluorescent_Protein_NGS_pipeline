
#!/usr/bin/env python3
"""
Collision consensus pipeline (Python port of collision.R).

For each library (C1, C2):
  - Merge per-read barcode+sequence (parents.C{X}bc_list.csv) with STARcode clusters
    (parents.C{X}_clustered_LD1.txt) and rename barcodes to their cluster representative.
  - Extract insert between NdeI (CATATG) and KpnI (GGTACC), keeping only correct orientation.
  - Build per-barcode consensus using staged tie resolution:
      1) Majority reads
      2) Tie → identical translation → keep one
      3) Tie → longest trimmed translation
      4) Tie → most common trimmed translation
  - Outputs per library:
      * parents.C{X}_consensus_summary_v2.csv : BC, consensus_ntseq, total_reads, consensus_call
      * parents.C{X}_all_info_v2.csv         : final BC→ntseq with AA fields
"""
import pandas as pd
from Bio.Seq import Seq

def count_occurrences(haystack: str, needle: str) -> int:
    return haystack.count(needle)

def locate_first(haystack: str, needle: str):
    idx = haystack.find(needle)
    return None if idx < 0 else idx  # 0-based

def locate_last(haystack: str, needle: str):
    idx = haystack.rfind(needle)
    return None if idx < 0 else idx  # 0-based

def translate_nt(seq_nt: str) -> str:
    try:
        return str(Seq(seq_nt).translate(to_stop=False))
    except Exception:
        return ""

def trim_at_first_stop(aa: str) -> str:
    pos = aa.find("*")
    return aa[:pos] if pos != -1 else aa

def expand_cluster_rows(bcs_df: pd.DataFrame) -> pd.DataFrame:
    bcs_df = bcs_df.copy()
    bcs_df["totalBCsCollapsed"] = bcs_df["collapedBCs"].str.count(",")
    bcs_df = bcs_df.assign(collapedBCs=bcs_df["collapedBCs"].str.split(","))
    bcs_df = bcs_df.explode("collapedBCs", ignore_index=True)
    bcs_df = bcs_df.rename(columns={"bc": "tempbc", "collapedBCs": "bc"})
    return bcs_df

def extract_insert(ntseq: str):
    ndei_ct = count_occurrences(ntseq, "CATATG")
    kpni_ct = count_occurrences(ntseq, "GGTACC")
    if ndei_ct == 1 and kpni_ct == 1:
        ndei = locate_first(ntseq, "CATATG")
        kpni = locate_last(ntseq, "GGTACC")
        if ndei is not None and kpni is not None and ndei < kpni:
            return ntseq[ndei + 6:kpni]
    return None

def extract_insert_rescue(ntseq: str):
    ndei_ct = count_occurrences(ntseq, "CATATG")
    kpni_ct = count_occurrences(ntseq, "GGTACC")
    if ndei_ct > 1 and kpni_ct == 1:
        ndei = locate_last(ntseq, "CATATG")
        kpni = locate_last(ntseq, "GGTACC")
        if ndei is not None and kpni is not None and ndei < kpni:
            return ntseq[ndei + 6:kpni]
    return None

def per_lib_process(x: int):
    bcl_path = f"parents.C{x}bc_list.csv"
    bcs_path = f"parents.C{x}_clustered_LD1.txt"
    bcl = pd.read_csv(bcl_path, header=None, names=["bc","ntseq"], sep=",", dtype=str)
    bcs = pd.read_csv(bcs_path, header=None, names=["bc","numo","collapedBCs"], sep="\t", dtype=str)

    bcs_exp = expand_cluster_rows(bcs)
    bcl = (bcl.merge(bcs_exp, on="bc", how="inner")
              .rename(columns={"bc":"oldbc", "tempbc":"bc"}))
    bcl = bcl[bcl["ntseq"].str.len() > 3].copy()

    nice_rows, rescue_rows = [], []
    for _, row in bcl.iterrows():
        ins = extract_insert(row["ntseq"])
        if ins is not None:
            nice_rows.append({"bc": row["bc"], "ntseq": ins, "oldbc": row["oldbc"],
                              "numo": row["numo"], "totalBCsCollapsed": row["totalBCsCollapsed"]})
            continue
        ins2 = extract_insert_rescue(row["ntseq"])
        if ins2 is not None:
            rescue_rows.append({"bc": row["bc"], "ntseq": ins2, "oldbc": row["oldbc"],
                                "numo": row["numo"], "totalBCsCollapsed": row["totalBCsCollapsed"]})
    bcl_nice = pd.DataFrame(nice_rows)
    bcl_resc = pd.DataFrame(rescue_rows)
    bclnk = pd.concat([bcl_nice, bcl_resc], ignore_index=True)
    bclnk = bclnk[bclnk["ntseq"].str.len() > 3].copy()

    grp = (bclnk.drop(columns=["oldbc","numo","totalBCsCollapsed"])
                 .groupby(["bc","ntseq"], as_index=False).size()
                 .rename(columns={"size":"n"}))

    max_n = grp.groupby("bc", as_index=False)["n"].max().rename(columns={"n":"nmax"})
    bcl_max_rd = grp.merge(max_n, on="bc")
    bcl_max_rd = bcl_max_rd[bcl_max_rd["n"] == bcl_max_rd["nmax"]].drop(columns=["nmax"]).copy()

    num_per_bc = bcl_max_rd.groupby("bc", as_index=False).size().rename(columns={"size":"num"})
    bcl_mx_nocoll = bcl_max_rd.merge(num_per_bc.query("num==1")[["bc"]], on="bc")
    bcl_mx_coll   = bcl_max_rd.merge(num_per_bc.query("num>1")[["bc"]], on="bc")

    if not bcl_mx_coll.empty:
        trans = []
        for _, r in bcl_mx_coll.iterrows():
            aa = translate_nt(r["ntseq"])
            aatrim = trim_at_first_stop(aa)
            trans.append({"bc": r["bc"], "ntseq": r["ntseq"], "n": r["n"],
                          "aaseq": aa, "aatrim": aatrim})
        bcl_mx_coll_trans = pd.DataFrame(trans)
    else:
        bcl_mx_coll_trans = pd.DataFrame(columns=["bc","ntseq","n","aaseq","aatrim"])

    if not bcl_mx_coll_trans.empty:
        tx_cnt = (bcl_mx_coll_trans[["bc","aatrim"]].drop_duplicates()
                  .groupby("bc", as_index=False).size().rename(columns={"size":"n"}))
    else:
        tx_cnt = pd.DataFrame(columns=["bc","n"])

    same_tx_bcs = tx_cnt.query("n==1")[["bc"]]
    if not same_tx_bcs.empty:
        bcl_mx_coll_filt = (bcl_mx_coll.merge(same_tx_bcs, on="bc")
                            .groupby("bc", as_index=False).head(1))
    else:
        bcl_mx_coll_filt = pd.DataFrame(columns=bcl_mx_coll.columns)

    remaining_bcs = set(bcl_mx_coll["bc"]) - set(same_tx_bcs["bc"]) if not bcl_mx_coll.empty else set()
    if remaining_bcs:
        rem = bcl_mx_coll_trans[bcl_mx_coll_trans["bc"].isin(remaining_bcs)].copy()
        rem["aatrim_len"] = rem["aatrim"].str.len()
        rem = rem.sort_values(["bc","aatrim_len"], ascending=[True, False])
        bcl_mx_coll_rescue_by_len = (rem.groupby("bc", as_index=False).head(1)
                                        [["bc","ntseq","n"]])
    else:
        bcl_mx_coll_rescue_by_len = pd.DataFrame(columns=["bc","ntseq","n"])

    resolved_bcs = set(bcl_mx_coll_filt["bc"]) | set(bcl_mx_coll_rescue_by_len["bc"])
    if not bcl_mx_coll_trans.empty:
        unresolved_bcs = set(bcl_mx_coll_trans["bc"]) - resolved_bcs
    else:
        unresolved_bcs = set()

    if unresolved_bcs:
        rem2 = bcl_mx_coll_trans[bcl_mx_coll_trans["bc"].isin(unresolved_bcs)].copy()
        cnt = (rem2.groupby(["bc","aatrim"], as_index=False).size()
                    .rename(columns={"size":"aatrim_count"}))
        cnt = cnt.sort_values(["bc","aatrim_count"], ascending=[True, False])
        top = cnt.groupby("bc", as_index=False).head(1)[["bc","aatrim"]]
        joined = rem2.merge(top, on=["bc","aatrim"])
        bcl_mx_coll_rescue_by_common_aatrim = joined.groupby("bc", as_index=False).head(1)[["bc","ntseq","n"]]
    else:
        bcl_mx_coll_rescue_by_common_aatrim = pd.DataFrame(columns=["bc","ntseq","n"])

    cols = ["bc","ntseq","n"]
    frames = []
    for df in [bcl_mx_nocoll, bcl_mx_coll_filt, bcl_mx_coll_rescue_by_len, bcl_mx_coll_rescue_by_common_aatrim]:
        if not df.empty:
            frames.append(df[cols].copy())
    bc_final = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=cols)

    reads_per_bc = bclnk.groupby("bc", as_index=False).size().rename(columns={"size":"total_reads"})

    def label_frame(df, method):
        out = df.copy()
        out["method"] = method
        return out[["bc","ntseq","n","method"]]

    consensus_frames = []
    if not bcl_mx_nocoll.empty:
        consensus_frames.append(label_frame(bcl_mx_nocoll, "Majority reads"))
    if not bcl_mx_coll_filt.empty:
        consensus_frames.append(label_frame(bcl_mx_coll_filt, "Tie-resolved by identical translation"))
    if not bcl_mx_coll_rescue_by_len.empty:
        consensus_frames.append(label_frame(bcl_mx_coll_rescue_by_len, "Tie-resolved by longest translation"))
    if not bcl_mx_coll_rescue_by_common_aatrim.empty:
        consensus_frames.append(label_frame(bcl_mx_coll_rescue_by_common_aatrim, "Tie-resolved by most common translation"))

    if consensus_frames:
        consensus_labeled = pd.concat(consensus_frames, ignore_index=True)
        consensus_labeled = consensus_labeled.groupby("bc", as_index=False).head(1)
        consensus_summary = (consensus_labeled.merge(reads_per_bc, on="bc", how="left")
                             .assign(method=lambda d: d.apply(lambda r: "Single read" if r["total_reads"]==1 else r["method"], axis=1))
                             .rename(columns={"bc":"BC","ntseq":"consensus_ntseq","method":"consensus_call"})
                             [["BC","consensus_ntseq","total_reads","consensus_call"]]
                             .sort_values("BC"))
    else:
        consensus_summary = pd.DataFrame(columns=["BC","consensus_ntseq","total_reads","consensus_call"])

    out_sum = f"parents.C{x}_consensus_summary_v2.csv"
    consensus_summary.to_csv(out_sum, index=False)

    if not bc_final.empty:
        bc_join = bc_final.merge(bclnk, on=["bc","ntseq"], how="left").copy()
        aa_list = []
        for s in bc_join["ntseq"]:
            aa = translate_nt(s or "")
            aatrim = trim_at_first_stop(aa)
            aa_list.append((aa, aatrim))
        bc_join["aaseq"] = [t[0] for t in aa_list]
        bc_join["aatrim"] = [t[1] for t in aa_list]
        if "oldbc" in bc_join.columns:
            bc_join = bc_join.drop(columns=["oldbc"])
        bc_join = bc_join.drop_duplicates()
    else:
        bc_join = pd.DataFrame(columns=["bc","ntseq","n","aaseq","aatrim","numo","totalBCsCollapsed"])

    out_info = f"parents.C{x}_all_info_v2.csv"
    desired_order = ["bc", "ntseq", "n", "aaseq", "aatrim", "numo", "totalBCsCollapsed"]
    bc_join = bc_join.reindex(columns=desired_order)
    bc_join.to_csv(out_info, index=False)

if __name__ == "__main__":
    per_lib_process(1)
    per_lib_process(2)
