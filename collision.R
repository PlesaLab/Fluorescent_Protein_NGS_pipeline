# collision.R
# -----------
# For each library (C1, C2), merge raw barcode+sequence reads with STARcode clusters,
# trim inserts between NdeI (CATATG) and KpnI (GGTACC), and call a consensus per barcode.
# Tie resolution (in order):
#   1) Majority reads
#   2) If tie → identical translation → keep one
#   3) If still tied → keep longest trimmed translation
#   4) If still tied → keep most common trimmed translation
#
# Outputs (per library):
#   - parents.C?_consensus_summary_v2.csv : BC, consensus_ntseq, total_reads, consensus_call
#   - parents.C?_all_info_v2.csv         : final BC→ntseq with AA translation fields
#
library(stringi)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(seqinr)

# Robust "script directory" detection
get_script_dir <- function() {
  # In RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  # When run via Rscript
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 1) {
    p <- sub("^--file=", "", file_arg)
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  # Fallback: current working directory
  getwd()
}

setwd(get_script_dir())


# loop over two files (C1, C2)
for (x in 1:2) {
  
  # Load raw per-read barcode+sequence and STARcode clusters
  bcl <- read.table(file = paste('parents.C',as.character(x),'bc_list.csv',sep=""), sep =',', header = FALSE)
  bcs <- read.table(file = paste('parents.C',as.character(x),'_clustered_LD1.txt',sep=""), sep ='	', header = FALSE)

  # Name columns
  colnames(bcl) <- c("bc","ntseq")
  colnames(bcs) <- c("bc","numo","collapedBCs")

  # Expand collapsed clusters to one row per (barcode member)
  bcs <- bcs %>%
    mutate(totalBCsCollapsed=str_count(collapedBCs, ",")) %>%
    separate_rows(collapedBCs) %>%
    dplyr::rename(tempbc=bc, bc=collapedBCs)

  # Merge reads with cluster representative, rename barcodes to representative
  bcl <- bcl %>%
    inner_join(bcs,by="bc") %>%
    dplyr::rename(oldbc=bc, bc=tempbc) %>%
    rowwise() %>%
    dplyr::filter(nchar(as.character(ntseq)) > 3)

  rm(bcs)

  # Extract DNA between NdeI and KpnI (one of each, correct order)
  bcl_nice_ndei_kpni <- bcl %>%
    mutate(ndei=str_count(ntseq,"CATATG"),
         kpni=str_count(ntseq,"GGTACC")) %>%
    dplyr::filter(ndei==1 & kpni==1) %>%
    mutate(ndei=stri_locate_first_fixed(ntseq,"CATATG")[,1],
           kpni=stri_locate_last_fixed(ntseq,"GGTACC")[,1]) %>%
    dplyr::filter(ndei<kpni) %>%
    mutate(ntseq2=substr(ntseq,ndei+6,kpni-1)) %>%
    dplyr::select(-ndei, -kpni, -ntseq) %>%
    dplyr::rename(ntseq=ntseq2)

  # Rescue: multiple NdeI but single KpnI → take last NdeI and last KpnI
  bcl_ndei_rescues <- bcl %>%
    mutate(ndei=str_count(ntseq,"CATATG"),
           kpni=str_count(ntseq,"GGTACC")) %>%
    dplyr::filter(ndei>1 & kpni==1) %>%
    mutate(ndei=stri_locate_last_fixed(ntseq,"CATATG")[,1],
           kpni=stri_locate_last_fixed(ntseq,"GGTACC")[,1]) %>%
    dplyr::filter(ndei<kpni) %>%
    mutate(ntseq2=substr(ntseq,ndei+6,kpni-1)) %>%
    dplyr::select(-ndei, -kpni, -ntseq) %>%
    dplyr::rename(ntseq=ntseq2)

  # Keep inserts only (between NdeI and KpnI)
  bclnk <- dplyr::bind_rows(bcl_nice_ndei_kpni, bcl_ndei_rescues) %>%
    dplyr::filter(!is.na(ntseq) & nchar(as.character(ntseq)) > 3)

  rm(bcl_nice_ndei_kpni,bcl_ndei_rescues)

  # Counts per (bc, ntseq) to identify majority calls
  bclnk_uniq_rd_counts_per_bc <- bclnk %>% 
    dplyr::select(-oldbc, -numo, -totalBCsCollapsed) %>% 
    dplyr::group_by_all() %>% 
    dplyr::count() %>% 
    ungroup()

  # Max count per barcode
  bclnk_max_for_bc <- bclnk_uniq_rd_counts_per_bc %>%
    group_by(bc) %>%
    summarise(ms=max(n)) %>%
    dplyr::rename(n=ms) %>%
    ungroup()

  # Majority per barcode (highest-count ntseq)
  bcl_max_rd <- bclnk_uniq_rd_counts_per_bc %>%
    semi_join(bclnk_max_for_bc,by=c("bc","n"))

  rm(bclnk_uniq_rd_counts_per_bc, bclnk_max_for_bc)

  # Split: no tie (n==1) vs tie (n>1) for the max count
  bcl_max_rd_num <- bcl_max_rd %>% 
    group_by(bc) %>%
    dplyr::count() %>%
    ungroup()

  bcl_mx_nocoll <- semi_join(bcl_max_rd, 
                             bcl_max_rd_num %>% dplyr::filter(n==1) %>% dplyr::select(-n),
                             by="bc")

  bcl_mx_coll <- semi_join(bcl_max_rd, 
                           bcl_max_rd_num %>% dplyr::filter(n>1) %>% dplyr::select(-n),
                           by="bc")

  rm(bcl_max_rd_num)

  # Translate tied ntseqs; trim to first stop
  bcl_mx_coll_trans <- bcl_mx_coll %>%
    rowwise() %>%
    mutate(aaseq=seqinr::c2s(seqinr::getTrans(seqinr::s2c(ntseq)))) %>%
    mutate(aatrim = if_else(str_detect(aaseq, "\\*"),
                            substr(aaseq, 1, stri_locate_first_fixed(aaseq, "*")[,1]-1),
                            aaseq))

  # Count distinct trimmed translations per barcode
  bcl_mx_coll_trans_cnt <- bcl_mx_coll_trans %>% 
    ungroup() %>% 
    dplyr::select(bc, aatrim) %>% 
    group_by(bc) %>% 
    distinct() %>% 
    dplyr::count()

  # Barcodes with a single translation among ties
  bcl_same_tx <- bcl_mx_coll_trans_cnt %>% dplyr::filter(n==1) %>% dplyr::select(bc)

  # Resolve ties with identical translation: keep one
  bcl_mx_coll_filt <- bcl_mx_coll %>%
    semi_join(bcl_same_tx, by="bc") %>%
    group_by(bc) %>%
    dplyr::filter(row_number()==1) %>%
    ungroup()

  # Resolve remaining ties by longest trimmed translation (if tie, keep one)
  bcl_mx_coll_rescue_by_len <- bcl_mx_coll_trans %>%
    anti_join(bcl_same_tx, by="bc") %>%
    group_by(bc) %>%
    mutate(aatrim_len = nchar(aatrim)) %>%
    dplyr::filter(aatrim_len == max(aatrim_len)) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::select(-aatrim_len) %>%
    ungroup()

  # Barcodes still unresolved after previous steps
  resolved_bcs <- bind_rows(bcl_mx_coll_filt, bcl_mx_coll_rescue_by_len) %>% pull(bc) %>% unique()
  unresolved_bcs <- bcl_mx_coll_trans %>% filter(!(bc %in% resolved_bcs)) %>% pull(bc) %>% unique()

  # Final rescue: choose most common trimmed translation (if tie, keep one)
  bcl_mx_coll_rescue_by_common_aatrim <- bcl_mx_coll_trans %>%
    dplyr::filter(bc %in% unresolved_bcs) %>%
    group_by(bc, aatrim) %>%
    mutate(aatrim_count = n()) %>%
    ungroup() %>%
    group_by(bc) %>%
    dplyr::filter(aatrim_count == max(aatrim_count)) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::select(-aatrim_count) %>%
    ungroup()

  # Combine all consensus-resolved barcodes
  bc_final <- bind_rows(
    bcl_mx_nocoll,
    bcl_mx_coll_filt,
    bcl_mx_coll_rescue_by_len,
    bcl_mx_coll_rescue_by_common_aatrim %>% dplyr::select(names(bcl_mx_coll))
  )

  # Per-BC consensus summary CSV
  reads_per_bc <- bclnk %>%
    group_by(bc) %>%
    summarise(total_reads = n(), .groups = "drop")

  consensus_majority <- bcl_mx_nocoll %>% dplyr::select(bc, ntseq, n) %>% mutate(method = "Majority reads")
  consensus_ident_tx <- bcl_mx_coll_filt %>% dplyr::select(bc, ntseq, n) %>% mutate(method = "Tie-resolved by identical translation")
  consensus_longest_tx <- bcl_mx_coll_rescue_by_len %>% dplyr::select(bc, ntseq, n) %>% mutate(method = "Tie-resolved by longest translation")
  consensus_common_tx <- bcl_mx_coll_rescue_by_common_aatrim %>% dplyr::select(bc, ntseq, n) %>% mutate(method = "Tie-resolved by most common translation")

  consensus_labeled <- dplyr::bind_rows(consensus_majority, consensus_ident_tx, consensus_longest_tx, consensus_common_tx) %>%
    dplyr::group_by(bc) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  consensus_summary <- consensus_labeled %>%
    left_join(reads_per_bc, by = "bc") %>%
    mutate(method = if_else(total_reads == 1, "Single read", method)) %>%
    transmute(
      BC = bc,
      consensus_ntseq = ntseq,
      total_reads = total_reads,
      consensus_call = method
    ) %>%
    arrange(BC)

  write.csv(
    consensus_summary,
    paste0("parents.C", as.character(x), "_consensus_summary_v2r.csv"),
    row.names = FALSE, quote = FALSE
  )

  # Attach translation fields and write final mapping
  bc_join <- bc_final %>%
    left_join(bclnk, by = c("bc","ntseq")) %>%
    rowwise() %>%
    mutate(aaseq = c2s(getTrans(s2c(ntseq)))) %>%
    mutate(aatrim = if_else(str_detect(aaseq, "\\*"),
                            substr(aaseq, 1, stri_locate_first_fixed(aaseq, "*")[,1]-1),
                            aaseq)) %>%
    dplyr::select(-oldbc) %>%
    distinct()


  rm(bclnk, bcl, bc_final)

  write.csv(bc_join, 
            paste("parents.C",as.character(x),"_all_info_v2r.csv", sep=""), 
            row.names = FALSE, quote=FALSE)

  rm(bc_join)
}
