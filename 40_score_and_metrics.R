#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape); library(phangorn); library(dplyr); library(tidyr)
  library(stringr); library(readr); library(purrr); library(ggplot2)
})

# ========= CONFIG (tutto in questa cartella) =========
outdir    <- "."               # dove salvare CSV/PNG
treedir   <- "."               # .treefile qui
aligndir  <- "."               # .phy qui
meta_fn   <- "expected.tsv"    # sample<TAB>species
pollen_fn <- "pollen.txt"      # lista pollini

# opzional
do_collapse <- TRUE
uf_cut      <- 95

# ========= META =========
meta   <- read_tsv(meta_fn, show_col_types = FALSE)
pollen <- read_lines(pollen_fn) %>% unique()
tip2sp <- setNames(meta$species, meta$sample)

# ========= HELPER =========
parse_uf <- function(lbl) {
  if (length(lbl)==0) return(NA_real_)
  x <- as.character(lbl)
  ifelse(grepl("/", x), as.numeric(sub(".*/","",x)), suppressWarnings(as.numeric(x)))
}
get_descendants <- function(tr, node){
  kids <- tr$edge[tr$edge[,1]==node,2]; des <- integer(0)
  for (k in kids) des <- c(des, if (k<=length(tr$tip.label)) k else get_descendants(tr,k))
  des
}
collapse_by_uf <- function(tr, cutoff){
  if (is.null(tr$node.label)) return(tr)
  uf <- parse_uf(tr$node.label); nTip <- length(tr$tip.label)
  bad <- which(!is.na(uf) & uf < cutoff) + nTip
  if (!length(bad)) return(tr)
  idx <- which(tr$edge[,2] %in% bad)
  if (!is.null(tr$edge.length)) tr$edge.length[idx] <- 0
  di2multi(tr, tol=0)
}
is_sister_expected <- function(tr, pol_tip, exp_species, map_sp){
  if (!(pol_tip %in% tr$tip.label)) return(NA)
  exp_tips <- intersect(names(map_sp)[map_sp==exp_species], tr$tip.label)
  if (!length(exp_tips)) return(NA)
  mr <- tryCatch(getMRCA(tr, c(pol_tip,exp_tips)), error=function(e) NA)
  if (is.na(mr)) return(NA)
  parent_edge <- which(tr$edge[,2]==mr); if (!length(parent_edge)) return(NA)
  parent <- tr$edge[parent_edge,1]
  childs <- tr$edge[tr$edge[,1]==parent,2]; if (length(childs)<2) return(NA)
  sister <- childs[childs!=mr]; sis_lab <- tr$tip.label[get_descendants(tr,sister)]
  all(sis_lab %in% exp_tips)
}
get_node_uf <- function(tr, pol_tip, exp_species, map_sp){
  nTip <- length(tr$tip.label); if (is.null(tr$node.label)) return(NA_real_)
  exp_tips <- intersect(names(map_sp)[map_sp==exp_species], tr$tip.label)
  if (!length(exp_tips)) return(NA_real_)
  mr <- tryCatch(getMRCA(tr, c(pol_tip,exp_tips)), error=function(e) NA)
  if (is.na(mr)) return(NA_real_)
  idx <- mr - nTip; if (idx<1 || idx>length(tr$node.label)) return(NA_real_)
  parse_uf(tr$node.label[idx])
}
pvs_from_phylip <- function(phyfile){
  if (!file.exists(phyfile)) return(NA_real_)
  aln <- read.phyDat(phyfile, format="phylip", type="DNA")
  X <- as.matrix(as.DNAbin(aln)); L <- ncol(X); if (L==0) return(NA_real_)
  var_cols <- 0
  for (j in seq_len(L)){
    col <- toupper(as.character(X[,j])); col <- col[col!="N"]
    if (length(unique(col))>1) var_cols <- var_cols + 1
  }
  var_cols / L
}
treeness <- function(tr){
  if (is.null(tr$edge.length)) return(NA_real_)
  nTip <- length(tr$tip.label); internal <- tr$edge[,2] > nTip
  sum(tr$edge.length[internal], na.rm=TRUE) / sum(tr$edge.length, na.rm=TRUE)
}
avg_patristic <- function(tr){
  C <- cophenetic(tr); mean(C[upper.tri(C)], na.rm=TRUE)
}
saturation_r2 <- function(phyfile,tr){
  if (!file.exists(phyfile)) return(NA_real_)
  aln <- read.phyDat(phyfile, format="phylip", type="DNA")
  X <- as.DNAbin(aln)
  pd  <- dist.dna(X, model="raw", pairwise.deletion=TRUE, as.matrix=TRUE)
  pat <- cophenetic(tr)
  cmn <- intersect(rownames(pd), rownames(pat))
  if (length(cmn)<3) return(NA_real_)
  v1 <- pd[cmn,cmn][upper.tri(pd[cmn,cmn])]; v2 <- pat[cmn,cmn][upper.tri(pat[cmn,cmn])]
  if (length(v1)!=length(v2)) return(NA_real_)
  suppressWarnings(summary(lm(v1 ~ v2))$r.squared)
}

# ========= SCAN FILES (no “pct” required) =========

phy_files <- list.files(aligndir, pattern="^sub\\..*\\.phy$", full.names=TRUE)
stopifnot(length(phy_files)>0)

# funzione per estrarre strategy e rep dal nome; pct se presente, altrimenti NA
parse_name <- function(path){
  f <- basename(path)
  # esempi: sub.random.p0.50.7.phy | sub.bias.250.7.phy | sub.random.K500.12.phy
  # strategy = token dopo 'sub.'
  tokens <- unlist(strsplit(f, "\\."))
  # expected: sub <strategy> <maybe-pXXXX or K or number> <rep> phy  (ma può variare)
  # strategy:
  strat <- if (length(tokens)>=2) tokens[2] else NA_character_
  # rep: ultimo numero prima di "phy"
  rep <- suppressWarnings(as.integer(tokens[length(tokens)-1]))
  # pct se c'è pattern pX o pX.Y
  m <- str_match(f, "\\.p([0-9]+\\.[0-9]+|[0-9]+)\\.")
  pct <- if (!is.na(m[1,2])) m[1,2] else NA_character_
  tibble(file=f, strategy=strat, rep=rep, pct=pct)
}

info <- map_dfr(phy_files, ~parse_name(.x)) %>%
  mutate(phy_path = phy_files)


find_tree_for_phy <- function(phy_row){
  f <- phy_row$file
  # prova stesso core cambiando "sub." con "iq."
  cand <- sub("^sub\\.", "iq.", f)
  path1 <- file.path(treedir, sub("\\.phy$", ".treefile", cand))
  if (file.exists(path1)) return(path1)
  # pattern più permissivo: cerca iq.<strategy>.*.<rep>.treefile
  patt <- paste0("^iq\\.", phy_row$strategy, ".*\\.", phy_row$rep, "\\.treefile$")
  hits <- list.files(treedir, pattern=patt, full.names=TRUE)
  if (length(hits)>=1) return(hits[1])
  return(NA_character_)
}

info$tree_path <- map_chr(seq_len(nrow(info)), ~find_tree_for_phy(info[.x,]))
info <- filter(info, !is.na(tree_path) & file.exists(tree_path))

len_from_phy <- function(p){
  aln <- tryCatch(read.phyDat(p, format="phylip", type="DNA"), error=function(e) NULL)
  if (is.null(aln)) return(NA_integer_)
  as.integer(ncol(as.matrix(as.DNAbin(aln))))
}
info$L <- map_int(info$phy_path, len_from_phy)

info <- info %>%
  group_by(strategy) %>%
  mutate(Lmax = max(L, na.rm=TRUE),
         pct_calc = ifelse(!is.na(L) & Lmax>0, round(L/Lmax, 2), NA_real_),
         pct_final = ifelse(!is.na(pct),
                            as.numeric(pct),  # es. "0.50" → 0.50
                            pct_calc         # altrimenti usa percentuale calcolata
                            )
         ) %>%
  ungroup()

ord <- c(0.10,0.25,0.50,0.75,1.00)
closest_level <- function(x, levels=ord){
  if (is.na(x)) return(NA_character_)
  lev <- levels[which.min(abs(x - levels))]
  sprintf("%.2f", lev)
}
info$pct_lab <- vapply(info$pct_final, closest_level, character(1))

# ========= CORE: scoring + metriche =========
rows <- list()
for (i in seq_len(nrow(info))) {
  phy <- info$phy_path[i]; tf <- info$tree_path[i]
  strat <- info$strategy[i]; repi <- info$rep[i]; pct_lab <- info$pct_lab[i]

  if (!file.exists(phy) || !file.exists(tf)) next
  tr <- tryCatch(read.tree(tf), error=function(e) NULL)
  if (is.null(tr)) next

  # metriche (raw)
  PVS <- tryCatch(pvs_from_phylip(phy), error=function(e) NA_real_)
  TRE <- tryCatch(treeness(tr),          error=function(e) NA_real_)
  APD <- tryCatch(avg_patristic(tr),     error=function(e) NA_real_)
  SAT <- tryCatch(saturation_r2(phy,tr), error=function(e) NA_real_)


  trc <- if (isTRUE(do_collapse)) tryCatch(collapse_by_uf(tr, uf_cut), error=function(e) NULL) else NULL

  for (pol in pollen) {
    exp_sp <- tip2sp[pol]
    ok_raw <- is_sister_expected(tr,  pol, exp_sp, tip2sp)
    uf_raw <- get_node_uf(tr,      pol, exp_sp, tip2sp)
    if (!is.null(trc)) {
      ok_col <- is_sister_expected(trc, pol, exp_sp, tip2sp)
      uf_col <- get_node_uf(trc,       pol, exp_sp, tip2sp)
    } else { ok_col <- NA; uf_col <- NA }

    rows[[length(rows)+1]] <- data.frame(
      strategy=strat, pct=pct_lab, rep=repi, pollen=pol, expected=exp_sp,
      sister_raw=ok_raw, ufboot_raw=uf_raw,
      sister_coll=ok_col, ufboot_coll=uf_col,
      PVS=PVS, Treeness=TRE, AvgPatristic=APD, SaturationR2=SAT,
      stringsAsFactors=FALSE
    )
  }
}

res <- bind_rows(rows)

# ========= SUMMARIES & PLOTS =========
summ <- res %>%
  group_by(strategy, pct) %>%
  summarise(
    n=n(),
    prop_sister_raw  = mean(sister_raw, na.rm=TRUE),
    prop_sister_coll = mean(sister_coll, na.rm=TRUE),
    mean_boot_raw    = mean(ufboot_raw, na.rm=TRUE),
    med_boot_raw     = median(ufboot_raw, na.rm=TRUE),
    mean_boot_coll   = mean(ufboot_coll, na.rm=TRUE),
    med_boot_coll    = median(ufboot_coll, na.rm=TRUE),
    mean_PVS         = mean(PVS, na.rm=TRUE),
    mean_Treeness    = mean(Treeness, na.rm=TRUE),
    mean_AvgPatristic= mean(AvgPatristic, na.rm=TRUE),
    mean_SaturationR2= mean(SaturationR2, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(pct=factor(pct, levels=sprintf("%.2f", c(0.05,0.10,0.25,0.50,0.75,1.00))))

write_csv(res,  file.path(outdir, "pct_results_per_run.csv"))
write_csv(summ, file.path(outdir, "pct_summary_by_strategy.csv"))

# Plot (raw + collapsed se attivo)
mk_line <- function(df, y, ttl, ylab){
  ggplot(df, aes(pct, !!sym(y), group=strategy, color=strategy)) +
    geom_line() + geom_point() +
    labs(x="Fraction of total SNPs", y=ylab, title=ttl) +
    theme_classic()
}
p1  <- mk_line(summ, "prop_sister_raw", "Correct assignment (raw trees)", "Proportion correct sister")
p2  <- mk_line(summ, "med_boot_raw",    "Bootstrap support (raw trees)",  "Median UFboot")
p3  <- mk_line(summ, "mean_Treeness",   "Treeness vs SNP fraction",       "Treeness")
p4  <- mk_line(summ, "mean_AvgPatristic","Avg patristic vs SNP fraction", "Avg pairwise patristic")
p5  <- mk_line(summ, "mean_SaturationR2","Saturation vs SNP fraction",    "R² (p-dist ~ patristic)")
ggsave(file.path(outdir,"pct_prop_sister_raw.png"), p1, width=6.5, height=4.2, dpi=300)
ggsave(file.path(outdir,"pct_median_bootstrap_raw.png"), p2, width=6.5, height=4.2, dpi=300)
ggsave(file.path(outdir,"pct_treeness.png"), p3, width=6.5, height=4.2, dpi=300)
ggsave(file.path(outdir,"pct_avg_patristic.png"), p4, width=6.5, height=4.2, dpi=300)
ggsave(file.path(outdir,"pct_saturation_r2.png"), p5, width=6.5, height=4.2, dpi=300)

if (isTRUE(do_collapse)) {
  p1c <- mk_line(summ, "prop_sister_coll",
                 sprintf("Correct assignment (collapsed; UF<%d)", uf_cut),
                 "Proportion correct sister")
  p2c <- mk_line(summ, "med_boot_coll",
                 sprintf("Bootstrap (collapsed; UF<%d)", uf_cut),
                 "Median UFboot")
  ggsave(file.path(outdir,"pct_prop_sister_collapsed.png"), p1c, width=6.5, height=4.2, dpi=300)
  ggsave(file.path(outdir,"pct_median_bootstrap_collapsed.png"), p2c, width=6.5, height=4.2, dpi=300)
}

message("[OK] Finish !!!!!")
