#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape); library(phangorn); library(dplyr); library(tidyr)
  library(stringr); library(readr); library(purrr); library(ggplot2)
  library(broom); library(mgcv); library(lme4); library(lmerTest); library(FSA)
  library(cowplot)  # per pannelli figure
})

# ================== CONFIG ==================
base_dir <- "subsampling_runs_pct"         
aligndir <- file.path(base_dir, "alignments")
treedir  <- file.path(base_dir, "trees")
meta_fn  <- "expected.tsv"                 
pollen_fn<- "pollen.txt"                  
prefer_contree <- FALSE                    
collapse_uf <- TRUE                        
uf_cut <- 95

dir.create(file.path(base_dir, "plots"), showWarnings = FALSE)
outplots <- file.path(base_dir, "plots")

# ============== CARICA META =================
stopifnot(file.exists(meta_fn), file.exists(pollen_fn))
meta   <- read_tsv(meta_fn, show_col_types = FALSE) %>% distinct(sample, species)
pollen <- read_lines(pollen_fn) %>% unique()
tip2sp <- setNames(meta$species, meta$sample)

# ============== HELPERS =====================
parse_uf <- function(lbl){
  if (length(lbl)==0) return(NA_real_)
  x <- as.character(lbl)
  ifelse(grepl("/", x), as.numeric(sub(".*/","",x)), suppressWarnings(as.numeric(x)))
}

get_descendants <- function(tr, node){
  kids <- tr$edge[tr$edge[,1]==node,2]; out <- integer(0)
  for (k in kids) out <- c(out, if (k<=length(tr$tip.label)) k else get_descendants(tr,k))
  out
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
  mr <- tryCatch(getMRCA(tr, c(pol_tip, exp_tips)), error=function(e) NA)
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
  mr <- tryCatch(getMRCA(tr, c(pol_tip, exp_tips)), error=function(e) NA)
  if (is.na(mr)) return(NA_real_)
  idx <- mr - nTip; if (idx<1 || idx>length(tr$node.label)) return(NA_real_)
  parse_uf(tr$node.label[idx])
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
  v1 <- pd[cmn,cmn][upper.tri(pd[cmn,cmn])]
  v2 <- pat[cmn,cmn][upper.tri(pat[cmn,cmn])]
  if (length(v1)!=length(v2)) return(NA_real_)
  suppressWarnings(summary(lm(v1 ~ v2))$r.squared)
}

# ============== SCAN FILES ==================
phy_files <- list.files(aligndir, pattern="\\.phy$", full.names=TRUE)
stopifnot(length(phy_files)>0)

# Strategy e pct dal nome se presenti, es: sub.random.p0.50.7.phy
parse_name <- function(p){
  f <- basename(p)
  tokens <- strsplit(f, "\\.")[[1]]
  strat <- if (length(tokens)>=2) tokens[2] else NA_character_
  rep   <- suppressWarnings(as.integer(tokens[length(tokens)-1]))
  m     <- str_match(f, "\\.p([0-9]+\\.[0-9]+|[0-9]+)\\.")
  pct   <- if (!is.na(m[1,2])) m[1,2] else NA_character_
  tibble(file=f, strategy=strat, rep=rep, pct=pct, phy_path=p)
}
info <- map_dfr(phy_files, parse_name)

# Trova il tree corrispondente
tree_for <- function(phy_row){
  base <- sub("\\.phy$","", phy_row$file)
  pre  <- file.path(treedir, paste0("iq.", base))
  tf   <- paste0(pre, ".treefile")
  ct   <- paste0(pre, ".contree")
  if (prefer_contree && file.exists(ct)) return(ct)
  if (file.exists(tf)) return(tf)
  if (file.exists(ct)) return(ct)
  NA_character_
}
info$tree_path <- map_chr(seq_len(nrow(info)), ~tree_for(info[.x,]))
info <- filter(info, !is.na(tree_path) & file.exists(tree_path))

# Lunghezza allineamento e pct calcolata per strategia
len_from_phy <- function(p){
  aln <- tryCatch(read.phyDat(p, format="phylip", type="DNA"), error=function(e) NULL)
  if (is.null(aln)) return(NA_integer_)
  as.integer(ncol(as.matrix(as.DNAbin(aln))))
}
info$L <- map_int(info$phy_path, len_from_phy)
info <- info %>%
  mutate(pct_num = suppressWarnings(as.numeric(pct))) %>%
  group_by(strategy) %>%
  mutate(Lmax = max(L, na.rm=TRUE),
         pct_calc = ifelse(is.na(pct_num), round(L/Lmax, 2), pct_num),
         pct_lab  = sprintf("%.2f", pmin(pct_calc, 1))) %>%
  ungroup()

# Ordine percentuali per i plot
ord <- sprintf("%.2f", c(0.10,0.25,0.50,0.75,1.00))
info$pct_lab <- factor(info$pct_lab, levels=ord)

# ============== CORE: scoring + metriche ============
rows <- vector("list", nrow(info))
for (i in seq_len(nrow(info))) {
  phy <- info$phy_path[i]; tf <- info$tree_path[i]
  strat <- info$strategy[i]; repi <- info$rep[i]; pct_lab <- as.character(info$pct_lab[i])

  tr <- tryCatch(read.tree(tf), error=function(e) NULL)
  if (is.null(tr)) next

  TRE <- tryCatch(treeness(tr),          error=function(e) NA_real_)
  APD <- tryCatch(avg_patristic(tr),     error=function(e) NA_real_)
  SAT <- tryCatch(saturation_r2(phy,tr), error=function(e) NA_real_)

  trc <- if (isTRUE(collapse_uf)) tryCatch(collapse_by_uf(tr, uf_cut), error=function(e) NULL) else NULL

  recs <- lapply(pollen, function(pol){
    exp_sp <- tip2sp[pol]
    ok_raw <- is_sister_expected(tr,  pol, exp_sp, tip2sp)
    uf_raw <- get_node_uf(tr,      pol, exp_sp, tip2sp)
    if (!is.null(trc)) {
      ok_col <- is_sister_expected(trc, pol, exp_sp, tip2sp)
      uf_col <- get_node_uf(trc,       pol, exp_sp, tip2sp)
    } else { ok_col <- NA; uf_col <- NA }
    data.frame(strategy=strat, pct=pct_lab, rep=repi, pollen=pol, expected=exp_sp,
               sister_raw=ok_raw, ufboot_raw=uf_raw,
               sister_coll=ok_col, ufboot_coll=uf_col,
               Treeness=TRE, AvgPatristic=APD, SaturationR2=SAT)
  })
  rows[[i]] <- bind_rows(recs)
}
res <- bind_rows(rows)

# ============== SUMMARIES ============================
summ <- res %>%
  group_by(strategy, pct) %>%
  summarise(
    n=n(),
    prop_sister_raw  = mean(sister_raw,  na.rm=TRUE),
    prop_sister_coll = mean(sister_coll, na.rm=TRUE),
    med_boot_raw     = median(ufboot_raw, na.rm=TRUE),
    med_boot_coll    = median(ufboot_coll, na.rm=TRUE),
    mean_Treeness    = mean(Treeness, na.rm=TRUE),
    mean_AvgPatristic= mean(AvgPatristic, na.rm=TRUE),
    mean_SaturationR2= mean(SaturationR2, na.rm=TRUE),
    .groups="drop"
  )

write_csv(res,  file.path(base_dir, "pct_results_per_run.csv"))
write_csv(summ, file.path(base_dir, "pct_summary_by_strategy.csv"))

# ============== MODELLI / TEST =======================
# GLM per Pr(correct sister) ~ pct (numerico) su dati raw
res_num <- res %>%
  mutate(pct_num = as.numeric(as.character(pct)),
         sister_raw = as.logical(sister_raw))

glm1 <- glm(sister_raw ~ pct_num, family=binomial, data=res_num)
glm1_tidy <- broom::tidy(glm1, exponentiate=TRUE, conf.int=TRUE)
write_csv(glm1_tidy, file.path(base_dir, "glm_OR_pct.csv"))

# Soglia K95 (frazione di SNP) dal GLM
pred_fun <- function(x) predict(glm1, newdata=data.frame(pct_num=x), type="response")
rng <- range(res_num$pct_num, na.rm=TRUE)
K95 <- tryCatch(uniroot(function(x) pred_fun(x)-0.95, interval=rng)$root, error=function(e) NA_real_)
write_lines(ifelse(is.na(K95), "NA", sprintf("%.3f", K95)), file.path(base_dir, "fraction_for_95pct_correct.txt"))

# GAM (nonlinearity)
gam1 <- tryCatch(mgcv::gam(sister_raw ~ s(pct_num, k=4), family=binomial, data=res_num), error=function(e) NULL)
if (!is.null(gam1)) saveRDS(gam1, file.path(base_dir, "gam_model.rds"))

# Kruskal–Wallis + Dunn per UFboot (raw)
kw <- kruskal.test(ufboot_raw ~ as.factor(pct), data=res_num)
write_lines(paste0("KW p-value = ", signif(kw$p.value,3)), file.path(base_dir, "kw_bootstrap_raw.txt"))
dunn_df <- tryCatch(FSA::dunnTest(ufboot_raw ~ as.factor(pct), data=res_num, method="bonferroni")$res, error=function(e) NULL)
if (!is.null(dunn_df)) write_csv(dunn_df, file.path(base_dir, "dunn_bootstrap_raw.csv"))

# GLMM opzionale (se più strategie)
if (length(unique(res_num$strategy)) > 1 && "pollen" %in% names(res_num)) {
  glmm1 <- lme4::glmer(sister_raw ~ pct_num * strategy + (1|pollen), family=binomial, data=res_num)
  glmm_tidy <- broom::tidy(glmm1, conf.int=TRUE, effects="fixed", exponentiate=TRUE)
  write_csv(glmm_tidy, file.path(base_dir, "glmm_strategy_effects.csv"))
}

# ============== PLOT LINEARI DI RIEPILOGO ================================
summ$pct <- factor(summ$pct, levels=ord)

mk <- function(df,y,ttl,yl) ggplot(df, aes(pct, !!sym(y), group=strategy, color=strategy)) +
  geom_line() + geom_point() + labs(x="Fraction of total SNPs", y=yl, title=ttl) +
  theme_classic()

p1 <- mk(summ,"prop_sister_raw","Correct assignment (raw trees)","Proportion correct sister")
p1c<- mk(summ,"prop_sister_coll",sprintf("Correct assignment (collapsed; UF<%d)", uf_cut),"Proportion correct sister")
p2 <- mk(summ,"med_boot_raw","Bootstrap support (raw trees)","Median UFboot")
p2c<- mk(summ,"med_boot_coll",sprintf("Bootstrap (collapsed; UF<%d)", uf_cut),"Median UFboot")
p3 <- mk(summ,"mean_Treeness","Treeness vs SNP fraction","Treeness")
p4 <- mk(summ,"mean_AvgPatristic","Average patristic vs SNP fraction","Avg pairwise patristic")
p5 <- mk(summ,"mean_SaturationR2","Saturation vs SNP fraction","R² (p-dist ~ patristic)")

ggsave(file.path(outplots,"prop_sister_raw.pdf"), p1, width=6.8, height=4.4, dpi=300)
ggsave(file.path(outplots,"prop_sister_collapsed.pdf"), p1c, width=6.8, height=4.4, dpi=300)
ggsave(file.path(outplots,"bootstrap_raw.pdf"), p2, width=6.8, height=4.4, dpi=300)
ggsave(file.path(outplots,"bootstrap_collapsed.pdf"), p2c, width=6.8, height=4.4, dpi=300)
ggsave(file.path(outplots,"treeness.pdf"), p3, width=6.8, height=4.4, dpi=300)
ggsave(file.path(outplots,"avg_patristic.pdf"), p4, width=6.8, height=4.4, dpi=300)
ggsave(file.path(outplots,"saturation_r2.pdf"), p5, width=6.8, height=4.4, dpi=300)

# ============== VIOLIN (UFboot separati) =================
res$pct <- factor(res$pct, levels=ord)

save_violin_simple <- function(df, y, ylab, filename, digits=0){
  d <- df %>% filter(!is.na(.data[[y]]))
  if(!nrow(d)) return(invisible(NULL))
  lab <- d %>% group_by(pct) %>% summarise(med=median(.data[[y]],na.rm=TRUE), .groups="drop")
  p <- ggplot(d, aes(pct, .data[[y]], fill=pct)) +
    geom_violin(trim=FALSE, width=0.9, color="grey30", alpha=0.9) +
    stat_summary(fun=median, geom="point", size=2, color="black") +
    geom_text(data=lab, aes(x=pct, y=med, label=round(med, digits)), vjust=-1, size=3) +
    labs(x="Fraction of total SNPs", y=ylab) + theme_classic() +
    theme(legend.position="none", axis.text.x = element_text(angle=45, hjust=1))
  ggsave(file.path(outplots, filename), p, width=7.2, height=4.6, dpi=300)
}
save_violin_simple(res, "ufboot_raw", "UFboot (raw)",       "violin_ufboot_raw_by_fraction.pdf", 0)
save_violin_simple(res, "ufboot_coll", sprintf("UFboot (collapsed; UF<%d)", uf_cut),
                   "violin_ufboot_collapsed_by_fraction.pdf", 0)

# ============== VIOLIN PANEL (stessa scala Y tra metriche) ==============
# Una riga per (strategy, pct, rep): metriche per-albero + proporzione sister_raw
tree_metrics <- res %>%
  group_by(strategy, pct, rep) %>%
  summarise(
    PropSister   = mean(sister_raw, na.rm = TRUE),
    Treeness     = dplyr::first(na.omit(Treeness)),
    AvgPatristic = dplyr::first(na.omit(AvgPatristic)),
    SaturationR2 = dplyr::first(na.omit(SaturationR2)),
    .groups = "drop"
  ) %>%
  mutate(pct = factor(pct, levels = ord))
# salva per riferimento
readr::write_csv(tree_metrics, file.path(base_dir, "pct_tree_metrics_with_prop_sister.csv"))

metrics_vec <- c("PropSister","Treeness","AvgPatristic","SaturationR2")

longM <- tree_metrics %>%
  tidyr::pivot_longer(cols = all_of(metrics_vec), names_to = "metric", values_to = "value") %>%
  dplyr::filter(!is.na(value))

# limiti globali comuni (stessa scala Y per tutti i pannelli)
ymin <- min(longM$value, na.rm=TRUE)
ymax <- max(longM$value, na.rm=TRUE)

mk_violin <- function(df_metric, metric_label){
  ggplot(df_metric, aes(x=pct, y=value, fill=pct)) +
    geom_violin(trim=FALSE, width=0.9, color="grey30", alpha=0.9) +
    stat_summary(fun=median, geom="point", size=2, color="black") +
    coord_cartesian(ylim=c(ymin, ymax)) +
    labs(x="Fraction of total SNPs", y=metric_label) +
    theme_classic() +
    theme(legend.position="none", axis.text.x = element_text(angle=45, hjust=1))
}

plots_list <- lapply(metrics_vec, function(m){
  lab <- dplyr::recode(m,
    "PropSister"   = "Proportion correct sister (per tree)",
    "Treeness"     = "Treeness",
    "AvgPatristic" = "Average pairwise patristic",
    "SaturationR2" = "R² (p-dist ~ patristic)"
  )
  mk_violin(dplyr::filter(longM, metric==m), lab)
})

panel_metrics <- cowplot::plot_grid(plotlist = plots_list, ncol = 2, labels = "AUTO")
ggsave(file.path(outplots, "violin_metrics_panel_sameY.pdf"),
       panel_metrics, width=10, height=8, dpi=300)
message("[OK] Pannello violin con stessa scala Y salvato: ",
        file.path(outplots, "violin_metrics_panel_sameY.pdf"))

# ============== ALBERI SEMPLICI (senza colori/annotazioni) ==============
suppressPackageStartupMessages({
  have_ggtree <- requireNamespace("ggtree", quietly = TRUE) && requireNamespace("treeio", quietly = TRUE)
})
trees_outdir <- file.path(outplots, "trees_plain")
dir.create(trees_outdir, showWarnings = FALSE)

plot_one_tree_plain <- function(tr, base_name, suffix = "raw") {
  out_pdf <- file.path(trees_outdir, paste0(base_name, "_", suffix, ".pdf"))
  if (isTRUE(have_ggtree)) {
    p <- ggtree::ggtree(tr) +
      ggtree::geom_tiplab(size = 2.4, align = TRUE, linetype = NA, offset = 0.002) +
      ggplot2::theme_classic()
    ggplot2::ggsave(out_pdf, p, width = 7.5, height = 8.5, limitsize = FALSE)
  } else {
    grDevices::pdf(out_pdf, width = 7.5, height = 9.0)
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mar = c(2,2,2,2))
    ape::plot.phylo(tr, type = "phylogram", cex = 0.55, label.offset = 0.002)
    grDevices::dev.off()
  }
}

# loop su tutti i tree trovati
for (i in seq_len(nrow(info))) {
  tf <- info$tree_path[i]
  base_name <- sub("\\.treefile$|\\.contree$", "", basename(tf))
  tr <- tryCatch(ape::read.tree(tf), error = function(e) NULL)
  if (is.null(tr)) next
  # RAW
  plot_one_tree_plain(tr, base_name, suffix = "raw")
  # COLLAPSED (se richiesto)
  if (isTRUE(collapse_uf)) {
    trc <- tryCatch(collapse_by_uf(tr, uf_cut), error = function(e) NULL)
    if (!is.null(trc)) plot_one_tree_plain(trc, base_name, suffix = paste0("collapsed_UF<", uf_cut))
  }
}
message("[OK] Alberi semplici salvati in: ", trees_outdir)

# ============== SUPERTREE for each pct ==================
supert_outdir <- file.path(outplots, "supertrees")
dir.create(supert_outdir, showWarnings = FALSE)

build_from_collapsed <- FALSE 
.collect_trees_for_pct <- function(pct_level, collapse = FALSE) {
  tf_vec <- info %>% dplyr::filter(pct_lab == pct_level) %>% dplyr::pull(tree_path)
  trs <- purrr::map(tf_vec, ~tryCatch(ape::read.tree(.x), error = function(e) NULL)) %>% purrr::compact()
  if (!length(trs)) return(NULL)
  # opzionale: collassa UF<cut prima di costruire il supertree
  if (isTRUE(collapse)) {
    trs <- purrr::map(trs, ~tryCatch(collapse_by_uf(.x, uf_cut), error=function(e) .x))
  }
  # rimuovi alberi con <3 tip o duplicati strani
  trs <- purrr::keep(trs, ~!is.null(.x$tip.label) && length(.x$tip.label) >= 3)
  if (length(trs) < 2) return(NULL)
  class(trs) <- "multiPhylo"
  trs
}

.build_supertree_one_pct <- function(pct_level, collapse = FALSE) {
  trs <- .collect_trees_for_pct(pct_level, collapse = collapse)
  if (is.null(trs)) return(NULL)
  st <- tryCatch(
    phangorn::superTree(trs, method = "MRP", rooted = FALSE),
    error = function(e) NULL
  )
  if (is.null(st)) return(NULL)
  tag <- if (collapse) sprintf("pct_%s_collUF%d", pct_level, uf_cut) else sprintf("pct_%s_raw", pct_level)
  out_newick <- file.path(supert_outdir, paste0("supertree_", tag, ".nwk"))
  ape::write.tree(st, file = out_newick)
  st
}

pct_levels <- levels(info$pct_lab)

# --- Supertree dai RAW ---
st_list_raw <- setNames(vector("list", length(pct_levels)), pct_levels)
for (k in seq_along(pct_levels)) {
  st_list_raw[[k]] <- .build_supertree_one_pct(pct_levels[k], collapse = FALSE)
}
st_list_raw <- st_list_raw[!vapply(st_list_raw, is.null, logical(1))]

st_list_col <- NULL
if (isTRUE(build_from_collapsed)) {
  tmp <- setNames(vector("list", length(pct_levels)), pct_levels)
  for (k in seq_along(pct_levels)) {
    tmp[[k]] <- .build_supertree_one_pct(pct_levels[k], collapse = TRUE)
  }
  st_list_col <- tmp[!vapply(tmp, is.null, logical(1))]
}

+make_panel_from_st_list <- function(st_list, title_suffix, file_suffix, ncol = 3) {
  if (!length(st_list)) return(invisible(NULL))
  if (requireNamespace("ggtree", quietly = TRUE)) {
    plots_st <- lapply(names(st_list), function(lbl){
      tr <- st_list[[lbl]]
      ggtree::ggtree(tr) +
        ggtree::geom_tiplab(size = 2.2, align = TRUE, linetype = NA, offset = 0.002) +
        ggplot2::ggtitle(paste("pct =", lbl, title_suffix)) +
        ggplot2::theme_classic()
    })
    panel_st <- cowplot::plot_grid(plotlist = plots_st, ncol = ncol, labels = NULL)
    ggplot2::ggsave(file.path(supert_outdir, paste0("supertrees_panel_", file_suffix, ".pdf")),
                    panel_st, width = 12, height = 8, dpi = 300)
  } else {
    out_pdf <- file.path(supert_outdir, paste0("supertrees_panel_", file_suffix, ".pdf"))
    grDevices::pdf(out_pdf, width = 8.5, height = 10)
    for (lbl in names(st_list)) {
      tr <- st_list[[lbl]]
      par(mar = c(2,2,2,2))
      plot.new(); title(main = paste("pct =", lbl, title_suffix))
      ape::plot.phylo(tr, type = "phylogram", cex = 0.55, label.offset = 0.002)
    }
    grDevices::dev.off()
  }
}

make_panel_from_st_list(st_list_raw, title_suffix = " (raw)", file_suffix = "raw", ncol = 3)

if (isTRUE(build_from_collapsed)) {
  make_panel_from_st_list(st_list_col, title_suffix = sprintf(" (UF<%d)", uf_cut),
                          file_suffix = sprintf("collapsedUF%d", uf_cut), ncol = 3)
}

# Report
if (length(st_list_raw)) {
  message("[OK]", paste(names(st_list_raw), collapse=", "),
          " | Newick in: ", supert_outdir)
} else {
  message("[WARN] .")
}
if (isTRUE(build_from_collapsed)) {
  if (length(st_list_col)) {
    message("[OK] SuperTree COLLAPSED: ", paste(names(st_list_col), collapse=", "),
            " | Newick in: ", supert_outdir)
  } else {
    message("[WARN] .")
  }
}


