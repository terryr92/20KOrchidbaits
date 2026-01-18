#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(vcfR)
})

args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]   # Multi-sample VCF for reference species (filtered, no het)
meta_csv <- args[2]   # sample,species mapping
out_vcf  <- args[3]   # Output VCF (species-specific)

vcf <- read.vcfR(vcf_file, verbose = FALSE)
gt  <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

meta <- read.csv(meta_csv, stringsAsFactors = FALSE)
stopifnot(all(colnames(gt) %in% meta$sample))
meta <- meta[match(colnames(gt), meta$sample), ]

# Helpers
is_missing <- function(x) x %in% c("./.", ".", NA)
is_hom_ref <- function(x) x %in% c("0/0", "0|0")
is_hom_alt <- function(x) x %in% c("1/1", "1|1")
is_hom     <- function(x) is_hom_ref(x) || is_hom_alt(x)

species <- unique(meta$species)

keep_any <- rep(FALSE, nrow(gt))
tag_sp   <- rep(NA_character_, nrow(gt))

# Strict taxon-aware definition:
# For a given focal species:
# - all individuals in focal species are fixed homozygous (0/0 or 1/1) and identical
# - all individuals in all other species are fixed for the opposite allele
# - no missing genotypes anywhere at that site
for (sp in species) {
  idx_sp <- which(meta$species == sp)
  idx_ot <- which(meta$species != sp)

  g_sp <- gt[, idx_sp, drop = FALSE]
  g_ot <- gt[, idx_ot, drop = FALSE]

  no_missing <- apply(gt, 1, function(r) all(!is_missing(r)))

  sp_fixed <- apply(g_sp, 1, function(r) {
    all(is_hom(r)) && length(unique(r)) == 1
  })

  # Determine the focal allele (0/0 or 1/1) for each site (only if fixed)
  sp_allele <- apply(g_sp, 1, function(r) unique(r)[1])

  ot_opposite <- mapply(function(focal, others) {
    if (!is_hom(focal)) return(FALSE)
    if (is_hom_ref(focal)) {
      return(all(others %in% c("1/1","1|1")))
    } else {
      return(all(others %in% c("0/0","0|0")))
    }
  }, focal = sp_allele, others = split(g_ot, row(g_ot)))

  keep_sp <- no_missing & sp_fixed & ot_opposite

  new_hits <- which(keep_sp & !keep_any)
  keep_any[new_hits] <- TRUE
  tag_sp[new_hits] <- sp
}

vcf_out <- vcf[keep_any, ]

# Annotate focal species in INFO (SPECIES=<name>)
info <- vcf_out@fix[, "INFO"]
info <- ifelse(info == ".", paste0("SPECIES=", tag_sp[keep_any]),
               paste0(info, ";SPECIES=", tag_sp[keep_any]))
vcf_out@fix[, "INFO"] <- info

write.vcf(vcf_out, file = out_vcf)
cat("DONE: species-specific sites =", sum(keep_any), "\n")
