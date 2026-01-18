#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(vcfR)
})

args <- commandArgs(trailingOnly = TRUE)
merged_vcf <- args[1]
diag_vcf   <- args[2]
out_vcf    <- args[3]

v_merged <- read.vcfR(merged_vcf, verbose = FALSE)
v_diag   <- read.vcfR(diag_vcf,   verbose = FALSE)

# Use CHROM:POS keys for matching
key <- function(v) paste0(v@fix[, "CHROM"], ":", v@fix[, "POS"])
k_m <- key(v_merged)
k_d <- key(v_diag)

keep <- k_m %in% k_d
v_sub <- v_merged[keep, ]

write.vcf(v_sub, file = out_vcf)
cat("DONE: subset sites =", sum(keep), "\n")
