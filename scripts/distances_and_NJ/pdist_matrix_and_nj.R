#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(vcfR)
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
out_tree <- args[2]
out_csv  <- args[3]

vcf <- read.vcfR(vcf_file, verbose = FALSE)
gt  <- extract.gt(vcf, "GT", as.numeric = FALSE)

# Encode genotypes: 0/0 -> 0, 1/1 -> 1, missing -> NA
encode <- function(x) {
  if (x %in% c("0/0", "0|0")) return(0)
  if (x %in% c("1/1", "1|1")) return(1)
  return(NA_real_)
}

X <- apply(gt, c(1,2), encode)
X <- t(X)  # samples x sites

# p-distance on comparable sites (genotype-aware for homozygous calls)
dist_pair <- function(a, b) {
  ok <- which(!is.na(a) & !is.na(b))
  if (length(ok) == 0) return(NA_real_)
  mean(a[ok] != b[ok])
}

n <- nrow(X)
D <- matrix(0, n, n, dimnames = list(rownames(X), rownames(X)))
for (i in seq_len(n)) {
  for (j in i:n) {
    d <- dist_pair(X[i, ], X[j, ])
    D[i, j] <- d
    D[j, i] <- d
  }
}

write.csv(D, out_csv, quote = FALSE)
tree <- nj(as.dist(D))
write.tree(tree, file = out_tree)

cat("DONE: NJ tree ->", out_tree, "\n")
