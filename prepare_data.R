#!/usr/bin/Rscript
library(reshape2)
# Little script that write to stdout the TBA matrix, starting from expression data and vcf_rider output.
# It gets the two filenames as positional arguments, the first for expression data ("expression.tsv" in the
# main text) and the second for vcf_rider output ("tba.tsv").

# Setup of input files
args <- commandArgs(TRUE)
expr_file <- args[1]
vcf_rider_file <- args[2]
expr <- read.table(expr_file, sep="\t")
tba <- read.table(vcf_rider_file, sep="\t")
colnames(tba) <- c("gene","b","e", "PWM", "sample", "allele", "tba")
tba$b <- NULL
tba$e <- NULL

# We sum the tba values of the two alleles (for the same PWM and individual) and get the log2.
recast_tba <- dcast(tba, sample+gene ~ PWM, value.var="tba", fun.aggregate=function(x) { log(sum(x), base=2) })
t_expr <- t(expr)

# We transform expression data in a format suitable to merge it with TBA using individual and gene
# ids as the join field.
m <- melt(t_expr)
m$id <- paste0(m$Var1, "_", m$Var2)
recast_tba$id <- paste0(recast_tba$sample, "_", recast_tba$gene)
merged <- merge(recast_tba, m, by="id")
colnames(merged)[ncol(merged)] <- "expression"
merged <- merged[order(c(merged$gene, merged$ind)),]
res <- merged
res$id <- NULL
res$Var1 <- NULL
res$Var2 <- NULL
res <- res[order(res$gene, res$sample),]
# Write to stdout the results.
write.table(res, file="", sep="\t", quote=FALSE, row.names = FALSE)