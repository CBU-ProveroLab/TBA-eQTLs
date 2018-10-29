#!/usr/bin/Rscript

# Little script that write to stdout the TBA matrix, starting from expression data, vcf_rider output
# and a file mapping regulatory regions to genes
# It gets the three filenames as positional arguments, the first for expression data ("expression.tsv" in the
# main text), the second for vcf_rider output ("tba.tsv") and the third for mapping file ("mapping.tsv")

library(reshape2)

# Setup of input files
args <- commandArgs(TRUE)
expr_file <- args[1]
vcf_rider_file <- args[2]
mapping_file <- args[3]

expr <- read.table(expr_file, sep="\t")
tba <- read.table(vcf_rider_file, sep="\t")
colnames(tba) <- c("region","b","e", "PWM", "sample", "allele", "tba")
tba$b <- NULL
tba$e <- NULL
mapping <- read.table(mapping_file, sep="\t")
colnames(mapping) <- c("region","gene")

# We sum the tba values of the two alleles (for the same PWM and individual) and get the log2.
recast_tba <- dcast(tba, sample+region ~ PWM, value.var="tba", fun.aggregate=function(x) { log(sum(x), base=2) })

# We transform expression data in a format suitable to merge it with TBA using individual and gene
# ids as the join field.
t_expr <- t(expr)
m <- melt(t_expr)
m$id <- paste0(m$Var1, "_", m$Var2)

# We associate target genes to the corresponding regulatory regions
# and then transform the tba values in a forma suitable to merge with gene expression values
recast_tba_withgenes <- merge(recast_tba, mapping, by="region")
recast_tba_withgenes$id <- paste0(recast_tba_withgenes$sample, "_", recast_tba_withgenes$gene)

merged <- merge(recast_tba_withgenes, m, by="id")
colnames(merged)[ncol(merged)] <- "expression"
merged <- merged[order(merged$gene, merged$region, merged$sample),]

res <- merged
res$id <- paste0(res$id, "_", res$region)
res$region <- NULL
res$sample <- NULL
res$gene <-NULL
res$Var1 <- NULL
res$Var2 <- NULL

# Write to stdout the results.
write.table(res, file="", sep="\t", quote=FALSE, row.names = FALSE)