# TBA-eQTLs
Supplementary material for "Mapping of eQTL using transcription factors affinity

```
$ vcf_rider -v wgs.vcf -p pwms.tsv -b regulatory.bed -r reference.fa -f bg.tsv -a regulatory_snps.tsv > tba.tsv
[W::vcf_parse] contig '1' is not defined in the header. (Quick workaround: index the file with tabix.)
$ ./prepare_data.R expression.tsv tba.tsv > tba_matrix.tsv
$. /pcr_linear_model_tba.R tba_matrix.tsv covariates.tsv 1 95 10 >
sample_models.tsv
```
