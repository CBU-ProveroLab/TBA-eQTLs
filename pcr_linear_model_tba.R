#!/usr/bin/Rscript
# Script that starting from the TBA matrix and the covariates fits the PCR models.

# Library that implements the ddply function
.libPaths("/mnt/red/elly/bioinfotree/binary/tungsteno/local/lib/R")
library("plyr")

# set.seed initializes the generation of random numbers setting the starting value
# it must be used when the pipeline includes randomizations in order to get reproducible results
# the value chosen as seed has no meaning and therefore one can freely choose one
set.seed(42)

# load all arguments in a vector
args<-commandArgs(TRUE)
# TBA matrix with predictors (TBA on which the PCA must be done) and predicted values (expression)
data_file<-args[1]
# matrix with covariates
cov_file<-args[2]
# number of cores for the parallelization (put 1 if you do not want to parallelize)
n_cores<-args[3]
# percentage of variance that must be explained by the selected PCs
target<-args[4]
# n of permutations wanted to compute the empiric pvalue
n_perm<-args[5]

# open files
# in the TBA matrix we expect to have an 'id' column with sample identifiers, gene names and region names, separated by underscores
# in addition we expect to have an 'expression' columns in the TBA matrix and a 'sample' column in the covariates matrix
# pay attention: in the TBA matrix each block refers to a gene_region pair and all of them must be sorted according to samples identifiers

data<-read.table(data_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char="")

# all the provided covariates must be included in the null-model that is the same for all comparisons
all.cov<-read.table(cov_file, header=TRUE, sep="\t", stringsAsFactors = FALSE, comment.char = "")
all.cov <- all.cov[order(all.cov$sample),]

# check that the order of the samples is the same in each block of the TBA matrix and in covariates matrix
ids<-as.data.frame(do.call(rbind, strsplit(data$id, split="_")))
gene1 <- as.character(ids[1,2])
samples_tba <- as.character(ids[ids$V2==gene1,]$V1)
if (!all(samples_tba == all.cov$sample)) {
  stop("We need to have the same sample identifiers for individuals in the TBA and covariate matrix!")
}

all.cov$sample <- NULL
data$id <- paste0(ids$V2,"_",ids$V3)

formula<-as.formula("expression~.")

# function that is evocated recursively by ddply to make models for all genes
all_models<-function(data, all.cov, formula, target){
  
  # elimination of the column with the id
  data$id <- NULL
  
  # separation of the predictors and the predicted values
  tba<-subset(data, select=-c(expression))
  y<-subset(data, select = expression)
  
  # run the PCA on the TBA values for the selected gene
  tryCatch(
    {
      prcomp<-prcomp(tba)
      pca.tba<-prcomp$x
      # compute the eigenvalue of each PC and the corresponding percentage of variance explained
      eigs<-prcomp$sdev^2
      varexpl<-eigs/sum(eigs)*100
      # number of PCs that explain at least the wanted percentage of variance
      wantedpc <- which.max(cumsum(varexpl) >= target)
    },
    error=function(cond){
      cond <- sub("\n","", cond)
      message(paste(cond, gene, sep="\t"))
      wantedpc<<-NULL
    }
  )
  
  if(length(wantedpc)!=0){
    
    # selection of the minimum number of PCs that explain the desidered percentage of variance
    pca.tba.selected<-as.data.frame(pca.tba[,1:wantedpc])
    colnames(pca.tba.selected)<-paste('tPC', c(1:wantedpc), sep='')
    
    tryCatch(
      {
        # costruction of different data.frame and fitting of models
        # the null model must includes only the common covariates
        null_model_data<-cbind(all.cov, y)
        null_model<-lm(formula=formula, data=null_model_data)
        
        # covariates+TBA
        tba_model_data<-cbind(pca.tba.selected,all.cov,y)
        tba_model<-lm(formula=formula, data=tba_model_data)
        
        # the covariates+TBA model must be compared with the null model (F-test implemented by the ANOVA function)
        anova_tba<-anova(null_model, tba_model)
        anova_tba_pvalue<-anova_tba$`Pr(>F)`[2]
        
        # function that uses anova to compute an empirical pvalue comparing the null and the TBA model on shuffled y values
        get_perm_apval<-function(pca.data.selected, y, cov, formula){
          resample <- function(x, ...) x[sample.int(length(x), ...)]
          y_perm<-as.data.frame(resample(y))
          colnames(y_perm)<-'expression'
          # TBA+covariates model
          perm_lm_complete_data<-cbind(pca.data.selected, cov, y_perm)
          perm_lm_complete<-lm(formula=formula, data=perm_lm_complete_data)
          # "null" model, only covariates
          perm_lm_reduced_data<-cbind(cov, y_perm)
          perm_lm_reduced<-lm(formula=formula, data=perm_lm_reduced_data)
          # anova pvalues
          perm_a<-anova(perm_lm_complete, perm_lm_reduced, test="F")
          perm_apval <- perm_a$'Pr(>F)'[2]
          return(perm_apval)
        }
        
        if(n_perm > 0){
          # replicate gives us a vector with all empirical pvalues
          perms<-replicate(n_perm, get_perm_apval(pca.tba.selected, y, all.cov, formula))
          # number of pvalues < that the "real" one found in permutations
          n_smaller<-sum(perms < anova_tba_pvalue)
          # empirical pvalue
          emp_apval<-(n_smaller+1)/(as.numeric(n_perm)+1)
        }
        
        if(n_perm > 0){
          out<-c(anova_tba_pvalue, emp_apval)
          return(out)
        } else {
          out<-c(anova_tba_pvalue)
          return(out)
        }
      },
      error=function(cond){
        cond<-sub("\n","",cond)
        message(paste(cond, gene, sep="\t"))
        out<<-c(NA,NA)
      }
    )  
  } else {
    text<-"length(wantedpc)==0"
    message(paste(text, gene, sep="\t"))
    out<-c(NA,NA)
  }
  return(out)
}

# exploit ddply in order to fit a model for each gene_region pair
if (n_cores > 1) {
  n_cores <- as.numeric(n_cores)
  library(doMC)
  registerDoMC(n_cores)
  res <- ddply(.data=data, .(id), .fun=all_models, all.cov, formula, target, .parallel = TRUE)
} else {
  res <- ddply(.data=data, .(id), .fun=all_models, all.cov, formula, target)
}

if(n_perm > 0){
  colnames(res) <- c("id","pvalue", "empirical_pvalue")
} else {
  colnames(res) <- c("id","pvalue")
}

write.table(res, quote=FALSE, sep="\t", row.names=FALSE)