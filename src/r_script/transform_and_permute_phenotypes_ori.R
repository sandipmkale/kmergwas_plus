library(MASS) #ginv
library(mvnpermute) #mvnpermute
library(matrixcalc) #is.positive.semi.definite

args_for_path <- commandArgs(trailingOnly = FALSE)
args <- commandArgs(trailingOnly = TRUE)
script.name <- sub("--file=", "", args_for_path[grep("--file=", args_for_path)])
path_emma <- file.path(dirname(script.name), "emma.R")
source(path_emma)

fn_phenotypes <- args[1]
fn_kinship <- args[2]
n_permute <- as.numeric(args[3])
fn_out_phenotypes <- args[4]
fn_out_trans_phenotypes <- args[5]
f_log <- file(args[6])

phenotypes <- read.csv(fn_phenotypes,sep='\t',header=TRUE)
av_pheno <- mean(phenotypes$phenotype_value)
phenotypes$phenotype_value <- phenotypes$phenotype_value-av_pheno
n_acc = nrow(phenotypes)

K <- as.matrix(read.csv(fn_kinship, sep='\t',header=FALSE))

if(!is.positive.semi.definite(K)) {
  writeLines('Kinship matrix is not positive semi-definite'); quit()
}

null <- emma.REMLE(phenotypes$phenotype_value,as.matrix(x = rep(1, n_acc),dim = c(n_acc,1)),K)
herit <- null$vg/(null$vg+null$ve)
COV_MATRIX <- null$vg*K+null$ve*diag(dim(K)[1])
CM_inv <- ginv(COV_MATRIX)

writeLines(c(paste('EMMA_n_permutation','=',n_permute), 
             paste('EMMA_n_accessions','=',n_acc), 
             paste('EMMA_vg','=',null$vg), 
             paste('EMMA_ve','=',null$ve), 
             paste('EMMA_herit','=',herit)),f_log)
close(f_log)

if(n_permute > 0) {
  permute_phenotype <- mvnpermute(phenotypes$phenotype_value, rep(1, n_acc), COV_MATRIX, nr=n_permute)
  for(i in 1:n_permute) { phenotypes[paste("P",i, sep = '')] <- permute_phenotype[,i] }
}

trans_phenotypes <- phenotypes
for(i in 2:(n_permute+2)) { trans_phenotypes[,i] <- CM_inv %*% phenotypes[,i] }

write.table(x=phenotypes, file=fn_out_phenotypes, eol = "\n",sep = "\t", quote=F, row.names = FALSE)
write.table(x=trans_phenotypes, file=fn_out_trans_phenotypes, eol = "\n",sep = "\t", quote=F, row.names = FALSE)
