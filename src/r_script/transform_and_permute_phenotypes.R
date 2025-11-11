##!/usr/bin/env Rscript
# transform_and_permute_phenotypes.R
# Robust Rscript that loads emma.R from the same folder, permutes phenotypes,
# and applies GRAMMAR-Gamma transformation.

suppressMessages({
   library(MASS)        # ginv
     library(mvnpermute)  # mvnpermute
     library(matrixcalc)  # is.positive.semi.definite
})

# -----------------------------
# Parse arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
   stop(paste(
	          "Usage:",
		      "Rscript transform_and_permute_phenotypes.R",
		      "<phenotypes.tsv> <kinship.tsv> <n_perm> <out_pheno.tsv> <out_trans.tsv> <logfile>",
		          sep = "\n  "
		        ))
}

fn_phenotypes            <- args[1]
fn_kinship               <- args[2]
n_permute                <- as.numeric(args[3])
fn_out_phenotypes        <- args[4]
fn_out_trans_phenotypes  <- args[5]
logfile                  <- args[6]

# -----------------------------
# Robustly locate and source emma.R
# -----------------------------
args_full   <- commandArgs(trailingOnly = FALSE)
script_path <- NULL

# Try sys.frames first
if (!is.null(sys.frames()) &&
        length(sys.frames()) >= 1 &&
	    !is.null(sys.frames()[[1]]$ofile)) {
   script_path <- sys.frames()[[1]]$ofile
} else {
   # Try from commandArgs
   file_arg <- grep("^--file=", args_full, value = TRUE)
  if (length(file_arg) > 0) {
       script_path <- sub("^--file=", "", file_arg[1])
    }
}

# Fallback: current working directory
if (is.null(script_path) || script_path == "" || is.na(script_path)) {
   script_path <- getwd()
}

script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
emma_path  <- file.path(script_dir, "emma.R")

if (!file.exists(emma_path)) {
   stop(paste("emma.R not found at", emma_path))
}
source(emma_path)

# -----------------------------
# Load data
# -----------------------------
phenotypes <- read.csv(fn_phenotypes, sep = "\t", header = TRUE, check.names = FALSE)
if (!"phenotype_value" %in% colnames(phenotypes)) {
   stop("phenotypes file must contain a column named 'phenotype_value'")
}
av_pheno <- mean(phenotypes$phenotype_value)
phenotypes$phenotype_value <- phenotypes$phenotype_value - av_pheno
n_acc <- nrow(phenotypes)

K <- as.matrix(read.csv(fn_kinship, sep = "\t", header = FALSE))
if (nrow(K) != ncol(K)) {
   stop("Kinship matrix must be square")
}
if (!is.positive.semi.definite(K)) {
   stop("Kinship matrix is not positive semi-definite")
}

# -----------------------------
# EMMA variance components
# -----------------------------
# Intercept-only design matrix
X_int <- matrix(1, nrow = n_acc, ncol = 1)
null   <- emma.REMLE(phenotypes$phenotype_value, X_int, K)
herit  <- null$vg / (null$vg + null$ve)

COV_MATRIX <- null$vg * K + null$ve * diag(nrow(K))
CM_inv     <- ginv(COV_MATRIX)

# -----------------------------
# Logging
# -----------------------------
con <- file(logfile, open = "w")
writeLines(c(
	       paste("Loaded emma.R from:", emma_path),
	         paste("EMMA_n_permutation =", n_permute),
	         paste("EMMA_n_accessions =", n_acc),
		   paste("EMMA_vg =", null$vg),
		   paste("EMMA_ve =", null$ve),
		     paste("EMMA_herit =", herit)
		   ), con)
close(con)

# -----------------------------
# Permutations (if requested)
# -----------------------------
if (n_permute > 0) {
permute_phenotype <- mvnpermute(phenotypes$phenotype_value, rep(1, n_acc), COV_MATRIX, nr=n_permute)
  for (i in seq_len(n_permute)) {
       phenotypes[[paste0("P", i)]] <- permute_phenotype[, i]
    }
}

# -----------------------------
# GRAMMAR-Gamma transformation
# -----------------------------
trans_phenotypes <- phenotypes
# Column 1 is IDs; phenotype is column 2; permutations (if any) follow
last_col <- 2 + n_permute
for (i in 2:last_col) {
   trans_phenotypes[, i] <- as.vector(CM_inv %*% phenotypes[, i])
}

# -----------------------------
# Write outputs
# -----------------------------
write.table(phenotypes,
	                file = fn_out_phenotypes,
			            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(trans_phenotypes,
	                file = fn_out_trans_phenotypes,
			            sep = "\t", quote = FALSE, row.names = FALSE)

