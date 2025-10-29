
# eigenMT-qvalue-significance.threshold.r - compute q-value for eigeneMT resul. and define the significance thresholds
# Usage: Rscript ~/eigenMT-qvalue-significance.threshold.r --input "Results.QTL.chr*" --output eigenMT.result_final_thresholds.txt.gz --fdr 0.05


suppressMessages(library(qvalue))
suppressMessages(library(tools))
suppressMessages(library(argparse))
suppressMessages(library(data.table))

# parse inputs
parser <- ArgumentParser(description='Annotates eigenMT output and runs qvalue')
parser$add_argument('--input', default = 'Results.QTL.chr*', help = 'Pattern to match files.')
parser$add_argument('--output', type='character', help='Output file path')
parser$add_argument('--lambda', type="double", help="", default=NULL) ## If you don't specify lambda, the qvalue function will internally choose a set of candidate lambda values and use them to estimate pi0. This is generally recommended for most users.
parser$add_argument("--fdr", type="double", help="")
args <- parser$parse_args()

#------------------------------------------------------------------------------#
# input file for eigenMT - Multiple-Testing Adjustment for eQTL Studies that Accounts for Linkage Disequilibrium between Variants
#------------------------------------------------------------------------------#
##
Files <- list.files(path = "../eigenMT/", pattern = args$input, full.names = T)
print(Files)

message ("# --------Combine results of Chr: ---------------#")
MergedData <- do.call(rbind, lapply(Files, fread))

MergedData$BF.FDR <- p.adjust(MergedData$BF, method="fdr")
MergedData$sig = ifelse(MergedData$BF.FDR <= 0.05, "TRUE", "FALSE")

print(table(MergedData$sig))
print(head(MergedData))

#------------------------------------------------------------------------------#
#  q-values (Storey) and p-value thresholds for all genes in the eigenMT results files. adapted from Francois Aguet - fastqtl/R/calculateSignificanceFastQTL.R
#------------------------------------------------------------------------------#
# remove genes w/o variants
miss <- MergedData[is.na(MergedData$BF),]
cat("  * Number of genes tested: ", nrow(MergedData), " (excluding ",
    nrow(miss), " genes w/o variants)\n", sep="")

# calculate q-values
if (is.null(args$lambda) || is.na(args$lambda)) {
  Q <- qvalue(MergedData[, 'BF'])
} else {
  cat("  * Calculating q-values with lambda = ", args$lambda, "\n", sep="")
  Q <- qvalue(MergedData[, 'BF'], lambda=args$lambda)
}

MergedData$qval <- signif(Q$qvalues, 6)
cat("  * Proportion of significant phenotypes (1-pi0): " , round((1 - Q$pi0), 2), "\n", sep="")
cat("  * eGenes @ FDR ", args$fdr, ":   ", sum(MergedData[, 'qval']<args$fdr), "\n", sep="")

# Calculate the global p-value threshold (based on sorting q-values)
ub <- min(MergedData$BF[MergedData$qval > args$fdr])
lb <- max(MergedData$BF[MergedData$qval <= args$fdr])
pthreshold <- (lb + ub) / 2  # global threshold based on p-values

cat("Global p-value threshold at FDR", args$fdr, ":", pthreshold, "\n")
# Calculate nominal p-value thresholds for each gene adjusted for Meff
MergedData$pval_nominal_threshold <- pthreshold / MergedData$TESTS
#
write.table(MergedData, gzfile(args$output), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


thresholds <- rbindlist(lapply(1:nrow(MergedData), function(g) {
  gene <- as.character(MergedData$gene[g])  # Extract gene name at row 'g'
  n.tests <- MergedData$TESTS[g]  # Extract the effective number of tests (Meff) at row 'g'
  threshold <- pthreshold / n.tests  # Calculate the nominal threshold for the gene
  return(data.frame(gene, n.tests, threshold))  # Return a data frame for each gene
}))
write.table(thresholds, "nominal_pval_thresholds.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

#------------------------------------------------------------------------------#
# Preprocess eQTL results and select SNPs based on pval_nominal_threshold “see eQTL.all.threshold.sig.txt [all significant variant-gene pairs associated with eGenes]”.
#------------------------------------------------------------------------------#
eQTL.full <- fread("../PC47_LMM.resul_with.se.N.maf.txt.gz", stringsAsFactors=FALSE)
eQTL.full$pval_nominal_threshold <- MergedData$pval_nominal_threshold[match(eQTL.full$gene_name, MergedData$gene)]

eQTL.full.sig <- eQTL.full[which(eQTL.full$p <= eQTL.full$pval_nominal_threshold), ]

write.table(eQTL.full.sig, "eQTL.all.threshold.sig.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


