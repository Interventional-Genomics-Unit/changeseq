# Import libraries
library("DESeq2")
library("dplyr")
library("bcbioRNASeq")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")

# DEBUG Input
rawCounts_in <- read.csv("reports/toy_raw_counts.csv", row.names = 1)

# Take variables from Snakemake
rawCounts_in <- read.csv(snakemake@input[["raw_counts"]], row.names = 1)
output_vst <- snakemake@output[["vst_counts"]]
output_tmm <- snakemake@output[["tmm_counts"]]

# Round up float numbers in raw counts created by Salmon
rawCounts_in <- rawCounts_in %>%
		mutate_if(is.numeric, round)
rawCounts_in <- rawCounts_in[rowSums(rawCounts_in != 0) > 0,
							   colSums(rawCounts_in != 0) > 0]
rowlabs <- rownames(rawCounts_in)
rawCounts_int <- as.data.frame(lapply(rawCounts_in,as.integer))
# Convert df into matrix
rawCounts_mat <- as.matrix(rawCounts_int)

# Normalize results
print("TMM Normalization")
rawCounts_tmm <- tmm(rawCounts_mat)
tmm_dataframe <- as.data.frame(rawCounts_tmm)
rownames(tmm_dataframe) <- rowlabs

# print("RLOG Normalization")
# rawCounts_rlog <- rlog(rawCounts_mat, fitType = "local", blind=FALSE)
# rlog_dataframe <- as.data.frame(rawCounts_rlog)
# rownames(rlog_dataframe) <- rowlabs

print("Variance Stabilizing transformation")
rawCounts_vst <- vst(rawCounts_mat, fitType = "parametric", blind=FALSE)
vst_dataframe <- as.data.frame(rawCounts_vst)
rownames(vst_dataframe) <- rowlabs

# Export output files
write.csv(vst_dataframe, file=output_vst, row.names = T)
write.csv(tmm_dataframe, file=output_tmm, row.names = T)
# write.csv(rlog_dataframe, file=output_vst, row.names = T)
