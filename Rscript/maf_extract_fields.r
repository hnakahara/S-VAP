# Rscripts maf_filter.r inputfile1 inputfile2 outputfile1 outputfile2

inputfile1 = commandArgs(trailingOnly=TRUE)[1]
inputfile2 = commandArgs(trailingOnly=TRUE)[2]
outputfile1 = commandArgs(trailingOnly=TRUE)[3]
outputfile2 = commandArgs(trailingOnly=TRUE)[4]

fields <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Tumor_Sample_Barcode", "HGVSc", "HGVSp_Short", "t_depth", "t_alt_count", "tumor_vaf", "CLIN_SIG", "ExAC_AF", "gnomAD_AF", "SIFT", "PolyPhen", "MUTATION_EFFECT")

matrix1 <- read.csv(inputfile1, header=T, sep="\t")
matrix_filter1 <- matrix1[, colnames(matrix1) %in% fields]
write.table(matrix_filter1, outputfile1, na="", row.names=F, quote=F, sep="\t")

matrix2 <- read.csv(inputfile2, header=T, sep="\t")
matrix_filter2 <- matrix2[, colnames(matrix2) %in% fields]
write.table(matrix_filter2, outputfile2, na="", row.names=F, quote=F, sep="\t")
