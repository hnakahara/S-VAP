# Rscripts maf_filter.r inputdir1 outputfile1 outputfile2 outputfile3

inputdir1 = commandArgs(trailingOnly=TRUE)[1]
outputfile1 = commandArgs(trailingOnly=TRUE)[2]
outputfile2 = commandArgs(trailingOnly=TRUE)[3]
outputfile3 = commandArgs(trailingOnly=TRUE)[4]

filepath <- list.files(inputdir1, recursive=T, full.names=T, pattern="_variant_class.maf")
files <- lapply(filepath, function(i){read.csv(i, header=T, sep="\t")})
mergefile <- do.call("rbind", files)
write.table(mergefile, outputfile1, na="", row.names=F, quote=F, sep="\t")

filepath2 <- list.files(inputdir1, recursive=T, full.names=T, pattern="_variant_class_vus_prior.maf")
files2 <- lapply(filepath2, function(i){read.csv(i, header=T, sep="\t")})
mergefile2 <- do.call("rbind", files2)
write.table(mergefile2, outputfile2, na="", row.names=F, quote=F, sep="\t")

filepath3 <- list.files(inputdir1, recursive=T, full.names=T, pattern="_merged_raw.maf")
files3 <- lapply(filepath3, function(i){read.csv(i, header=T, sep="\t")})
mergefile3 <- do.call("rbind", files3)
write.table(mergefile3, outputfile3, na="", row.names=F, quote=F, sep="\t")
