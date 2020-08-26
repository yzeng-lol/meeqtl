
# this script will analysis each case -------------------------------------

probeid <- "chr8:23415796"
probegr <- GRanges(Rle("chr8"), IRanges(23415796, width = 1, names = probeid))
liftOver_wrapper(probegr, convertFrom = "Hg38", chain = "~/Documents/DB/hg38ToHg19.over.chain")

plot_blood_cor(probeid, outfile = paste0("plots/blood_",probeid,".pdf"))
