#' make_tracks
#'
#' This function will output tracks to load in UCSC browser
#' @param x dataframe of eqtl result
#' @return will generate plot
#' @import assertthat
#' @export
#'
make_tracks <- function(x){
    eqtl_sig <- x %>% filter(sigdelta,  classSimplified!="None")
    bigint <- data.frame(
        chrom = eqtl_sig$seqnames,
        chromStart = eqtl_sig$start,
        chromEnd = eqtl_sig$end,
        name = eqtl_sig$comb,
        score = 0,
        value = abs(eqtl_sig$estimatedelta),
        exp = "sigeQTL",
        color = ifelse(eqtl_sig$siglow & !eqtl_sig$sighigh, "#EE6A50",
                       ifelse(eqtl_sig$sighigh & !eqtl_sig$siglow, "#68838B", "#4D4D4D")),
        sourceChrom = eqtl_sig$seqnames,
        sourceStart = eqtl_sig$snpPos,
        sourceEnd = eqtl_sig$snpPos +1,
        sourceName = unlist(lapply(eqtl_sig$snp, function(x) paste(rev(rev(strsplit(x,"_")[[1]])[-1]), collapse = ":"))),
        sourceStrand = ".",
        targetChrom = eqtl_sig$seqnames,
        targetStart = eqtl_sig$tss,
        targetEnd = eqtl_sig$tss + 1,
        targetName = eqtl_sig$gene,
        targetStrand = "."
    )
    write("track type=interact name=\"Significant eQTLs\" description=\"eQTLs that are different between low and high 5mc groups\" interactDirectional=true maxHeightPixels=200:100:50 visibility=full", file = "data/sig_qtl_bigInteract.txt")
    write("browser position chr16:89590592-89519592", file = "data/sig_qtl_bigInteract.txt", append = T)
    write.table(bigint, file = "data/sig_qtl_bigInteract.txt", row.names = F, quote = F, col.names = F, sep = "\t", append = T)

    t <- eqtl_sig %>% group_by(snp) %>% arrange(pdelta) %>% dplyr::slice(1)
    t <- as.data.frame(t)
    snpbed <- data.frame(t$seqnames, t$snpPos, t$snpPos + 1,
                         unlist(lapply(t$snp, function(x) paste(rev(rev(strsplit(x,"_")[[1]])[-1]), collapse = ":"))),
                         abs(t$estimatedelta)*1000)
    write("track type=bed name=\"eQTLs\" visibility=full", file = "data/sig_qtl_snps.bed")
    write("browser position chr16:89590592-89519592", file = "data/sig_qtl_snps.bed", append = T)
    write.table(snpbed, file = "data/sig_qtl_snps.bed", row.names = F, quote = F, sep = "\t", col.names = F, append = T)


    t <- eqtl_sig %>% group_by(gene) %>% arrange(pdelta) %>% dplyr::slice(1)
    t <- as.data.frame(t)
    genebed <- data.frame(t$seqnames, t$tss, t$tss + 1,
                          t$gene,
                          abs(t$estimatedelta)*1000)
    write("track type=bed name=\"eGenes\" visibility=full", file = "data/sig_qtl_genes.bed")
    write("browser position chr16:89590592-89519592", file = "data/sig_qtl_genes.bed", append = T)
    write.table(genebed,  file = "data/sig_qtl_genes.bed", row.names = F, col.names = F, append =T, quote = F, sep = "\t")

    t <- eqtl_sig %>% group_by(probe) %>% arrange(pdelta) %>% dplyr::slice(1)
    t <- as.data.frame(t)
    cpgbed <- data.frame(t$seqnames, t$probePos, t$probePos + 1,
                         t$probe,
                         abs(t$estimatedelta)*1000)
    write("track type=bed name=\"eQTLs CpGs\" visibility=full", file = "data/sig_qtl_cpgs.bed")
    write("browser position chr16:89590592-89519592", file = "data/sig_qtl_cpgs.bed", append = T)
    write.table(cpgbed, file = "data/sig_qtl_cpgs.bed", row.names = F, quote = F, sep = "\t", col.names = F, append = T)
}
