#' correlate_blood_prostate
#'
#' This function will calculate the correlation between blood and prostate methylation
#' @param blood_catalog "blood/catalog.N.sorted.merged.bed"
#' @param prostate_catalog "encode_correlation/catalog.N.sorted.merged.all_samples.bed"
#' @param mdata dataframe of mdata
#' @import GenomicRanges
#' @importFrom logging loginfo
#' @return a list of three elements - prostate mdata, blood mdata and correlation summary
#' @examples
#' file_blood_catalog <- "/Users/musaahmed/Google\ Drive/ctcf/blood/catalog.N.sorted.merged.bed"
#' file_prostate_catalog <- "/Users/musaahmed/Google\ Drive/ctcf/encode_correlation/catalog.N.sorted.merged.all_samples.bed"
#' file_mdata <- "/Users/musaahmed/Google\ Drive/ctcf/CPGEA_mdata.normal.txt"
#' mdata <- read.table(file_mdata, header=T, row.names = 1, sep="\t")
#' blood_prostate_cor <- correlate_blood_prostate(file_blood_catalog, file_prostate_catalog, mdata)
#' @export
#'
correlate_blood_prostate <- function(blood_catalog, prostate_catalog, mdata){
    blood <- read.table(blood_catalog, header = F, sep ="\t", stringsAsFactors = F)
    blood_loci.hg19 <- GRanges(Rle(blood$V1), IRanges(blood$V2, width = 1, names = paste0(blood$V1,":",blood$V2)),
                               cells = blood$V4, m=blood$V5, um=blood$V6)
    blood_loci.hg38 <- liftOver_wrapper(blood_loci.hg19, convertFrom = "hg19", chain = "~/Documents/DB/hg19ToHg38.over.chain")
    blood_loci.hg38 <- blood_loci.hg38[intersect(names(blood_loci.hg38), row.names(mdata))]

    t <- list()
    for (i in 1:length(blood_loci.hg38)){
        a1 <- unlist(strsplit(mcols(blood_loci.hg38)$m[i], ","))
        a2 <- unlist(strsplit(mcols(blood_loci.hg38)$um[i], ","))
        t[[i]] <- as.numeric(a1)/(as.numeric(a1)+as.numeric(a2))
        t[[i]] <- round(as.numeric(t[[i]]), 2)
        names(t[[i]]) <- unlist(strsplit(mcols(blood_loci.hg38)$cells[i], ","))
    }
    blood_samples <- names(t[[1]])
    blood <- matrix(NA, nrow = length(t), ncol = length(blood_samples))
    for (i in 1:length(blood_samples)){
        blood[,i] <- as.numeric(unlist(lapply(t, function(x) x[blood_samples[i]])))
    }
    blood_samples <- gsub("T","N",blood_samples)
    blood <- as.data.frame(blood)
    colnames(blood) <- blood_samples
    row.names(blood) <- names(blood_loci.hg38)

    prostate <- read.table(prostate_catalog, header = F, sep ="\t", stringsAsFactors = F)
    meth_loci <- paste0(prostate$V1,":",prostate$V2)
    row.names(prostate) <- meth_loci
    prostate <- prostate[intersect(row.names(mdata), row.names(prostate)),]
    meth_loci <- row.names(prostate)
    t <- list()
    for (i in 1:nrow(prostate)){
        t[[i]] <- unlist(strsplit(prostate$V5[i], ","))
        t[[i]] <- round(as.numeric(t[[i]]), 2)
        names(t[[i]]) <- unlist(strsplit(prostate$V4[i], ","))
    }
    prostate <- matrix(NA, nrow = length(t), ncol = length(blood_samples))
    for (i in 1:length(blood_samples)){
        prostate[,i] <- as.numeric(unlist(lapply(t, function(x) x[blood_samples[i]])))
    }
    prostate <- as.data.frame(prostate)
    colnames(prostate) <- blood_samples
    row.names(prostate) <- meth_loci
    mdata_matched_with_blood <- prostate[row.names(blood), ]
    mdata_matched_with_blood[is.na(mdata_matched_with_blood)] <- 0
    all(row.names(blood) == row.names(mdata_matched_with_blood))
    all(colnames(blood) == colnames(mdata_matched_with_blood))
    cors <- foreach(i=1:nrow(blood), .combine = "c") %do%{
        try(as.numeric(cor.test(as.numeric(blood[i,]), as.numeric(mdata_matched_with_blood[i,]))$estimate))
    }
    cors <- as.numeric(cors)
    names(cors) <- row.names(blood)

    cors.p <- foreach(i=1:nrow(blood), .combine = "c") %do%{
        try(cor.test(as.numeric(blood[i,]), as.numeric(mdata_matched_with_blood[i,]))$p.value)
    }
    cors.p <- as.numeric(cors.p)
    names(cors.p) <- row.names(blood)
    blood_prostate_cor <- data.frame(pcc = cors, p = cors.p)
    row.names(blood_prostate_cor) <- names(cors.p)
    blood_prostate_cor <- list(prostate = mdata_matched_with_blood,
                               blood = blood,
                               corsum = blood_prostate_cor)
    return(blood_prostate_cor)
}
