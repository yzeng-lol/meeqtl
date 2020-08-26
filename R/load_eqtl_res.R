#' load_eqtl_res
#'
#' This function will load eqtl result and retail one case per atac-ctcf-gene combination
#' @param x eqtl result file
#' @param cordat "encode_correlation/methylation_correlated_with_ctcf.Padj_filtered.bed"
#' @param atac "CPGEA_gdata.ATAC.txt"
#' @import tidyverse
#' @import logging
#' @importFrom plyr mapvalues
#' @keywords TCGA,metadata,boxplot
#' @examples
#' res<-load_eqtl_res("/Users/musaahmed/Google\ Drive/ctcf/out_all.noNA.normalized.txt", "/Users/musaahmed/Google\ Drive/ctcf/encode_correlation/methylation_correlated_with_ctcf.Padj_filtered.bed","/Users/musaahmed/Google\ Drive/ctcf/CPGEA_gdata.ATAC.txt")
#' @export
#'
load_eqtl_res <- function(x, cordat, atac){
    logging::loginfo("Reading file")
    eqtl_res <- read.table(x, header = T, sep = "\t", stringsAsFactors = F)
    eqtl_res <- eqtl_res %>% filter(!is.na(as.numeric(pmhigh)) & !is.na(as.numeric(pmlow)))
    logging::loginfo("Reading encode cor data")
    cpg_pca_ctcf <- read.table(cordat, header = T, sep = "\t", stringsAsFactors = F)
    logging::loginfo("Reading SNP ATAC file")
    atac <- read.table(atac, header = F, sep = "\t", stringsAsFactors = F)
    eqtl_res$ctcf <- NA
    eqtl_res$ctcf <- plyr::mapvalues(eqtl_res$probe, cpg_pca_ctcf$cpg, cpg_pca_ctcf$ctcf, warn_missing = F)
    eqtl_res$atac <- NA
    eqtl_res$atac <- plyr::mapvalues(eqtl_res$snp, atac$V12, atac$V4, warn_missing = F)
    eqtl_res$comb <- paste0(eqtl_res$atac,":",eqtl_res$ctcf,":",eqtl_res$gene)
    test <- eqtl_res %>% group_by(probe) %>% mutate(fdrhigh = p.adjust(pmhigh, "fdr"),
                                                    fdrlow = p.adjust(pmlow, "fdr"),
                                                    fdrdelta = p.adjust(pdelta, "fdr"))
    eqtl_res <- as.data.frame(test)
    t <- eqtl_res %>% group_by(comb) %>% arrange(pdelta) %>% slice(1)
    eqtl_res <- as.data.frame(t)
    row.names(eqtl_res) <- eqtl_res$comb
    eqtl_res$encodeCorEst <- NA
    eqtl_res$encodeCorPadj <- NA
    logging::loginfo("Adding ENCODE correlation data")
    for (i in 1:nrow(eqtl_res)){
        eqtl_res$encodeCorEst[i] <- cpg_pca_ctcf[cpg_pca_ctcf$cpg==eqtl_res$probe[i] & cpg_pca_ctcf$ctcf==eqtl_res$ctcf[i],'est']
        eqtl_res$encodeCorPadj[i] <- cpg_pca_ctcf[cpg_pca_ctcf$cpg==eqtl_res$probe[i] & cpg_pca_ctcf$ctcf==eqtl_res$ctcf[i],'padj']
    }
    return(eqtl_res)
}
