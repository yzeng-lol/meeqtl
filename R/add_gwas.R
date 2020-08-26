#' add_gwas
#'
#' This function will add GWAS info
#' @param x dataframe of eqtl result
#' @param fdr fdr cutoff; default = 0.1
#' @param p p cutoff; default = 0.05
#' @return Will add two columns - sigdelta and class
#' @import assertthat
#' @importFrom plyr mapvalues
#' @examples
#' load("data/res.Rdata")
#' res <- add_siglabel_class(res)
#' @export
#'
add_gwas <- function(x, file_gwas){
    gwas <- read.table(file_gwas, header = F, sep = "\t", stringsAsFactors = F, quote = "")
    x$gwasTrait <- NA
    x[which(x$snp %in% gwas$V4), 'gwasTrait'] <- plyr::mapvalues(
        x[which(x$snp %in% gwas$V4), 'snp'], gwas$V4, gwas$V18)
    x$gwasP <- NA
    x[which(x$snp %in% gwas$V4), 'gwasP'] <- plyr::mapvalues(
        x[which(x$snp %in% gwas$V4), 'snp'], gwas$V4, gwas$V19)
    return(x)
}
