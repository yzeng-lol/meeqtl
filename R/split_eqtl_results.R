#' split_eqtl_results
#'
#' This function will split the res data frame into a list of data frames.
#' The list will contain main result, cpg unique, snp unique, gene unique
#' @param x dataframe of eqtl result
#' @return will output a list of data.frames
#' @import tidyverse
#' @importFrom dplyr slice
#' @export
#'
split_eqtl_results <- function(x){
    t <- list(ori = x,
              cpg = as.data.frame(x %>% group_by(probe) %>% arrange(pdelta) %>% dplyr::slice(1) ),
              snp = as.data.frame(x %>% group_by(snp) %>% arrange(pdelta) %>% dplyr::slice(1) ),
              gene = as.data.frame(x %>% group_by(gene) %>% arrange(pdelta) %>% dplyr::slice(1) )
    )
    row.names(t$cpg) <- t$cpg$probe
    row.names(t$snp) <- t$snp$snp
    row.names(t$gene) <- t$gene$gene
    return(t)
}
