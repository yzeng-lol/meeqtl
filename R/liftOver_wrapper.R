#' liftOver_wrapper
#'
#' This is a wrapper for liftover function
#' @param x GRanges object
#' @param chain string link to chain file
#' @param convertFrom string to identify assembly version converting from
#' @param changeName boolean, if T (default), the lifted over ranges will have chr:pos naming; if F, will retain original names as in x
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom rtracklayer liftOver import.chain
#' @export
#'
liftOver_wrapper <- function(x,convertFrom="hg19", chain, changeName = T){
    #x=Granges object
    genome(x) <- convertFrom
    GenomeInfoDb::seqlevelsStyle(x) = "UCSC"
    ch = import.chain(chain)
    res = liftOver(x, ch)
    if (changeName){
        res <- unlist(res)
        names(res) <- paste0(as.character(seqnames(res)),":",as.character(start(res)))
        return(res)
    } else {
        return(unlist(res, use.names = T))
    }
}
