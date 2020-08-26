#' find_ctcf_in_prostate
#'
#' This function will split the res data frame into a list of data frames.
#' The list will contain main result, cpg unique, snp unique, gene unique
#' @param x dataframe of eqtl result
#' @param ctcf_dir string path to ctcf files
#' @param ctcf_metadata string path to ctcf metadata file
#' @return will output a list of data.frames
#' @import GenomicRanges
#' @import IRanges
#' @importFrom rtracklayer import.bed
#' @export
#'
find_ctcf_in_prostate <- function(x, ctcf_dir, ctcf_metadata){
    cpg_gr <- GenomicRanges::GRanges(x$seqnames,
                      IRanges(x$probePos, width = 1, names = x$probe))
    metadata <- read.table(ctcf_metadata, header = T, sep = "\t")
    df <- matrix(0, nrow = length(cpg_gr), ncol = nrow(metadata))
    for (i in 1:nrow(metadata)){
        logging::loginfo(paste("loading", metadata[i,'sample']))
        t <- rtracklayer::import(paste0(ctcf_dir,"/",metadata[i,'bedid'],".bed.gz"), format="narrowPeak")
        ol <- as.data.frame(findOverlaps(cpg_gr, t))
    df[ol$queryHits,i] <- 1
    }
    df <- as.data.frame(df)
    colnames(df) <- metadata$sample
    row.names(df) <- names(cpg_gr)
    return(df)
}
