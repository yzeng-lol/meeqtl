#' get_histone_signals
#'
#' This function will extract histone levels in probe
#' @name get_histone_signals
#' @rdname get_histone_signals
#' @title Extract histone ChIP-seq signals for probe regions
#' @param x dataframe of eqtl result
#' @return will output a list of data.frames
#' @import GenomicRanges
#' @import IRanges
#' @import rtracklayer
#' @importFrom ChIPpeakAnno reCenterPeaks
#' @importFrom logging loginfo
#' @import assertthat
#' @export
#'
get_histone_signals <- function(x){
    cpg_gr <- GenomicRanges::GRanges(x$seqnames,
                      IRanges(x$probePos, width = 1, names = x$probe))
    cpg_gr <- ChIPpeakAnno::reCenterPeaks(cpg_gr, 500)
    cpg_gr.hg19 <- liftOver_wrapper(cpg_gr, convertFrom = "hg38", chain = "~/Documents/DB/hg38ToHg19.over.chain", changeName = F)
    cpg_gr.hg19 <- cpg_gr.hg19[intersect(names(cpg_gr.hg19), names(cpg_gr))]

    metadata <- read.table("/Volumes/Seagate_Backup_Plus/chipatlas/selected.tab", header = F, sep = "\t", stringsAsFactors = F)
    metadata <- rbind(metadata, c("GSM3058128_22rv1-h3k4me3_acttgatg_l006_r1.fastq.bam_spmr", "hg19","Histone","H3K4me3","Prosatte","22Rv1"))
    metadata <- metadata[metadata$V1!="SRX539657",]

    get_bw <- function(what_histone, what_cell){
        bwfile <- paste0("/Volumes/Seagate_Backup_Plus/chipatlas/bigwigs/",metadata[metadata$V4==what_histone & metadata$V6==what_cell,'V1'],".bw")
        bw <- BigWigFile(bwfile[1])
        return(bw)
    }

    histone_score <- list()
    for (histone in unique(metadata$V4)){
        assertthat::is.string(histone)
        cells <- unique(metadata[metadata$V4==histone,'V6'])
        assertthat::is.string(cells)
        histone_score[[histone]] <- matrix(NA, nrow = length(cpg_gr.hg19), ncol = length(cells))
        for (i in 1:length(cells)){
            cell <- cells[i]
            logging::loginfo(paste("loading", histone, "of", cell))
            bw <- get_bw(histone,cell)
            out.gr <- unlist(try(summary(bw, which=cpg_gr.hg19, type="max")))
            histone_score[[histone]][,i] <- as.numeric(mcols(out.gr)[['score']])
        }
        histone_score[[histone]] <- as.data.frame(histone_score[[histone]])
        colnames(histone_score[[histone]]) <- cells
        row.names(histone_score[[histone]]) <- names(cpg_gr.hg19)
    }
    return(histone_score)
}

#' @rdname get_histone_signals
#' @param histone_list output from get_histone_signals
#' @param histone which histone
#' @param cell which cell
#' @export

plot_histone_density <- function(x,histone_list, histone, cell, ...){
    t <- x[row.names(histone_list[[1]]),]
    assertthat::are_equal(row.names(t), row.names(histone_list[[1]]))
    assertthat::assert_that(histone_list %has_name% histone)
    assertthat::assert_that(histone_list[[histone]] %has_name% cell)
    plot(density(histone_list[[histone]][t$classSimplified=="Insulating",cell], na.rm = T), cex.axis = 0.8, cex.lab = 1,
         cex.main = 1, main=paste(cell, histone), frame = F, col = 2, ...)
    lines(density(histone_list[[histone]][t$classSimplified=="Activating",cell], na.rm = T))
    legend("topright", bty = "n", cex = 0.8, lty = 1, col=c(2,1), legend = c("Insulating","Activating"))
}
