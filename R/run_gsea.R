#' run_gsea
#'
#' This will estimate gsea
#' @name run_gsea
#' @param x eqtl dataframe, e.g., res$gene
#' @param gmt path to gmt file
#' @param plot boolean to plot barplot, default F
#' @param pdffile path to output pdf file
#' @param siggenes optional character vector of genes to test
#' @param bggenes optional character vector of genes to use as background
#' @import clusterProfiler
#' @importFrom logging loginfo
#' @export
#'
run_gsea <- function(x,gmt, plot = F, siggenes, bggenes, pdffile, ...){
    if (missing(siggenes)){
        siggenes <- unique(x[x$sigdelta & x$classSimplified != "None", "gene"])
    }
    if (missing(bggenes)){
        bggenes <- unique(x$gene)
    }
    h <- read.gmt(gmt)
    enr <- enricher(siggenes, TERM2GENE=h, universe = bggenes, pvalueCutoff = 1, qvalueCutoff = 1, ...)
    if (plot){
        pdf(pdffile, height = 4, width = 4, useDingbats = F)
        barplot(enr, font.size = 8)
        plot(enr@result$Count, -log2(enr@result$p.adjust),
                cex.lab = 1, cex.axis = 0.8, cex.lab = 1,
             xlab = "Gene counts", frame=F,
                ylab = expression(-log[2]*" P"[adj]), ...)
        text(enr@result[enr@result$p.adjust<=0.05, 'Count'],
             -log2(enr@result[enr@result$p.adjust<=0.05, 'p.adjust']),
            labels = enr@result[enr@result$p.adjust<=0.05,'Description'], cex = 0.8)
        abline(h=-log2(0.05), lty = 2, col = "grey30")
        dev.off()

    }
    return(enr@result)
}

#' @rdname run_gsea
#' @export
#' @importFrom logging loginfo
#'
check_all_gsea <- function(x){
    gmt_path <- "~/Documents/DB/"
    gmt_list <- list.files(gmt_path, pattern = ".gmt")
    df <- data.frame(gmt = gmt_list, nsigp = 0, nsigpad = 0)
    for (i in 1:length(gmt_list)){
        logging::loginfo(paste("Reading", gmt_list[i]))
        t <- run_gsea(x, paste0(gmt_path, gmt_list[i]))
        df$nsigp[i] <- length(t$pvalue[t$pvalue<=0.05])
        df$nsigpad[i] <- length(t$pvalue[t$p.adjust<=0.1])
    }
    return(df)
}
