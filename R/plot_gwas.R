#' plot_gwas
#'
#' This function will plot gwas data
#' @param x dataframe of eqtl result
#' @param pdfname path to pdf file to save
#' @return will generate plot
#' @import assertthat
#' @import RColorBrewer
#' @export
#'
plot_gwas <- function(x, pdfname){
    assertthat::assert_that(x %has_name% "gwasTrait")
    x <- x[order(x$fdrdelta),]
    pdf(pdfname, height = 4, width = 4, useDingbats = F)
    par(ps = 12, mgp=c(1.5,0.5,0))
    plot(seq(1,nrow(x)), -log2(x$fdrdelta), type = "p", frame = F, xlab = "eQTLs ranked by effect size", ylab = expression(-log[2]*" P"[adj]),
         cex.lab = 1, cex.axis = 0.8, cex = 0.2)
    label_idx <- grep("cancer", x$gwasTrait, ignore.case = T)
    label_idx <- intersect(label_idx, grep("response", x$gwasTrait, ignore.case = T, invert = T))
    label_idx <- intersect(label_idx, grep("gene", x$gwasTrait, ignore.case = T, invert = T))
    label_idx <- intersect(label_idx, grep("levels", x$gwasTrait, ignore.case = T, invert = T))
    label_idx <- intersect(label_idx, grep("treatment", x$gwasTrait, ignore.case = T, invert = T))
    t <- as.factor(x[label_idx,'gwasTrait'])
    jColors <- data.frame(cancer = levels(t),
                        color = I(brewer.pal(nlevels(t), name = 'Set3')))
    points(label_idx, -log2(x[label_idx,'fdrdelta']), cex = 0.4,
           col = jColors$color[match(x[label_idx,'gwasTrait'], jColors$cancer)])
    abline(h=-log2(0.1), lty = 2)
    legend(x = 'topright',
           legend = as.character(jColors$cancer),
           col = jColors$color, bty = 'n', pch = 16, xjust = 1, cex = 0.5)
    dev.off()
}
