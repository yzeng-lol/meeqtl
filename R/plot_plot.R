#' plot_plot
#'
#' This function will plot from the eqtl result dataframe
#' @param x dataframe of eqtl result
#' @param class_col string column name to used as separtors
#' @param y_col string column name (must be numeric column)
#' @param nlog2 boolean if y_col should be -log2 transformed. Default F
#' @param method what plot to draw, e.g., boxplot, density, ecdf
#' @param bw numeric density bw
#' @return will generate plot
#' @import assertthat
#' @examples
#' load("data/res.Rdata")
#' res <- add_siglabel_class(res)
#' @export
#'
plot_plot <- function(x, class_col = 'class', y_col, method, nlog2 = F, bw = 10000, ...){
    assertthat::assert_that(!missing(y_col))
    assertthat::assert_that(!missing(method))
    x[[y_col]] <- as.numeric(x[[y_col]])
    if(nlog2) x[[y_col]] <- -log2(x[[y_col]])
    classes <- unique(x[[class_col]])
    switch(method,
           boxplot = {
               boxplot(x[[y_col]] ~ x[[class_col]], cex.axis = 0.8, cex.lab = 1, ...)
           },
           density = {
               max_y <- 0
               for (i in 1:length(classes)){
                   if (max(density(x[[y_col]][x[[class_col]]==classes[i]], bw = bw)$y) > max_y){
                       max_y <- max(density(x[[y_col]][x[[class_col]]==classes[i]], bw = bw)$y)
                   }
               }
               plot(density(x[[y_col]][x[[class_col]]==classes[1]], bw = bw), cex.axis = 0.8,
                    ylim=c(0,max_y), cex.lab = 1, ...)
               for (i in 2:length(classes)){
                   lines(density(x[[y_col]][x[[class_col]]==classes[i]], bw = bw), col = i)
               }
               legend("bottomright", bty = "n", lty = 1, lwd = 1,
                      col=1:length(classes),
                      legend = classes, cex = 0.8)
           },
           ecdf = {
                plot(ecdf(x[[y_col]][x[[class_col]]==classes[1]]), cex.axis = 0.8, cex.lab = 1, ...)
               for (i in 2:length(classes)){
                   lines(ecdf(x[[y_col]][x[[class_col]]==classes[i]]), col = i)
               }
               legend("bottomright", bty = "n", lty = 1, lwd = 1,
                      col=1:length(classes),
                      legend = classes, cex = 0.8)
           },
           stop("Incorrect method")
    )
}

#' @rdname plot_plot
#' @param probeid string probeid
#' @export
#'
plot_blood_cor <- function(probeid, outfile){
    pdf(outfile, height = 4, width = 4, useDingbats = F)
    par(mgp=c(1.5,0.5,0), ps = 12)
    plot(as.numeric(blood_prostate_cor$prostate[probeid,]),
         as.numeric(blood_prostate_cor$blood[probeid,]), frame = F,
         cex.lab = 1, cex.axis = 0.8, cex.main = 1,
         main = probeid, xlab = "Prostate", ylab = "Blood",
         pch = 16, ylim=c(0.6,1))
    abline(lm(as.numeric(blood_prostate_cor$prostate[probeid,]) ~
                  as.numeric(blood_prostate_cor$blood[probeid,])), col = 2)
    t <- cor.test(as.numeric(blood_prostate_cor$prostate[probeid,]),
                  as.numeric(blood_prostate_cor$blood[probeid,]))
    if(t$estimate<0){
        legend("bottomleft", bty = "n", cex = 0.8,
           legend = paste("PCC =", round(t$estimate, 2),"\nP-value =",
                          round(t$p.value, 2)))
    } else {
        legend("topright", bty = "n", cex = 0.8,
               legend = paste("PCC =", round(t$estimate, 2),"\nP-value =",
                              round(t$p.value, 2)))
    }
    dev.off()
}
