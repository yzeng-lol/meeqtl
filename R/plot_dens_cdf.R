#' plot_dens_cdf
#'
#' This function will generate densiy ot CDF plot from the eqtl result dataframe
#' @param x dataframe of me_eQTL result
#' @param class_col string column name to used as separtors
#' @param y_col string column name (must be numeric column)
#' @param nlog2 boolean if y_col should be -log2 transformed. Default F
#' @param method what plot to draw, e.g.,  density, ecdf
#' @param bw numeric density bw
#' @return will generate plot
#' @import assertthat
#' @export
#'
#'
plot_dens_cdf <- function(x, class_col = 'class', y_col, method, nlog2 = F, bw = 10000, ...)
{
    ## checking arguments
    assertthat::assert_that(!missing(y_col))
    assertthat::assert_that(!missing(method))

    x[[y_col]] <- as.numeric(x[[y_col]])
    if(nlog2) x[[y_col]] <- -log2(x[[y_col]])
    classes <- unique(x[[class_col]])

    switch(method,

           ## draw the density plot
           density = {
                        max_y <- 0
                        for (i in 1:length(classes)){
                            if (max(density(x[[y_col]][x[[class_col]]==classes[i]], bw = bw)$y) > max_y){
                                max_y <- max(density(x[[y_col]][x[[class_col]]==classes[i]], bw = bw)$y)
                            }
                        }
                        plot(density(x[[y_col]][x[[class_col]]==classes[1]], bw = bw),
                            cex.axis = 0.8,ylim=c(0,max_y), cex.lab = 1, ...)

                       ## adding extra class lines
                       for (i in 2:length(classes)){
                            lines(density(x[[y_col]][x[[class_col]]==classes[i]], bw = bw), col = i)
                       }

                      legend("bottomright", bty = "n", lty = 1, lwd = 1,
                              col=1:length(classes), legend = classes, cex = 0.8)
           },

           ## draw the ecdf plot
           ecdf = {
                    plot(ecdf(x[[y_col]][x[[class_col]]==classes[1]]), cex.axis = 0.8, cex.lab = 1, ...)
                    for (i in 2:length(classes)){
                        lines(ecdf(x[[y_col]][x[[class_col]]==classes[i]]), col = i)
                    }
                    legend("bottomright", bty = "n", lty = 1, lwd = 1, col=1:length(classes),
                            legend = classes, cex = 0.8)
                  },

         ## Stop for worong method input
         stop("Incorrect method")
    )
}

