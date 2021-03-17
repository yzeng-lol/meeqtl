#' add_siglabel_class
#'
#' This function will label me_eQTLs results for their significance or activator/insulator class
#' @param x dataframe of eqtl result
#' @param fdr fdr cutoff; default = 0.1
#' @param p p cutoff; default = 0.05
#' @return Will add two columns - sigdelta and class
#' @import assertthat
#' @examples
#' load("data/res.Rdata")
#' res <- add_siglabel_class(res)
#' @export
#'
add_siglabel_class <- function(x, fdr = 0.1, p=0.05)
{
    assertthat::assert_that(is.numeric(x$pmhigh),
                            is.numeric(x$pmlow),
                            is.numeric(x$fdrlow),
                            is.numeric(x$fdrhigh),
                            is.numeric(x$pdelta),
                            is.numeric(x$fdrdelta))

    ## adding
    x$sigdelta <- ifelse(x$pdelta <= p & x$fdrdelta <= fdr, T, F)
    x$class <- ifelse(x$fdrhigh <= 0.1 & x$pmhigh <= 0.05 & x$fdrlow > 0.1, "Insulating",
                      ifelse(x$fdrhigh > 0.1 & x$pmlow <= 0.05 & x$fdrlow <= 0.1, "Activating","None"))
    x[x$fdrhigh <= 0.1 & x$pmhigh <= 0.05 & x$fdrlow <= 0.1 & x$pmlow <= 0.05,'class'] <- "Both"

    ## checking cor(CpG , CTCF)
    x[which(x$class == "Insulating" & x$encodeCorEst > 0),'class'] <- "ActivatingByPosCor"
    x[which(x$class == "Activating" & x$encodeCorEst > 0),'class'] <- "InsulatingByPosCor"
    x$class[x$sigdelta] <- paste0(x$class[x$sigdelta], "_sigdelta")
    x$classSimplified <- "None"
    x[x$class %in% c("Insulating_sigdelta", "InsulatingByPosCor_sigdelta"), 'classSimplified'] = "Insulating"
    x[x$class %in% c("Activating_sigdelta", "ActivatingByPosCor_sigdelta"), 'classSimplified'] = "Activating"

    return(x)
}
