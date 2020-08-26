#' add_blood_cols
#'
#' This function will add two columns to eqtl result
#' @param x data frame eqtl result
#' @param corlist output from correlate_blood_prostate
#' @importFrom plyr mapvalues
#' @import assertthat
#' @return x with two added columns - bloodPCC and bloodP
#' @seealso correlate_blood_prostate
#' @examples
#'load("data/res.Rdata")
#'load("data/blood_prostate_cor.Rdata")
#'res <- add_blood_cols(res, blood_prostate_cor)
#' @export
#'
add_blood_cols <- function(x, corlist){
    assertthat::has_name(corlist, "corsum")
    assertthat::not_empty(corlist$corsum)
    x$bloodPCC <- NA
    x$bloodPCC <- plyr::mapvalues(x$probe, row.names(corlist$corsum), corlist$corsum$pcc, warn_missing = F)
    x$bloodP <- NA
    x$bloodP <- plyr::mapvalues(x$probe, row.names(corlist$corsum), corlist$corsum$p, warn_missing = F)
    return(x)
}
