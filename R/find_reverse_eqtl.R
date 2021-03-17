#' find_reverse_eqtl
#'
#' This will find eQTLs for a subset with a particular methylation probe
#' which showing oppostie direction
#' @name find_reverse_eqtl
#' @param x eqtl dataframe, e.g., res$gene
#' @param file_edata path to edata file
#' @param file_gdata path to gdata file
#' @param file_mdata path to mdata file
#' @param probeid string probeid
#' @param snpid string snpid
#' @param cutAtMedian boolean to indicaite if methylation should be stratified at median, defailt T
#' @param writetable boolean to indicate if track files be saved or not
#' @param methylation_status to calculate eQTLs in either High or Low methylation level of the probeid
#' @importFrom RNOmni rankNorm
#' @export
#'
#'
find_reverse_eqtl <- function(x, file_edata, file_gdata, file_mdata, probeid, snpid, methylation_status, cutAtMedian = T, writetable = F)
{
    snppos <- x[x$snp == snpid, 'snpPos'][1]
    chr <- x[x$snp == snpid, 'seqnames'][1]

    ## read in data
    e <- read.table(file_edata, header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
    g <- read.table(file_gdata, header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
    m <- read.table(file_mdata, header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
    m <- as.numeric(m[probeid,])

    if (cutAtMedian){
        m <- ifelse(m<=median(m), "Low","High")
    } else {
        m <- ifelse(m<=summary(m)[2], "Low",
                ifelse(m>=summary(m)[5], "High", "Mid"))
    }
    g <- as.numeric(g[snpid,])

    tssfile <- read.table("data/CPGEA_edata.normal.TSS.bed", header = F, sep = "\t", stringsAsFactors = F)
    row.names(tssfile) <- tssfile$V4

    neighbor_genes <- tssfile %>%
                      filter(V1==chr, V2 >= (snppos-500000), V2 <= (snppos+500000)) %>%
                      pull(V4)

    e <- e[neighbor_genes,m==methylation_status]
    g <- g[m==methylation_status]

    bigint <- data.frame(chrom = chr,
                         chromStart = NA,
                         chromEnd = NA,
                         name = ".",
                         score = 0,
                         value = NA,
                         exp = "sigeQTL",
                         color = NA,
                         sourceChrom = chr,
                         sourceStart = snppos,
                         sourceEnd = snppos +1,
                         sourceName = snpid,
                         sourceStrand = ".",
                         targetChrom = NA,
                         targetStart = NA,
                         targetEnd = NA,
                         targetName = NA,
                         targetStrand = ".")


    for (i in 1:nrow(e))
    {
        enorm <- RNOmni::rankNorm(as.numeric(e[i,]))
        f <- try(lm(enorm ~ g))
        if (class(f) != "lm") next
        p <- as.numeric(summary(f)$coefficients[2,4])
        est <- as.numeric(coef(f)[2])
        message(paste(p, est))
        bigint <- rbind(bigint,
                        data.frame(chrom = chr,
                                   chromStart = min(snppos, tssfile[neighbor_genes[i], 'V2']),
                                   chromEnd = max(snppos, tssfile[neighbor_genes[i], 'V2']),
                                   name = ".",
                                   score = ceiling(-log2(p)),
                                   value = -log2(p),
                                   exp = "sigeQTL",
                                   color = ifelse (est>=0, "#EE6A50","#68838B"),
                                   sourceChrom = chr,
                                   sourceStart = snppos,
                                   sourceEnd = snppos +1,
                                   sourceName = paste(rev(rev(strsplit(snpid,"_")[[1]])[-1]), collapse = ":"),
                                   sourceStrand = ".",
                                   targetChrom = chr,
                                   targetStart = tssfile[neighbor_genes[i], 'V2'],
                                   targetEnd = tssfile[neighbor_genes[i], 'V2'] + 1,
                                   targetName = neighbor_genes[i],
                                   targetStrand = "."))
    }
    bigint <- bigint[-1,]


    if (writetable){
        write(paste("track type=interact name=\"eQTLs for", snpid,"in",methylation_status,"methylation level of",probeid,"\" interactDirectional=false maxHeightPixels=200:100:50 visibility=full"), file = paste0("data/int_",probeid,"_",methylation_status,"_bigInteract.txt"))
        write("browser position chr16:89590592-89519592", file = paste0("data/int_",probeid,"_",methylation_status,"_bigInteract.txt"), append = T)
        write.table(bigint, file = paste0("data/int_",probeid,"_",methylation_status,"_bigInteract.txt"), row.names = F, quote = F, col.names = F, sep = "\t", append = T)
    }

    return(bigint)
}
