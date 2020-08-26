#' plot_eqtl
#'
#' This function will draw plots for individual cases
#' @param x dataframe of eqtl result
#' @param snpid string snp
#' @param probeid string cpg id
#' @param geneid string gene
#' @param plot_what string what to plot, one of encode_correlation, all
#' @return will generate plot
#' @import tidyverse
#' @import assertthat
#' @importFrom logging loginfo
#' @importFrom RNOmni rankNorm
#' @importFrom sjPlot plot_model
#' @import Sushi
#' @examples
#' load("data/res.Rdata")
#' res <- add_siglabel_class(res)
#' @export
#'
plot_eqtl <- function(x, snpid, probeid, geneid, plot_what="all"){
    if (plot_what == "all"){
        outfile_name <- paste0("summary_",snpid,":",probeid,":",geneid,".pdf")
        pdf(outfile_name, useDingbats = F)
        par(mfrow=c(2,3))
    }
    ctcfid <- x %>% filter(probe == probeid) %>%  dplyr::slice(1) %>% pull(ctcf)
    logging::loginfo(paste("ctcfid:", ctcfid))
    extr_encode_data <- function(probeid, file){
        t <- system(paste("fgrep -w", probeid, file), intern = T)
        return(as.numeric(unlist(strsplit(t, "\t"))[2:length(unlist(strsplit(t, "\t")))]))
    }
    if (plot_what == "all" | plot_what == "encode_correlation"){
        corest <- x %>% filter(probe == probeid) %>%  dplyr::slice(1) %>% pull(encodeCorEst)
        corest <- as.numeric(corest)
        assertthat::assert_that(is.character(ctcfid))
        m <- extr_encode_data(probeid, "data/encode_methylation_matrix.txt")
        tr <- extr_encode_data(probeid, "data/encode_totalRead_matrix.txt")
        ctcf <- extr_encode_data(ctcfid, "data/encode_ctcf_matrix.txt")
        plot(m[tr>=20], ctcf[tr>=20], cex.axis = 0.8, cex.lab = 1, cex.main = 1, xlab = probeid, ylab = ctcfid, frame = F)
        legend("topright", bty = "n", lty = 0,
               legend = paste0("PCC:", round(corest, 3)), cex = 0.8)
    }
    if (plot_what == "all" | plot_what == "pos_cpg"){
        t <- x %>% filter(ctcf==ctcfid, sigdelta, classSimplified!="None")
        t$snpDist <- t$snpPos - t$probePos
        t$snpDist <- t$snpDist/1000
        t$geneDist <- (t$tss - t$probePos)/1000
        maxDist <- max(c(t$snpDist, t$geneDist))
        minDist <- min(c(t$snpDist, t$geneDist))
        t$color <- ifelse(t$classSimplified=="Insulating", "coral2","lightblue4")
        plot(t$geneDist, t$snpDist, col = t$color, cex = -log2(t$pdelta)/max(-log2(t$pdelta)),
             pch = 3, main = ctcfid, xlab = "TSS distance from CpG (Kbp)", ylab = "SNP distance from CpG (Kbp)",
             cex.main = 1, cex.lab = 1, cex.axis = 0.8, frame = F,
             xlim=c(minDist,maxDist), ylim=c(minDist, maxDist))
        #text(t$geneDist, t$snpDist, labels = t$gene,cex = 0.8)
        abline(h=0, lty = 2, col = "grey60")
        abline(v=0, lty = 2, col = "grey60")
        legend("topright", bty = "n", col=(c("coral2","lightblue4")), legend = c("eQTL in high 5mc","eQTL in low 5mc"), cex = 0.8, lty = c(1,1))
        legend("bottomleft", bty = "n", legend = c("*Size refers to",expression(-log[2]*" P value")), cex = 0.6, lty = 0)
    }
    if (plot_what == "all" | plot_what == "pos_gene"){
        t <- x %>% filter(ctcf==ctcfid, sigdelta, classSimplified!="None")
        t$snpDist <- t$snpPos - t$tss
        t$snpDist <- t$snpDist/1000
        t$probeDist <- t$probePos - t$tss
        t$probeDist <- t$probeDist/1000
        maxDist <- max(c(t$snpDist, t$probeDist))
        minDist <- min(c(t$snpDist, t$probeDist))
        t$color <- ifelse(t$classSimplified=="Insulating", "coral2","lightblue4")
        plot(t$probeDist, t$snpDist, col = t$color, cex = -log2(t$pdelta)/max(-log2(t$pdelta)),
             pch = 3, main = geneid, xlab = "CpG distance from TSS (Kbp)", ylab = "SNP distance from TSS (Kbp)",
             cex.main = 1, cex.lab = 1, cex.axis = 0.8, frame = F,
             xlim=c(minDist,maxDist), ylim=c(minDist, maxDist))
        #text(t$probeDist, t$snpDist, labels = t$gene,cex = 0.8)
        abline(h=0, lty = 2, col = "grey60")
        abline(v=0, lty = 2, col = "grey60")
        legend("topright", bty = "n", col=(c("coral2","lightblue4")), legend = c("eQTL in high 5mc","eQTL in low 5mc"), cex = 0.8, lty = c(1,1))
        legend("bottomleft", bty = "n", legend = c("*Size refers to",expression(-log[2]*" P value")), cex = 0.6, lty = 0)
    }
    if (plot_what == "all" | plot_what == "eqtl"){
        m <- extr_encode_data(probeid, "data/CPGEA_mdata.normal.txt")
        Methylation <- ifelse(m<=median(m), "Low","High")
        e <- extr_encode_data(geneid, "data/CPGEA_edata.normal.txt")
        e <- RNOmni::rankNorm(e)
        g <- extr_encode_data(snpid, "data/CPGEA_gdata.txt")
        info <- x %>% filter(probe==probeid, gene==geneid, snp==snpid) %>% dplyr::select(pmlow,fdrlow,pmhigh,fdrhigh,pdelta)
        t <- boxplot(e[Methylation=="Low"] ~ g[Methylation=="Low"], xaxt = "n", frame = F, ylab = geneid, xlab = snpid,
                     cex.lab = 1, cex.axis = 0.8, cex.main = 1, main = "Low 5mc", sub = paste("p =", round(info$pmlow[1], 3), "FDR =", round(info$fdrlow[1], 3)))
        axis(1, at = as.numeric(t$names)+1, labels = paste(t$names, paste0("(n=",t$n,")"), sep = "\n"), cex = 0.8, tick = F)
        t <- boxplot(e[Methylation=="High"] ~ g[Methylation=="High"], xaxt = "n", frame = F, ylab = geneid, xlab = snpid,
                     cex.lab = 1, cex.axis = 0.8, cex.main = 1, main = "High 5mc", sub = paste("p =", round(info$pmhigh[1], 3), "FDR =", round(info$fdrhigh[1],3)))
        axis(1, at = as.numeric(t$names)+1, labels = paste(t$names, paste0("(n=",t$n,")"), sep = "\n"), cex = 0.8, tick = F)
        #f <- lm(e ~ g*Methylation)
        #plot(effects::allEffects(f, multiline = T, ci.style="bands"))
        #sjPlot::plot_model(f, type = "pred", terms = c("g", "Methylation"), axis.title = c(snpid, geneid),
        #                   cex.main = 1, cex.axis = 0.8, cex.lab = 1, title = "me-eqtl")
    }
    if (plot_what == "all" | plot_what == "meth_density"){
        m <- extr_encode_data(probeid, "data/CPGEA_mdata.normal.txt")
        plot(density(m), xlim=c(0,1), frame = F,
             cex.lab = 1, cex.axis = 0.8, cex.main = 1, main = probeid)
        abline(v=median(m), lty = 2, col = 2)
    }

    if(plot_what == "all"){
        dev.off()
    }
}
