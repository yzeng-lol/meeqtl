#' plot_me_eQTLs
#'
#' This function will draw combined plots inclduing me_eQTLs, eQTL and meQTL for snp_cpg_gene combinations
#' @param snp_cpg_gene vector of snp_cpg_gene comination
#' @param file_edata path to edata file
#' @param file_gdata path to gdata file
#' @param file_mdata path to mdata file
#' @param file_tss_loci path to TSS bed file
#' @param file_snp_loci path to SNP bed file
#' @param file_res path to me_eqtl_res.Rdata
#' @param name output plot file name prefix
#' @import ggplot2
#' @import RColorBrewer
#' @importFrom logging loginfo
#' @examples
#' snp_cpg_gene <- c("rs28520918_chr6:31402506_POU5F1", "rs1039750_chr11:94603277_SESN3")
#' plot_me_eQTLs(snp_cpg_gene, file_edata, file_gdata, file_mdata, file_snp_loci, file_res, "me_eQTLs_test")

####################
## plotting function
####################

.plot_cmb <- function(snp_geno, cpg_meth, gene_expr, effs, pvals, name)
{
    ## split to 5mC high and low group
    ## rmove sample equal to median cpg_meth
    idx_mh <- cpg_meth > median(cpg_meth)
    idx_mh[is.na(idx_mh)] <- FALSE

    idx_ml <- cpg_meth < median(cpg_meth)
    idx_ml[is.na(idx_ml)] <- FALSE

    ## expr range to make sure using equal y axis
    ymax <- max(log2(gene_expr + 1))
    ymin <- min(log2(gene_expr + 1))

    y_range <- ymax - ymin
    ymax <- ymax + y_range * 0.05
    ycnt <- ymin - y_range * 0.05       # geno_cnt y position
    ymin <- ymin - y_range * 0.1

    #############################
    ## me_eQTL for High 5mC group
    {
        dat <- data.frame(geno = snp_geno[idx_mh], trait = log2(gene_expr[idx_mh] + 1))

        ## number of sample for each genotyp
        geno_cnt <- table(dat$geno)

        gh <- ggplot(dat, aes(x = geno, y = trait, color = geno)) + geom_boxplot(outlier.shape = NA)
        gh <- gh + geom_jitter() + scale_color_brewer(palette="Dark2") + ylim(ymin, ymax)
        gh <- gh + labs(title = name, y = "log2(FPKM +1)", x = "SNP Genotype",
                        subtitle = paste0("High_5mC: Effect size = ", effs[1], "; P_val = ", pvals[1]))
        gh <- gh + geom_text(x = 1, y = ycnt , label = paste0("N = ", geno_cnt[1]), cex = 3, color = "gray50")
        gh <- gh + geom_text(x = 2, y = ycnt , label = paste0("N = ", geno_cnt[2]), cex = 3, color = "gray50")
        gh <- gh + geom_text(x = 3, y = ycnt , label = paste0("N = ", geno_cnt[3]), cex = 3, color = "gray50")
        gh <- gh + theme_bw() +  theme(legend.position = "none")
    }

    #############################
    ## me_eQTL for Low 5mC group
    {
        dat <- data.frame(geno = snp_geno[idx_ml], trait = log2(gene_expr[idx_ml] + 1))

        ## number of sample for each genotyp
        geno_cnt <- table(dat$geno)

        gl <- ggplot(dat, aes(x = geno, y = trait, color = geno)) + geom_boxplot(outlier.shape = NA)
        gl <- gl + geom_jitter() + scale_color_brewer(palette="Dark2") + ylim(ymin, ymax)
        gl <- gl + labs(title =  name, y = "log2(FPKM +1)", x = "SNP Genotype",
                        subtitle = paste0("Low_5mC: Effect size = ", effs[2], "; P_val = ", pvals[2]))
        gl <- gl + geom_text(x = 1, y = ycnt , label = paste0("N = ", geno_cnt[1]), cex = 3, color = "gray50")
        gl <- gl + geom_text(x = 2, y = ycnt , label = paste0("N = ", geno_cnt[2]), cex = 3, color = "gray50")
        gl <- gl + geom_text(x = 3, y = ycnt , label = paste0("N = ", geno_cnt[3]), cex = 3, color = "gray50")
        gl <- gl + theme_bw() +  theme(legend.position = "none")
    }

    ###########
    ## for eQTL
    {
        dat <- data.frame(geno = snp_geno, trait = log2(t(gene_expr) + 1))
        colnames(dat) <- c("geno", "trait")

        ## number of sample for each genotyp
        geno_cnt <- table(dat$geno)

        ge <- ggplot(dat, aes(x = geno, y = trait, color = geno)) + geom_boxplot(outlier.shape = NA)
        ge <- ge + geom_jitter() + scale_color_brewer(palette="Dark2") + ylim(ymin, ymax)
        ge <- ge + labs(title = name, y = "log2(FPKM +1)", x = "SNP Genotype",
                        subtitle = paste0("eQTL: Effect size = ", effs[3], "; P_val = ", pvals[3]))
        ge <- ge + geom_text(x = 1, y = ycnt , label = paste0("N = ", geno_cnt[1]), cex = 3, color = "gray50")
        ge <- ge + geom_text(x = 2, y = ycnt , label = paste0("N = ", geno_cnt[2]), cex = 3, color = "gray50")
        ge <- ge + geom_text(x = 3, y = ycnt , label = paste0("N = ", geno_cnt[3]), cex = 3, color = "gray50")
        ge <- ge + theme_bw() +  theme(legend.position = "none")
    }

    ############
    ## for meQTL
    {
        dat <- data.frame(geno = snp_geno, trait = t(cpg_meth))
        colnames(dat) <- c("geno", "trait")

        ## number of sample for each genotyp
        geno_cnt <- table(dat$geno)

        ## beta range to make sure using equal y axis
        ymax <- max(cpg_meth, na.rm = T)
        ymin <- min(cpg_meth, na.rm = T)

        y_range <- ymax - ymin
        ymax <- ymax + y_range * 0.05
        ycnt <- ymin - y_range * 0.05       # geno_cnt y position
        ymin <- ymin - y_range * 0.1


        gme <- ggplot(dat, aes(x = geno, y = trait, color = geno)) + geom_boxplot(outlier.shape = NA)
        gme <- gme + geom_jitter() + scale_color_brewer(palette="Dark2") + ylim(ymin, ymax)
        gme <- gme + labs(title = name, y = "CpG Methylation Beta value", x = "SNP Genotype",
                          subtitle = paste0("meQTL: Effect size = ", effs[4], "; P_val = ", pvals[4]))
        gme <- gme + geom_text(x = 1, y = ycnt , label = paste0("N = ", geno_cnt[1]), cex = 3, color = "gray50")
        gme <- gme + geom_text(x = 2, y = ycnt , label = paste0("N = ", geno_cnt[2]), cex = 3, color = "gray50")
        gme <- gme + geom_text(x = 3, y = ycnt , label = paste0("N = ", geno_cnt[3]), cex = 3, color = "gray50")
        gme <- gme + theme_bw() +  theme(legend.position = "none")
    }

    #####################
    ## put them together

    g_pool <- ggarrange(gh, gl, ge, gme, ncol = 2, nrow = 2)
    ggsave(paste0(name, "_me_eQTL.pdf"), width = 8, height = 8)

}


############################
## loading data and plotting
############################
{
    ## load the genotype, methylation, and gene expression data and me_eQTLs results
    edata <- read.table(file_edata, header = T)
    gdata <- read.table(file_gdata, header = T)
    mdata <- read.table(file_mdata, header = T)
    snp_bed <- read.table(file_snp_loci, as.is = T)

    load("/Users/Yong/Yong/CTCF_methy_QTL/changhai/normal/1_run_meeqtl/me_eqtl_res.Rdata")

    ## putll out data for snp_cpg_gene
    L <- length(snp_cpg_gene)
    for(i in 1:L)
    {
        ## res from matrixeQTL
        name <- snp_cpg_gene[i]
        idx_res <- match(name, rownames(me_eqtl_res))

        ## effect sizes
        effs <- c(round(me_eqtl_res$me_high_eqtl_beta[idx_res], 2),
                  round(me_eqtl_res$me_low_eqtl_beta[idx_res], 2),
                  round(me_eqtl_res$eqtl_beta[idx_res], 2),
                  round(me_eqtl_res$meqtl_beta[idx_res], 2))

        ## QTLs pvalues
        pvals <- c(formatC(me_eqtl_res$me_high_eqtl_pvalue[idx_res], format = "e", digits = 2),
                   formatC(me_eqtl_res$me_low_eqtl_pvalue[idx_res], format = "e", digits = 2),
                   formatC(me_eqtl_res$eqtl_pvalue[idx_res], format = "e", digits = 2),
                   formatC(me_eqtl_res$meqtl_pvalue[idx_res], format = "e", digits = 2))

        ## SNP, CpG and Gene IDs
        snp_id <- strsplit(snp_cpg_gene[i], "_")[[1]][1]
        cpg_id <- strsplit(snp_cpg_gene[i], "_")[[1]][2]
        gene_id <- strsplit(snp_cpg_gene[i], "_")[[1]][3]

        ## to genotype
        idx_snp <- match(snp_id, rownames(gdata))
        idx_snp_bed <- match(snp_id, snp_bed$V4)
        homo_ref <- paste0(snp_bed$V7[idx_snp_bed], snp_bed$V7[idx_snp_bed])
        homo_alt <- paste0(snp_bed$V8[idx_snp_bed], snp_bed$V8[idx_snp_bed])
        hete  <- paste0(snp_bed$V7[idx_snp_bed], snp_bed$V8[idx_snp_bed])

        #snp_geno <- factor(gdata[idx_snp, ], labels = c(homo_ref, hete, homo_alt))'
        ## could be only two levels
        snp_geno <- factor(gdata[idx_snp, ])
        levels(snp_geno)[levels(snp_geno) == "0"] <- homo_ref
        levels(snp_geno)[levels(snp_geno) == "1"] <- hete
        levels(snp_geno)[levels(snp_geno) == "2"] <- homo_alt

        ## methylation level
        idx_cpg <- match(cpg_id, rownames(mdata))
        cpg_meth <- as.matrix(mdata[idx_cpg, ])

        ## gene expression level
        idx_gene <- match(gene_id, rownames(edata))
        gene_expr <- as.matrix(edata[idx_gene, ])

        ## prefix of output file
        out_prefix <- name

        ## ploting
        .plot_cmb(snp_geno, cpg_meth, gene_expr, effs, pvals, out_prefix)

    }

}

