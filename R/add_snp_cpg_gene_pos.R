#' add_snp_cpg_gene_pos
#'
#' This function will add POS information of SNP/ATCA, CpG/CTCF and genes/TSS
#' And split the res to two parts sign of cor(CpG, CTCF)
#' @param file_res path to me_eqtl_res.Rdata
#' @param file_cpg_ctcf file of significant correlated CpG and CTCF
#' @param file_snp_loci path to SNP bed file
#' @param file_atac_peak file of ATAC-seq peaks with SNPs
#' @param file_tss_loci path to TSS bed file
#' @param name prefix of output files
#' @import dplyr
#'

add_snp_cpg_gene_pos <- function(file_res, file_cpg_ctcf, file_snp_loci, file_atac_peak, file_tss_loci, name)
{

    load(file_res)   ## load me_eqtl_res.Rdata with variable me_eqtl_res

    ## add cor(CpG, CTCF) info
    cor_cpg_ctcf <- read.table(file_cpg_ctcf, header = T, as.is = T)
    idx_cpg <- match(me_eqtl_res$cpg, cor_cpg_ctcf$cpg)
    cor_cpg_ctcf_add <- cor_cpg_ctcf[idx_cpg, ] %>%  dplyr::select(ctcf_id, cor_cpg_ctcf = pcc_2, cor_cpg_ctcf_fdr = qval_global,

                                                                                                                               chr, cpg_pos, ctcf_mid)
    ## add snp info
    snp_info <- read.table(file_snp_loci, header = F, as.is = T)
    idx_snp <- match(me_eqtl_res$snp, snp_info$V4)
    snp_info_add <- snp_info[idx_snp, ] %>%  dplyr::select(snp_pos = V3, snp_ref = V7, snp_alt = V8)

    ## add ATAC-seq peak info
    atac_info <- read.table(file_atac_peak, as.is = T)
    idx_atac <-  match(me_eqtl_res$snp, atac_info$V4)
    atac_info_add <- atac_info[idx_atac, ] %>%
                     mutate(snp2atac_mid = V3 - (V10 + V11)/2) %>%
                     dplyr::select(atac_id = V12, atac_start = V10, atac_end = V11, snp2atac_mid)

    ## add gene info
    gene_info <- read.table(file_tss_loci, header = F, as.is = T)
    idx_gene <- match(me_eqtl_res$gene, gene_info$V4)
    gene_info_add <- gene_info[idx_gene, ] %>%  dplyr::select(tss_pos = V3, gene_type = V5, gene_strand = V6)

    ## merging
    res <- cbind(me_eqtl_res, cor_cpg_ctcf_add, snp_info_add, atac_info_add,  gene_info_add)

    ###################################################
    # split res based on pos and neg cor(CpG, CTCF)
    ###################################################
    idx_pos <- res$cor_cpg_ctcf > 0
    idx_neg <- res$cor_cpg_ctcf < 0

    res_pos <- res[idx_pos, ]
    res_neg <- res[idx_neg, ]

    write.csv(res_pos, file = paste0(name, "_with_pos_cor_cpg_ctcf.csv"))
    write.csv(res_neg, file = paste0(name, "_with_neg_cor_cpg_ctcf.csv"))
}
