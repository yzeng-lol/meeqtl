#' This will filtering CpG methylation matrix and pull out corresponding CpG and CTCF bed file
#' @name init_cpg_meth
#' @param file_sample_id file of sample IDs
#' @param file_cpg_meth file of merged CpG sits methylation
#' @param file_cpg_ctcf file of significant correlated CpG and CTCF
#' @param meth_l cut-off of low methylation level,  default 0.25
#' @param meth_h cut-off of high methylation level,  default 0.75
#' @param name prefix of output files
#' @import dplyr

init_cpg_meth <- function(file_sample_id, file_cpg_meth, file_cpg_ctcf, meth_l = 0.25, meth_h = 0.75,  name)

{
    ####################################################
    ## tranfer merged CpG methyaltion files to matirx
    {
    ## sample IDs
    samples <- read.table(file_sample_id, as.is = T)
    sample_id <- samples$V1
    M <- length(sample_id)

    ## cpg_meth
    cpg_meth_raw <- read.table(file_cpg_meth, as.is = T)
    cpg_id <- paste(cpg_meth_raw$V1, cpg_meth_raw$V2, sep = ":")
    N <- length(cpg_id)

    me_beta <- matrix(0, N, M)            ## cpg methylation ratio
    me_tr <- matrix(0, N, M)              ## cpg total reads
    rownames(me_beta) <- rownames(me_tr) <- cpg_id
    colnames(me_beta) <- colnames(me_tr) <- sample_id

    ## spliting samples, methylaton levels
    sample_sp <- strsplit(cpg_meth_raw$V4, ",")
    beta_sp <- strsplit(cpg_meth_raw$V5, ",")
    tr_sp <- strsplit(cpg_meth_raw$V6, ",")

    for(i in 1:N)
    {
        idx <- match(sample_sp[[i]], sample_id)
        idx_s <- idx[!is.na(idx)]
        me_beta[i, idx_s] <- round(as.numeric(beta_sp[[i]])[!is.na(idx)], 2)         ## round 2 digits
        me_tr[i, idx_s] <- as.numeric(tr_sp[[i]])[!is.na(idx)]
    }

    write.table(me_beta, file = paste0(name, "_mdata.txt"), row.names = T, col.names = T, quote = F, sep = "\t")
    }

    ##############
    ### filtering
    {
    ### cpg_CTCF information in encoden
    ## rm cpgs pulled out while no in sigCor
    cpg_ctcf <- read.table(file_cpg_ctcf, header = T)
    idx_f <- match(rownames(me_beta), cpg_ctcf$cpg)
    idx_me <- r_median > meth_l & r_median < meth_h       ## filter out less vary CpGs
    idx_k <- idx_me & !is.na(idx_f)
    write.table(me_beta[idx_k, ], file = paste0(name, "_mdata_filtered.txt"), row.names = T, col.names = T, quote = F, sep = "\t")

    ## distribution of CpG methylation
    r_median <- apply(me_beta, 1, median)
    r_mean <- rowMeans(me_beta)

    pdf(paste0(name, "_mdata_filtered_median_density.pdf"), width = 4, height = 4)
    plot(density(r_median[idx_k]), main = "Median_CpG_methylation_level")
    dev.off()

    pdf(paste0(name, "_mdata_filtered_mean_density.pdf"), width = 4, height = 4)
    plot(density(r_mean[idx_k]), main = "Mean_CpG_methylation_level")
    dev.off()


    ### cpg_CTCF information in encoden
    me_f <- rownames(me_beta[idx_k, ])
    idx_out <- match(me_f, cpg_ctcf$cpg)
    cpg_ctcf_f <- cpg_ctcf[idx_out, ]

    ## cpg bed file
    cpg_bed <- cpg_ctcf_f %>%  select(chr, cpg_start, cpg_end, cpg, cpg_pos) %>%  mutate(strand = ".")
    ctcf_bed <- cpg_ctcf_f %>% select(chr, ctcf_start, ctcf_end, ctcf_id, ctcf_mid) %>% mutate(strand = ".")

    colnames(cpg_bed)[1] <- paste0("#", colnames(cpg_bed)[1])
    colnames(ctcf_bed)[1] <- paste0("#", colnames(ctcf_bed)[1])

    write.table(cpg_bed, file =  paste0(name, "_mdata_filtered_CpG.bed"), col.names = T, row.names = F, quote = F, sep = "\t")
    write.table(ctcf_bed, file = paste0(name, "_mdata_filtered_CTCF.bed"), col.names = T, row.names = F, quote = F, sep = "\t")



    }
}
