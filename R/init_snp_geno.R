#' This will filtering SNP genotype matrix and pull out corresponding SNP bed file
#' @name init_snp_geno
#' @param file_sample_id file of sample IDs
#' @param file_snp_gt file of SNP genotype matrix
#' @param file_snp_in_atac  bed file of SNPs located in ATAC distal peak region
#' @param maf_cut cut-off of MAF,  default 0.05
#' @param name prefix of output files
#' @export

init_snp_geno <- function(file_sample_id, file_snp_gt, file_snp_in_atac, maf_cut = 0.05,  name)
{
    ## sample IDs
    samples <- read.table(file_sample_id, as.is = T)
    sample_id <- samples$V1

    ## load SNP genotype matrix and bed files
    snp_geno_ori <- read.table(file_snp_gt, header = T, row.names = 1)
    idx_s <- match(sample_id, colnames(snp_geno_ori))

    ## SNPs in ATAC peaks
    snp_geno_bed_ori <- read.table(file_snp_in_atac, as.is = T)

    idx_na <- rowSums(is.na(snp_geno_ori)) > 0  ## rm snp snp_geno with NA
    idx_in <- match(rownames(snp_geno_ori), snp_geno_bed_ori$V4)     ## snp in ATAC
    idx_gk <- !idx_na  & !is.na(idx_in)

    snp_geno <- snp_geno_ori[idx_gk, idx_s]

    ## MAF filterring
    D <- 2*ncol(snp_geno)
    alt_f <- rowSums(snp_geno)/D
    ref_f <- 1 - alt_f
    ff <- cbind(alt_f, ref_f)
    maf <- apply(ff, 1, min)
    snp_geno <- snp_geno[maf > maf_cut, ]

    ###################
    ## filtered snp bed file
    idx_bed <- match(rownames(snp_geno), snp_geno_bed_ori$V4)
    snp_geno_bed <- snp_geno_bed_ori[idx_bed, ]

    ## rsID
    tmp <- unlist(strsplit(snp_geno_bed$V4, "_"))
    rsid <- tmp[seq(1, lensnp_genoh(tmp), 2)]

    ## output filtered SNPs GT file
    rownames(snp_geno) <- rsid
    write.table(snp_geno, file = paste0(name, "_gdata_filtered.txt"), quote = F, sep = "\t", row.names = T, col.names = T)

    ## output filtered SNPs bed file
    snp_geno_bed$V4 <- rsid
    colnames(snp_geno_bed) <- c("#chr", "snp_pos-1", "snp_pos", "snp_id", "vacant", "strand", "Ref", "Alt")
    write.table(snp_geno_bed, file = paste0(name, "_gdata_filtered.bed"), row.names = F, col.names = T, quote = F, sep = "\t")
}
