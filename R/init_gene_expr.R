#' This will filtering expression matrix and pull out corresponding gene TSS bed file
#' @name init_gene_expr
#' @param file_sample_id file of sample IDs
#' @param file_gene_expr file of gene expression matrix
#' @param file_gene_gtf  file of gene annotation
#' @param expr_med cut-off of median gene expression,  default 1
#' @param name prefix of output files
#' @export

init_gene_expr <- function(file_sample_id, file_gene_expr, file_gene_gtf, expr_med = 1,  name)
{
    ## sample IDs
    samples <- read.table(file_sample_id, as.is = T)
    sample_id <- samples$V1

    ## reading gene expression and gene annotation
    gene_expr <- read.table(file_gene_expr, header = T, as.is = T)
    gene_gtf <- read.table(file_gene_gtf, as.is = T)
    idx_s <- match(sample_id, colnames(gene_expr))

    ## remove gene_name with tow ensemble ID
    idx_dup <- table(gene_expr[, 2]) >= 2
    gene_dup <- names(table(gene_expr[, 2]))[idx_dup]
    idx_rm <- match(gene_expr[, 2], gene_dup)

    gene_expr_m <- gene_expr[is.na(idx_rm), idx_s]
    rownames(gene_expr_m) <- gene_expr[, 2][is.na(idx_rm)]

    ## checking gene types and keep protein_coding, lincRNA
    idx_m <- match(rownames(gene_expr_m), gene_gtf$V4)
    gene_type <- gene_gtf$V7[idx_m]
    idx_tk <- gene_type == "protein_coding" | gene_type == "lincRNA"
    idx_tk[is.na(idx_tk)] <- FALSE

    ## filtering based on the median expressioin level > expr_med and gene type
    idx_re <- apply(gene_expr_m, 1, median) > expr_med & idx_tk
    gene_expr_f <- gene_expr_m[idx_re, ]

    write.table(gene_expr_f, file = paste0(name, "_edata_filtered.txt"), quote = F, sep = "\t", row.names = T, col.names = T)

    ##########################################
    ##### TSS site bed file for filtered genes
    idx_g <- match(rownames(gene_expr_f), gene_gtf$V4)
    tss_bed <- gene_gtf[idx_g, c(1:4, 7, 6)]
    idx_plus <- tss_bed$V6 == "+"

    ## plus strand
    tss_bed[idx_plus, 3] <-  tss_bed[idx_plus, 2]
    tss_bed[idx_plus, 2] <-  tss_bed[idx_plus, 2] - 1

    ## minus strand
    tss_bed[!idx_plus, 2] <-  tss_bed[!idx_plus, 3] - 1

    colnames(tss_bed) <- c("#chr", "TSS_pos-1", "TSS_pos", "Gene", "Gene_type", "strand")
    write.table(tss_bed, file = paste0(name, "_edata_filtered_TSS.bed"), quote = F, sep = "\t", row.names = F, col.names = T)
}
