#' find_eqtl_tests
#'
#' This will run matrixeQTL for methylation-dependent eQTL
#' @name run_matrixeqtl
#' @param x eqtl dataframe, e.g., res$gene
#' @param file_edata path to edata file
#' @param file_gdata path to gdata file
#' @param file_mdata path to mdata file
#' @param file_snp_loci path to snp bed file
#' @param file_tss path to TSS bed file
#' @param geneSnpMaxDistance max distance between snp and gene, default 500000
#' @param cpgSnpMaxDistance max distance between snp and probe, default 500000
#' @param file_atac path to overlappiing file of snp and ATAC-seq peaks
#' @param file_ctcf path to overlappiing file of cpgs and ctcf peaks in ENCODE
#' @import GenomicRanges
#' @import bigstatsr
#' @import assertthat
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom rtracklayer import
#' @importFrom ChIPpeakAnno reCenterPeaks
#' @importFrom logging loginfo
#' @export
#'
#'
run_matrixeqtl <- function(file_edata, file_gdata, file_mdata, file_snp_loci,
                           file_tss, file_atac, file_ctcf,
                           geneSnpMaxDistance = 500000, cpgSnpMaxDistance = 500000)
{
###########################################
## combination for snp-CpG-gene for testing
###########################################
    {
    ## read in expression, genotype, methylation matrix ...
    logging::loginfo("Reading in files ...")
    edata <- read.table(file_edata, header = T, sep = "\t", stringsAsFactors = F)
    gdata <- read.table(file_gdata, header = T, sep = "\t", stringsAsFactors = F)
    mdata <- read.table(file_mdata, header = T, sep = "\t", stringsAsFactors = F)
    snp_atac <- read.table(file_atac, sep = "\t", header = F, stringsAsFactors = F)
    cpg_ctcf <- read.table(file_ctcf, sep = "\t", header = T, stringsAsFactors = F)

    ## CpG sits filtering based on the median beta values
    t <- mdata[,-1]
    logging::loginfo(paste("CpGs found:", nrow(t)))
    mdata <- mdata[apply(t, 1, median)>0.25 & apply(t, 1, median)<0.75, ]

    t <- mdata[,-1]
    row.names(t) <- mdata$cpgid
    logging::loginfo(paste("CpGs retained:", nrow(mdata)))

    ## Extenting to CpG site centered "cpgSnpMaxDistance * 2" region
    cpg_loci <- GRanges(unlist(lapply(mdata$cpgid, function(x) strsplit(x, ":")[[1]][1])),
                        IRanges(as.numeric(unlist(lapply(mdata$cpgid, function(x) strsplit(x, ":")[[1]][2]))),
                        width = 1, names = mdata$cpgid))
    cpg_loci <- ChIPpeakAnno::reCenterPeaks(cpg_loci, cpgSnpMaxDistance * 2)
    start(cpg_loci)[start(cpg_loci)<0] <- 0

    ## read in genotype information
    snp_loci <- read.table(file_snp_loci, sep = "\t", header = F, stringsAsFactors = F)
    snp_loci <- GRanges(snp_loci$V1, IRanges(snp_loci$V2, snp_loci$V3, names = snp_loci$V4))
    logging::loginfo(paste("SNPs found:", length(snp_loci)))

    ## intersect SNPs with extented CpG region
    snp_loci <- subsetByOverlaps(snp_loci, cpg_loci)
    snppos <- data.frame(snp = names(snp_loci), chr = seqnames(snp_loci), pos = start(snp_loci)) #for matrixeqtl
    logging::loginfo(paste("SNPs within extented CpG regions:", length(snp_loci)))


    tss <- read.table(file_tss, sep = "\t", header = F, stringsAsFactors = F)
    tss <- GRanges(tss$V1, IRanges(tss$V2, tss$V3, names = tss$V4))
    geneSnpWindow <- ChIPpeakAnno::reCenterPeaks(snp_loci, geneSnpMaxDistance * 2)
    start(geneSnpWindow)[start(geneSnpWindow)<0] <- 0
    logging::loginfo(paste(length(tss),"genes found"))
    tss_within_dist <- as.data.frame(findOverlaps(tss, geneSnpWindow))
    snps_genes_assoc <- split(tss_within_dist$queryHits, as.character(tss_within_dist$subjectHits))
    names(snps_genes_assoc) <- names(snp_loci[as.numeric(names(snps_genes_assoc))])
    snps_genes_assoc <- lapply(snps_genes_assoc, function(x) names(tss[x]))

    genepos <- data.frame(geneid = names(tss), chr = seqnames(tss), s1 = start(tss), s2 = end(tss))
    row.names(genepos) <- genepos$geneid
    logging::loginfo(paste(length(unique(tss_within_dist$subjectHits)),"genes within SNPs"))
    }


    logging::loginfo("Finding SNP CpG overlaps")
    snp_cpg_combo <- findOverlaps(query = snp_loci, subject = cpg_loci)

    sample_info_mat_high <- matrix(NA, ncol = 64, nrow = length(unique(to(snp_cpg_combo))))
    sample_info_mat_low <- matrix(NA, ncol = 64, nrow = length(unique(to(snp_cpg_combo))))
    assoc_snps <- list()
    logging::loginfo("Preparing high and low datagroups")
    for (i in 1:length(unique(to(snp_cpg_combo)))){
        j <- unique(to(snp_cpg_combo))[i]
        assoc_snps[[i]] <- names(snp_loci[from(snp_cpg_combo[to(snp_cpg_combo) == j])])
        m <- as.numeric(t[names(cpg_loci[j]), ])
        ids <- colnames(t[names(cpg_loci[j]), m<median(m)])
        length(ids) <- 64
        sample_info_mat_low[i,] <- ids
        ids <- colnames(t[names(cpg_loci[j]), m>median(m)])
        length(ids) <- 64
        sample_info_mat_high[i,] <- ids
    }
    row.names(sample_info_mat_low) <- names(cpg_loci[unique(to(snp_cpg_combo))])
    row.names(sample_info_mat_high) <- names(cpg_loci[unique(to(snp_cpg_combo))])
    names(assoc_snps) <- names(cpg_loci[unique(to(snp_cpg_combo))])

    logging::loginfo("Running eQTLs in high group")
    eqtls_high <- list()
    length(eqtls_high) <- length(assoc_snps)

    cl <- parallel::makeCluster(16)
    doParallel::registerDoParallel(cores = 6)

    #eqtls_high <- foreach::foreach (i = 1:length(assoc_snps), .combine = rbind, .packages = "MatrixEQTL") %dopar%

    eqtls_high <- foreach::foreach (i = 1:16, .combine = rbind, .packages = "MatrixEQTL") %dopar% {
    #eqtls_high <- rbind(eqtls_high,foreach (i = 1:1000, .combine = rbind, .packages = "MatrixEQTL") %dopar% {
    #eqtls_high <- foreach (i = 1:1000, .combine = rbind, .packages = "MatrixEQTL") %dopar% {
        logging::loginfo(i)
        message(i)
        j <- names(assoc_snps[i])
        snps = SlicedData$new()
        snps$initialize(as.matrix(gdata[assoc_snps[[i]], na.omit(as.character(sample_info_mat_high[j,]))]))
        gene_mat <- edata[unique(unlist(snps_genes_assoc[snps$GetAllRowNames()])),na.omit(as.character(sample_info_mat_high[j,]))]
        gene_mat <- t(apply(gene_mat, 1, RNOmni::rankNorm))
        if(assertthat::not_empty(gene_mat)){
            gene = SlicedData$new()
            gene$initialize(gene_mat)
            meh <- Matrix_eQTL_main(
                snps = snps,
                gene = gene,
                snpspos = snppos,
                genepos = genepos,
                cvrt = SlicedData$new(),
                output_file_name = "",
                output_file_name.cis  = "",
                pvOutputThreshold = 0.05,
                pvOutputThreshold.cis = 1,
                useModel = modelLINEAR,
                errorCovariance = numeric(),
                verbose = TRUE,
                pvalue.hist = F,
                cisDist = 500000)
            meh$cis$eqtls$cpg <- j
            return(meh$cis$eqtls)
        }
    }

    save(eqtls_high, file = "data/eqtls_high.Rdata")
    parallel::stopCluster(cl)

    logging::loginfo("Running eQTLs in low group")
    eqtls_low <- foreach::foreach (i = 1:length(assoc_snps), .combine = rbind, .packages = "MatrixEQTL") %dopar% {
        logging::loginfo(i)
        j <- names(assoc_snps[i])
        snps = SlicedData$new()
        snps$initialize(as.matrix(gdata[assoc_snps[[i]], na.omit(as.character(sample_info_mat_low[j,]))]))
        gene_mat <- edata[unique(unlist(snps_genes_assoc[snps$GetAllRowNames()])),na.omit(as.character(sample_info_mat_low[j,]))]
        gene_mat <- t(apply(gene_mat, 1, RNOmni::rankNorm))
        if(assertthat::not_empty(gene_mat)){
            gene = SlicedData$new()
            gene$initialize(gene_mat)
            meh <- Matrix_eQTL_main(
                snps = snps,
                gene = gene,
                snpspos = snppos,
                genepos = genepos,
                cvrt = SlicedData$new(),
                output_file_name = "",
                output_file_name.cis  = "",
                pvOutputThreshold = 0.05,
                pvOutputThreshold.cis = 1,
                useModel = modelLINEAR,
                errorCovariance = numeric(),
                verbose = TRUE,
                pvalue.hist = F,
                cisDist = 500000)
            meh$cis$eqtls$cpg <- j
            return(meh$cis$eqtls)
        }
    }
    save(eqtls_low, file = "data/eqtls_low.Rdata")
    parallel::stopCluster(cl)


    logging::loginfo("eQTLs done. Preparing result file")
    th <- eqtls_high %>% arrange(cpg, snps,gene)
    tl <- eqtls_low %>% arrange(cpg, snps,gene)
    assertthat::assert_that(assertthat::are_equal(th$snps, tl$snps))
    assertthat::assert_that(assertthat::are_equal(th$gene, tl$gene))
    assertthat::assert_that(assertthat::are_equal(th$cpg, tl$cpg))
    colnames(th) <- paste0("high_", colnames(th))
    colnames(tl) <- paste0("low_", colnames(tl))
    t <- cbind(th[,c(1,2,7,3,4,5,6)],tl[,3:6])
    colnames(t)[1:3] <- c("snp","gene","cpg")
    t$ctcf <- plyr::mapvalues(t$cpg, cpg_ctcf$cpg, cpg_ctcf$ctcf, warn_missing = F)
    t$atac <- plyr::mapvalues(t$snp, snp_atac$V12, snp_atac$V4, warn_missing = F)
    t$comb <- paste0(t$atac,":",t$ctcf,":",t$gene)
    matrixeqtl_res <- t
    matrixeqtl_res$siglow <- ifelse(matrixeqtl_res$low_pvalue<=0.05 & matrixeqtl_res$low_FDR<=0.1, TRUE, FALSE)
    matrixeqtl_res$sighigh <- ifelse(matrixeqtl_res$high_pvalue<=0.05 & matrixeqtl_res$high_FDR<=0.1, TRUE, FALSE)
    matrixeqtl_res$sigdelta <- ifelse((matrixeqtl_res$sighigh & !matrixeqtl_res$siglow) | (matrixeqtl_res$siglow & !matrixeqtl_res$sighigh), TRUE, FALSE)
    matrixeqtl_res$classSimplified <- "None"
    matrixeqtl_res[matrixeqtl_res$sighigh & !matrixeqtl_res$siglow, 'classSimplified'] <- "Insulating"
    matrixeqtl_res[!matrixeqtl_res$sighigh & matrixeqtl_res$siglow, 'classSimplified'] <- "Activating"
    matrixeqtl_res$gene <- as.character(matrixeqtl_res$gene)
    matrixeqtl_res$snp <- as.character(matrixeqtl_res$snp)
    matrixeqtl_res$estimatedelta <- matrixeqtl_res$high_beta - matrixeqtl_res$low_beta
    row.names(snppos) <- snppos$snp
    matrixeqtl_res$seqnames <- snppos[matrixeqtl_res$snp, 'chr']
    matrixeqtl_res$snpPos <- snppos[matrixeqtl_res$snp, 'pos']
    matrixeqtl_res$tss <- genepos[matrixeqtl_res$gene, 's1']
    matrixeqtl_res$cpgPos <- unlist(lapply(matrixeqtl_res$cpg, function(x) strsplit(x, ":")[[1]][2]))
    matrixeqtl_res$end <- pmax(matrixeqtl_res$snpPos, matrixeqtl_res$tss)
    matrixeqtl_res$start <- pmin(matrixeqtl_res$snpPos, matrixeqtl_res$tss)
    #save(matrixeqtl_res, file = "data/matrixeqtl_res.Rdata")
    return(matrixeqtl_res)
}
