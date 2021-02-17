#' This will run matrixeQTL for eQTL, meQTL and methylation-dependent eQTL
#' @name run_me_eQTL
#' @param x eqtl dataframe, e.g., res$gene
#' @param file_edata path to edata file
#' @param file_gdata path to gdata file
#' @param file_mdata path to mdata file
#' @param file_tss_loci path to TSS bed file
#' @param file_snp_loci path to SNP bed file
#' @param file_cpg_loci path to CpG bed file
#' @param geneSnpMaxDistance max distance between snp and gene, default 500000
#' @param cpgSnpMaxDistance max distance between snp and probe, default 500000
#' @param cluster_core doParallel core, default 2 for local
#' @import GenomicRanges
#' @import GenomicFeatures
#' @import bigstatsr
#' @import assertthat
#' @import MatrixEQTL
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom rtracklayer import
#' @importFrom ChIPpeakAnno reCenterPeaks
#' @importFrom logging loginfo
#' @export


run_me_eQTL <- function(file_edata, file_gdata, file_mdata, file_tss_loci, file_snp_loci,file_cpg_loci,
                           geneSnpMaxDistance = 500000, cpgSnpMaxDistance = 500000, cluster_core = 2)
{
###########################################
## combination for snp-CpG-gene for testing
###########################################
    {

        #######################################################
        ## read in expression, genotype, methylation matrix ...
        {
            logging::loginfo("Reading in files ...")
            edata <- read.table(file_edata, header = T, as.is = T)
            gdata <- read.table(file_gdata, header = T, as.is = T)
            mdata <- read.table(file_mdata, header = T, as.is = T)

            ## loci bed files without header
            tss_loci <- read.table(file_tss_loci, header = F, as.is = T)
            snp_loci <- read.table(file_snp_loci, header = F, as.is = T)
            cpg_loci <- read.table(file_cpg_loci, header = F, as.is = T)

            ## checking whether input files meet requirements and consistent
            ## ...

        }

        #######################################################
        ##  CpG sites
        {
            logging::loginfo(paste("Number of CpGs:", nrow(mdata)))

            ## Extenting to CpG site centered "cpgSnpMaxDistance * 2" region
            cpg_loci_ext <- GRanges(seqnames = cpg_loci$V1,
                                    ranges = IRanges(start = cpg_loci$V5 - cpgSnpMaxDistance * 2,
                                                     end = cpg_loci$V5 + cpgSnpMaxDistance * 2,
                                                     names = cpg_loci$V4))

            start(cpg_loci_ext)[start(cpg_loci_ext) < 0] <- 0

            cpg4qtl <- data.frame(cpg = cpg_loci$V4, chr = cpg_loci$V1, start = cpg_loci$V2, start = cpg_loci$V3)

        }

        #######################################################
        ##  SNPs
        {
            logging::loginfo(paste("Number of SNPs:", nrow(gdata)))
            snp_loci_gr <- GRanges(snp_loci$V1, IRanges(snp_loci$V2, snp_loci$V3, names = snp_loci$V4))

            ## intersect SNPs with extented CpG region
            snp_in_cpg_ext <- subsetByOverlaps(snp_loci_gr, cpg_loci_ext)
            logging::loginfo(paste("Number of SNPs within extented CpG regions:", length(snp_in_cpg_ext)))

            ## snp_in_cpg_ext for matrieqtl
            snp4qtl <- data.frame(snp = names(snp_in_cpg_ext), chr = seqnames(snp_in_cpg_ext), pos = end(snp_in_cpg_ext))

            ## Extenting to SNP site centered "geneSnpMaxDistance * 2" region
            snp4qtl_ext <- ChIPpeakAnno::reCenterPeaks(snp_in_cpg_ext, geneSnpMaxDistance * 2)
            start(snp4qtl_ext)[start(snp4qtl_ext) < 0] <- 0
        }

        #######################################################
        ##  SNP and genes combination
        {
            logging::loginfo(paste("Number of Genes:", nrow(edata)))

            ## intersect Genes with extented SNP region
            tss_loci_gr <- GRanges(tss_loci$V1, IRanges(tss_loci$V2, tss_loci$V3, names = tss_loci$V4))

            ## snp list with all gene in extented region
            tss_in_snp_ext <- as.data.frame(findOverlaps(query = tss_loci_gr, subject = snp4qtl_ext))
            snps_genes_assoc <- split(tss_in_snp_ext$queryHits, as.character(tss_in_snp_ext$subjectHits))      ## genes in extended SNP region

            ## SNP rsIDs
            snp_idx <- as.numeric(names(snps_genes_assoc))
            names(snps_genes_assoc) <- names(snp4qtl_ext)[snp_idx]

            ## gene names
            snps_genes_assoc <- lapply(snps_genes_assoc, function(x) names(tss_loci_gr[x]))

            logging::loginfo(paste("Number of SNPs with gene in extented region:", length(snps_genes_assoc)))
            logging::loginfo(paste("Number of Genes with in SNP extedted region:", length(unique(tss_in_snp_ext$queryHits))))

            gene4qtl <- data.frame(gene_name = names(tss_loci_gr), chr = seqnames(tss_loci_gr), start = start( tss_loci_gr), end = end(tss_loci_gr))
            row.names(gene4qtl) <- gene4qtl$gene_name
        }


        #######################################################
        ##  SNP and CpG combination
        {
            snp_cpg_combo <- findOverlaps(query = snp_in_cpg_ext, subject = cpg_loci_ext)

            cpg_num <- length(unique(to(snp_cpg_combo)))     ## number of CpG sits with SNP in extended region
            cpg_idx <- to(snp_cpg_combo)

            sample_info_mat_high <- matrix(NA, ncol = ncol(mdata)/2, nrow = cpg_num)
            sample_info_mat_low <-  matrix(NA, ncol = ncol(mdata)/2, nrow = cpg_num)
            cpgs_snps_assoc <- list()

            logging::loginfo("Spliting samples into methlylation high and low groups")
            for (i in 1:cpg_num){

                j <- unique(cpg_idx)[i]
                cpgs_snps_assoc[[i]] <- names(snp_in_cpg_ext[from(snp_cpg_combo[cpg_idx == j])])
                names(cpgs_snps_assoc)[i] <- names(cpg_loci_ext)[j]

                ## methylation for a single CpG site
                cpg_m <- as.numeric(mdata[j, ])

                idx_ml <- cpg_m < median(cpg_m)
                s_ml <- colnames(mdata)[idx_ml]
                length(s_ml) <- ncol(mdata)/2             ## adding sample "NA" to make sure length is the same
                sample_info_mat_low[i,] <- s_ml

                idx_mh <- cpg_m > median(cpg_m)
                s_mh <- colnames(mdata)[idx_mh]
                length(s_mh) <- ncol(mdata)/2
                sample_info_mat_high[i,] <- s_mh
            }

            rownames(sample_info_mat_low) <- rownames(sample_info_mat_high) <- names(cpgs_snps_assoc)
            logging::loginfo(paste("Number of CpGs with SNPs in extented region:", length(cpgs_snps_assoc)))

        }


    }

###########################################
## running Matrixqtl
###########################################
    {

        ###########################################
        ##  Running cis eQTLs
        {
            logging::loginfo("Running eQTLs  ...")

            ## doParallel
            cl <- parallel::makeCluster(cluster_core)
            doParallel::registerDoParallel(cores = cluster_core)

            snps =  MatrixEQTL::SlicedData$new()
            idx_snp <- match(snp4qtl$snp, rownames(gdata))
            snps$initialize(as.matrix(gdata[idx_snp, ]))

            gene_mat <- t(apply(edata[gene4qtl$gene_name, ], 1, RNOmni::rankNorm))
            gene = SlicedData$new()
            gene$initialize(gene_mat)

            eqtl <- Matrix_eQTL_main(
                snps = snps,
                gene = gene,
                snpspos = snp4qtl,
                genepos = gene4qtl,
                cvrt = SlicedData$new(),
                output_file_name = "",
                output_file_name.cis  = "",
                pvOutputThreshold = 0,                  ## ignore transQTLs
                pvOutputThreshold.cis = 1,
                useModel = modelLINEAR,
                errorCovariance = numeric(),
                verbose = TRUE,
                pvalue.hist = F,
                cisDist = 500000)
            cis_eqtls <- eqtl$cis$eqtls
            save(cis_eqtls, file = "cis_eqtls.Rdata")
            parallel::stopCluster(cl)

        }

        ###########################################
        ##  Running cis meQTLs
        {
            logging::loginfo("Running meQTLs  ...")

            ## doParallel
            cl <- parallel::makeCluster(cluster_core)
            doParallel::registerDoParallel(cores = cluster_core)

            snps = MatrixEQTL::SlicedData$new()
            idx_snp <- match(snp4qtl$snp, rownames(gdata))
            snps$initialize(as.matrix(gdata[idx_snp, ]))

            me_mat <- t(apply(mdata[cpg4qtl$cpg, ], 1, RNOmni::rankNorm))
            cpg_me = SlicedData$new()
            cpg_me$initialize(me_mat)

            meqtl <- Matrix_eQTL_main(
                snps = snps,
                gene = cpg_me,                          ## CpG methylation level
                snpspos = snp4qtl,
                genepos = cpg4qtl,                      ## CpG postion
                cvrt = SlicedData$new(),
                output_file_name = "",
                output_file_name.cis  = "",
                pvOutputThreshold = 0,                  ## ignore transQTLs
                pvOutputThreshold.cis = 1,
                useModel = modelLINEAR,
                errorCovariance = numeric(),
                verbose = TRUE,
                pvalue.hist = F,
                cisDist = 500000)
            cis_meqtls <- meqtl$cis$eqtls
            colnames(cis_meqtls)[2] <- "cpg"              ## replacing "gene" to "cpg"
            save(cis_meqtls, file = "cis_meqtls.Rdata")
            parallel::stopCluster(cl)

        }

        ###########################################
        ##  Running eQTLs in methylation high group
        {
        logging::loginfo("Running eQTLs in methylation high group ...")

        ## doParallel
        cl <- parallel::makeCluster(cluster_core)
        doParallel::registerDoParallel(cores = cluster_core)

        eqtls_me_high <- foreach::foreach (i = 1:length(cpgs_snps_assoc), .combine = rbind, .packages = "MatrixEQTL") %dopar%
                      {
                            ## Testing processing
                            if (i %% 10 == 0) {
                                logging::loginfo(paste("Run tests: ", i ))
                                message(paste("Run tests: ", i ))
                            }

                            j <- names(cpgs_snps_assoc[i])         ## testing CpG

                            ## genotype of snps
                            snps =  MatrixEQTL::SlicedData$new()
                            snps$initialize(as.matrix(gdata[cpgs_snps_assoc[[i]], na.omit(as.character(sample_info_mat_high[j,]))]))

                            ## gene expression for correspoding snps
                            gene_mat <- edata[unique(unlist(snps_genes_assoc[snps$GetAllRowNames()])),na.omit(as.character(sample_info_mat_high[j,]))]
                            gene_mat <- t(apply(gene_mat, 1, RNOmni::rankNorm))

                            ## run MatrixQTL
                            if(assertthat::not_empty(gene_mat)){
                                gene = SlicedData$new()
                                gene$initialize(gene_mat)
                                meh <- Matrix_eQTL_main(
                                                        snps = snps,
                                                        gene = gene,
                                                        snpspos = snp4qtl,
                                                        genepos = gene4qtl,
                                                        cvrt = SlicedData$new(),
                                                        output_file_name = "",
                                                        output_file_name.cis  = "",
                                                        pvOutputThreshold = 0,
                                                        pvOutputThreshold.cis = 1,     ## output all results to make sure me_high and me_low matched
                                                        useModel = modelLINEAR,
                                                        errorCovariance = numeric(),
                                                        verbose = TRUE,
                                                        pvalue.hist = F,
                                                        cisDist = 500000)
                                meh$cis$eqtls$cpg <- j
                                return(meh$cis$eqtls)
            }
        }

        save(eqtls_me_high, file = "eqtls_me_high.Rdata")
        parallel::stopCluster(cl)
        }

        ###########################################
        ##  Running eQTLs in methylation low group
        {
            logging::loginfo("Running eQTLs in methylation low group ...")

            ## doParallel
            cl <- parallel::makeCluster(cluster_core)
            doParallel::registerDoParallel(cores = cluster_core)

            eqtls_me_low <- foreach::foreach (i = 1:length(cpgs_snps_assoc), .combine = rbind, .packages = "MatrixEQTL") %dopar%
            {
                ## Testing processing
                if (i %% 10 == 0) {
                    logging::loginfo(paste("Run tests: ", i ))
                    message(paste("Run tests: ", i ))
                }

                j <- names(cpgs_snps_assoc[i])         ## testing CpG

                ## genotype of snps
                snps =  MatrixEQTL::SlicedData$new()
                snps$initialize(as.matrix(gdata[cpgs_snps_assoc[[i]], na.omit(as.character(sample_info_mat_low[j,]))]))

                ## gene expression for correspoding snps
                gene_mat <- edata[unique(unlist(snps_genes_assoc[snps$GetAllRowNames()])),na.omit(as.character(sample_info_mat_low[j,]))]
                gene_mat <- t(apply(gene_mat, 1, RNOmni::rankNorm))

                ## run MatrixQTL
                if(assertthat::not_empty(gene_mat)){
                    gene = SlicedData$new()
                    gene$initialize(gene_mat)
                    meh <- Matrix_eQTL_main(
                        snps = snps,
                        gene = gene,
                        snpspos = snp4qtl,
                        genepos = gene4qtl,
                        cvrt = SlicedData$new(),
                        output_file_name = "",
                        output_file_name.cis  = "",
                        pvOutputThreshold = 0,
                        pvOutputThreshold.cis = 1,       ## output all results to make sure me_high and me_low matched
                        useModel = modelLINEAR,
                        errorCovariance = numeric(),
                        verbose = TRUE,
                        pvalue.hist = F,
                        cisDist = 500000)
                    meh$cis$eqtls$cpg <- j
                    return(meh$cis$eqtls)
                }
            }

            save(eqtls_me_low, file = "eqtls_me_low.Rdata")
            parallel::stopCluster(cl)
        }
    }

######################################################
##  Combining eQTLs in methylation hihg and low group
#####################################################
    {
        logging::loginfo("Combining all QTLs results")

        ## merge methylation dependent eQTLs results
        ## arrange by CpG, SNP and Gene
        th <- eqtls_me_high %>% arrange(cpg, snps, gene)
        tl <- eqtls_me_low %>% arrange(cpg, snps, gene)
        colnames(th) <- paste0("me_high_eqtl_", colnames(th))
        colnames(tl) <- paste0("me_low_eqtl_", colnames(tl))

        ## check whether eqtl_me_hi and eqtl_me_low are in same order
        assertthat::assert_that(assertthat::are_equal(th$snps, tl$snps))
        assertthat::assert_that(assertthat::are_equal(th$gene, tl$gene))
        assertthat::assert_that(assertthat::are_equal(th$cpg, tl$cpg))

        ## combine th and tl
        res <- cbind(th[, c(1, 7, 2, 3:6)], tl[, 3:6])
        colnames(res)[1:3] <- c("snp","cpg", "gene")

        ## only remain me_eqtl either in me_hihg or me_low or both groups
        idx_rm = res$me_high_eqtl_FDR >= 0.05 & res$me_low_eqtl_FDR >= 0.05
        idx_rm[is.na(idx_rm)] <- TRUE
        res <- res[!idx_rm, ]

        ## adding eqtl information
        snp_gene <- paste(res$snp, res$gene, sep = "_")
        snp_gene_ex <- paste(cis_eqtls$snp, cis_eqtls$gene, sep = "_")
        idx_snp_gene <- match(snp_gene, snp_gene_ex)           ## cis_meqtl migth not avaliable

        match_eqtl <- cbind(snp_gene, cis_eqtls[idx_snp_gene, 3:6])
        colnames(match_eqtl) <- paste0("eqtl_", colnames(match_eqtl))

        ## adding  meqtl information
        snp_cpg <- paste(res$snp, res$cpg, sep = "_")
        snp_cpg_me <- paste(cis_meqtls$snp, cis_meqtls$cpg, sep = "_")
        idx_snp_cpg <- match(snp_cpg, snp_cpg_me)           ## cis_meqtl migth not avaliable

        match_meqtl <- cbind(snp_cpg, cis_meqtls[idx_snp_cpg, 3:6])
        colnames(match_meqtl) <- paste0("meqtl_", colnames(match_meqtl))

        ## merging all
        res <- cbind(res, match_eqtl, match_meqtl)
        rownames(res) <- paste(res$snp, res$cpg, res$gene, sep = "_")
        return(res)

    }

}
