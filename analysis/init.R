
# input file list ---------------------------------------------------------
file_eqtl_res <- "/Users/musaahmed/Google\ Drive/ctcf/out_all.noNA.normalized.txt"
file_encode_cor <- "/Users/musaahmed/Google\ Drive/ctcf/encode_correlation/methylation_correlated_with_ctcf.Padj_filtered.bed"
file_atac <- "/Users/musaahmed/Google\ Drive/ctcf/CPGEA_gdata.ATAC.txt"
file_gdata <- "/Users/musaahmed/Google\ Drive/ctcf/CPGEA_gdata.txt"
file_edata <- "/Users/musaahmed/Google\ Drive/ctcf/CPGEA_edata.normal.txt"
file_mdata <- "/Users/musaahmed/Google\ Drive/ctcf/me-eqtl/meeqtl/data/CPGEA_mdata.normal.txt"
file_tss <- "/Users/musaahmed/Google Drive/ctcf/me-eqtl/meeqtl/data/CPGEA_edata.normal.TSS.bed"
file_snp_loci <- "/Users/musaahmed/Google Drive/ctcf/me-eqtl/meeqtl/data/CPGEA_gdata.bed"
file_blood_catalog <- "/Users/musaahmed/Google\ Drive/ctcf/blood/catalog.N.sorted.merged.bed"
file_prostate_catalog <- "/Users/musaahmed/Google\ Drive/ctcf/encode_correlation/catalog.N.sorted.merged.all_samples.bed"
file_ctcf_prostate_metadata <- "/Volumes/Seagate_Backup_Plus/ENCODE/TF/CTCF/prostate/metadata.txt"
file_gwas <- "~/Google Drive/ctcf/CPGEA_gdata.GWAS_all.txt"
file_atac <- "/Users/musaahmed/Google Drive/ctcf/me-eqtl/meeqtl/data/CPGEA_gdata.ATAC.txt"
file_ctcf <- "/Users/musaahmed/Google Drive/ctcf/me-eqtl/meeqtl/data/methylation_correlated_with_ctcf.Padj_filtered.bed"
path_ctcf_prostate_beds <- "/Volumes/Seagate_Backup_Plus/ENCODE/TF/CTCF/prostate/"

gdata <- read.table(file_gdata, header=T, row.names = 1, sep="\t")
edata <- read.table(file_edata, header=T, row.names = 1, sep="\t")
mdata <- read.table(file_mdata, header=T, row.names = 1, sep="\t")


# matrixeqtl --------------------------------------------------------------

res <- run_matrixeqtl(file_edata = file_edata,
                    file_gdata = file_gdata,
                    file_mdata = file_mdata,
                    file_snp_loci = file_snp_loci,
                    file_tss = file_tss,
                    file_atac = file_atac,
                    file_ctcf = file_ctcf)
save(res, file = "data/matrixeqtl_res.Rdata")
load("data/matrixeqtl_res.Rdata")
length(unique(matrixeqtl_res %>% pull(cpg)))
length(unique(matrixeqtl_res %>% pull(snp)))
length(unique(matrixeqtl_res %>% pull(gene)))
length(unique(matrixeqtl_res %>% filter(sigdelta) %>% pull(cpg)))
length(unique(matrixeqtl_res %>% filter(sigdelta) %>% pull(snp)))
length(unique(matrixeqtl_res %>% filter(sigdelta) %>% pull(gene)))



# this is the main script -------------------------------------------------
res<-load_eqtl_res(file_eqtl_res, file_encode_cor, file_atac)
save(res,file = "data/res.Rdata")
load("data/res.Rdata")
res <- split_eqtl_results(res)
res$ori <- add_siglabel_class(res$ori)
res$snp <- add_siglabel_class(res$snp)
res$snp <- add_gwas(res$snp, file_gwas = file_gwas)
res$gene <- add_siglabel_class(res$gene)
res$cpg <- add_siglabel_class(res$cpg)
make_tracks(res$ori)

res$cpg$probeDistFromTSS <- res$cpg$probePos - res$cpg$tss
res$cpg$probeDistFromSNP <- res$cpg$snpPos - res$cpg$probePos
