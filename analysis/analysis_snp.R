# gwas --------------------------------------------------------------------

write.table(table(res$snp$gwasTrait[!is.na(res$snp$gwasTrait)]), "plots/gwas_traits.txt",
            sep = "\t", quote = F, row.names = F)
pdf("plots/gwas_count.pdf")
barplot(sort(table(
    res$snp %>% filter(sigdelta, classSimplified!="None", !is.na(gwasTrait)) %>% pull(gwasTrait))),
    las = 2, horiz = T, cex.axis = 0.8, cex.lab = 0.8,
        xlab = "GWAS count")
dev.off()
plot_gwas(res$snp, "plots/gwas_rank.pdf")

t <- res$snp %>% filter(grepl(pattern = "prostate", gwasTrait, ignore.case = T),
                        sigdelta, classSimplified!="None")
plot_eqtl(res$ori, snpid = t$snp[1], probeid = t$probe[1], geneid = t$gene[1], plot_what = "all")
plot_eqtl(res$ori, snpid = "rs1512270_42550", probeid = t$probe[1], geneid = t$gene[1], plot_what = "all")
find_reverse_eqtl(x = res$ori, file_edata = file_edata, file_gdata = file_gdata, file_mdata = file_mdata, probeid = t$probe[1], snpid = t$snp[1], methylation_status = "High")
find_reverse_eqtl(x = res$ori, file_edata = file_edata, file_gdata = file_gdata, file_mdata = file_mdata, probeid = t$probe[1], snpid = t$snp[1], methylation_status = "Low")
find_reverse_eqtl(x = res$ori, file_edata = file_edata, file_gdata = file_gdata, file_mdata = file_mdata, probeid = t$probe[1], snpid = "rs1512270_42550", methylation_status = "High", cutAtMedian = T)
find_reverse_eqtl(x = res$ori, file_edata = file_edata, file_gdata = file_gdata, file_mdata = file_mdata, probeid = t$probe[1], snpid = "rs1512270_42550", methylation_status = "Low", cutAtMedian = T)

snpid <- "rs1512270_42550"
probeid <- t$probe[1]
