# this is the main script -------------------------------------------------
check_all_gsea(res$gene)
check_all_gsea(matrixeqtl_res)
gsea_kegg <- run_gsea(matrixeqtl_res, "~/Documents/DB/c2.cp.kegg.v7.1.symbols.gmt", plot = T,
                      pdffile = "plots/gsea_kegg_barplot.matrixeqtl.pdf")
gsea_pos <- run_gsea(matrixeqtl_res, "~/Documents/DB/c1.all.v7.1.symbols.gmt", plot = T,
                     pdffile = "plots/gsea_pos_barplot.matrixeqtl.pdf")
siggenes <- unlist(strsplit(gsea_pos$geneID[1], "\\/"))
bggenes <- unique(matrixeqtl_res %>% filter(sigdelta) %>% pull(gene))
gsea_bp <- run_gsea(matrixeqtl_res, "~/Documents/DB/c5.bp.v6.2.symbols.gmt", siggenes = siggenes, bggenes = bggenes, plot = T,
                    pdffile = "plots/gsea_bp_sigchr6.pdf")
run_gsea(matrixeqtl_res, "~/Documents/DB/c5.bp.v6.2.symbols.gmt", plot = F, minGSSize = 20)
chemokine_genes <- strsplit(gsea_kegg$geneID[1], "\\/")[[1]]
matrixeqtl_res
hotprobe <- names(sort(table(matrixeqtl_res$gene[chemokine_genes,'probe']), decreasing = T)[1])
res$cpg[hotprobe,]
View(res$gene[chemokine_genes,])
plot_eqtl(res$ori, snpid = "rs4694176_22740", probeid = hotprobe, geneid = "CXCL1", plot_what = "all")
plot_eqtl(res$ori, snpid = "rs4694176_22740", probeid = hotprobe, geneid = "CXCL5", plot_what = "all")
plot_eqtl(res$ori, snpid = "rs4694176_22740", probeid = hotprobe, geneid = "CXCL3", plot_what = "all")

find_reverse_eqtl(res$ori, file_edata = file_edata,
                file_gdata = file_gdata, file_mdata = file_mdata,
                probeid = hotprobe, snpid = "rs4694176_22740", methylation_status = "High")


# checking reverse eqtl for activating probes -----------------------------

t <- res$cpg %>% filter(classSimplified!="None")
reverse_eqtl_list_low <- list()
reverse_eqtl_list_high <- list()
for (i in 1:nrow(t)){
    logging::loginfo(i)
    reverse_eqtl_list_low[[i]] <- find_reverse_eqtl(res$ori, file_edata = file_edata,
                                              file_gdata = file_gdata, file_mdata = file_mdata,
                                              probeid = t$probe[i], snpid = t$snp[i],
                                              methylation_status = "Low", writetable = F)
    reverse_eqtl_list_high[[i]] <- find_reverse_eqtl(res$ori, file_edata = file_edata,
                                                  file_gdata = file_gdata, file_mdata = file_mdata,
                                                  probeid = t$probe[i], snpid = t$snp[i],
                                                  methylation_status = "High", writetable = F)
}
barplot(unlist(lapply(reverse_eqtl_list_low, function(x) length(x$value[x$value > -log2(0.05)]))))
barplot(unlist(lapply(reverse_eqtl_list_high, function(x) length(x$value[x$value > -log2(0.05)]))))
nlow <- unlist(lapply(reverse_eqtl_list_low, function(x) length(x$value[x$value > -log2(0.05)])))
nhigh <- unlist(lapply(reverse_eqtl_list_high, function(x) length(x$value[x$value > -log2(0.05)])))
res$gene$rsid <- unlist(lapply(res$gene$snp, function(x) paste(rev(rev(strsplit(x,"_")[[1]])[-1]), collapse = ":")))
for (i in 1:nrow(t)){
    genelist <- res$gene %>% filter(snp == t$snp[[i]],
                                    probe == t$probe[[i]]) %>% pull(gene)
    hgenes <- NA
    lgenes <- NA
    hgenes <- reverse_eqtl_list_high[[i]] %>% filter(value > -log2(0.05)) %>% pull(targetName)
    hgenes <- if(length(hgenes)>0) setdiff(hgenes, genelist)
    lgenes <- reverse_eqtl_list_low[[i]] %>% filter(value > -log2(0.05)) %>% pull(targetName)
    lgenes <- if(length(lgenes)>0) setdiff(lgenes, genelist)
    nhigh[[i]] <- ifelse(length(hgenes)>0, length(setdiff(hgenes, lgenes)), 0)
    nlow[[i]] <- ifelse(length(lgenes)>0, length(setdiff(lgenes, hgenes)), 0)
}
plot(seq(1,sum(t$classSimplified=="Activating")),
     as.numeric(nhigh[t$classSimplified=="Activating"]), type = "h",
     ylab = "eGenes", ylim=c(-10,10), yaxt = "n", frame = F, col = "coral2", xlab = "CpGs with activating function",
     cex.lab = 1, cex.axis = 0.8, cex.main = 1, main = "eGenes CpGs with activating function")
lines(seq(1,sum(t$classSimplified=="Activating")),
     -as.numeric(nlow[t$classSimplified=="Activating"]), type = "h",
     col = "lightblue4")
axis(2, at = c(-10,0,10), labels = c(10,0,10), cex = 0.8)
text(0,-10, labels = "eGene in low methylation", cex = 0.8, adj = 0, col = "lightblue4")
text(0,10, labels = "eGene in high methylation", cex = 0.8, adj = 0, col = "coral2")

plot(seq(1,sum(t$classSimplified=="Insulating")),
     as.numeric(nhigh[t$classSimplified=="Insulating"]), type = "h",
     ylab = "eGenes", ylim=c(-10,10), yaxt = "n", frame = F, col = "coral2", xlab = "CpGs with insulating function",
     cex.lab = 1, cex.axis = 0.8, cex.main = 1, main = "eGenes CpGs with insulating function")
lines(seq(1,sum(t$classSimplified=="Insulating")),
      -as.numeric(nlow[t$classSimplified=="Insulating"]), type = "h",
      col = "lightblue4")
axis(2, at = c(-10,0,10), labels = c(10,0,10), cex = 0.8)
text(0,-10, labels = "eGene in low methylation", cex = 0.8, adj = 0, col = "lightblue4")
text(0,10, labels = "eGene in high methylation", cex = 0.8, adj = 0, col = "coral2")

