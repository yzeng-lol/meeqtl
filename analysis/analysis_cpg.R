plot(res)
plot_plot(res$cpg, class_col = "class", y_col = "encodeCorEst", method = "boxplot", las = 2, outline = F)
plot_plot(res, class_col = "class", y_col = "pdelta", nlog2 = T, method = "ecdf")


# cpg characterisation ----------------------------------------------------

plot_plot(res$cpg, class_col = "classSimplified", y_col = "probeDistFromTSS", method = "density", bw = 5000)
plot_plot(res$cpg, class_col = "classSimplified", y_col = "probeDistFromSNP", method = "density", bw = 1)

#prostate ctcf variability
ctcf_prostate <- find_ctcf_in_prostate(res$cpg, ctcf_metadata = file_ctcf_prostate_metadata, ctcf_dir = path_ctcf_prostate_beds)
assertthat::assert_that(all(row.names(ctcf_prostate) == res$cpg$probe))
plot(ecdf(rowSums(ctcf_prostate)[rowSums(ctcf_prostate)>0 & res$cpg$class=="Insulating_sigdelta"]), verticals = T, cex = 0,
     xlab = "No. of prostate samples", ylab = "Cumulative distribution", cex.axis = 0.8, cex.lab = 1, main = "CTCF Binding")
lines(ecdf(rowSums(ctcf_prostate)[rowSums(ctcf_prostate)>0 & res$cpg$class=="Activating_sigdelta"]), verticals = T, cex = 0, col = 2)
legend("bottomright", cex = 0.8, bty = "n", lwd = 1, col=1:2, legend = c("Insulating","Activating"))

#histone marks
histone_prostate <- get_histone_signals(res$cpg)
plot_histone_density(res$cpg, histone_prostate, histone = "H3K27me3", cell = "LNCAP", xlim=c(0,0.4))
plot_histone_density(res$cpg, histone_prostate, histone = "H3K27me3", cell = "VCaP", xlim=c(0,0.4))
plot_histone_density(res$cpg, histone_prostate, histone = "H3K27me3", cell = "RWPE-1", xlim=c(0,0.4))
plot_histone_density(res$cpg, histone_prostate, histone = "H3K27ac", cell = "LNCAP", xlim=c(0,0.4))
plot_histone_density(res$cpg, histone_prostate, histone = "H3K27ac", cell = "VCaP", xlim=c(0,0.4))
plot_histone_density(res$cpg, histone_prostate, histone = "H3K27ac", cell = "RWPE-1")
plot_histone_density(res$cpg, histone_prostate, histone = "H3K27me3", cell = "RWPE-1")


# blood -------------------------------------------------------------------


blood_prostate_cor <- correlate_blood_prostate(file_blood_catalog, file_prostate_catalog, mdata)
save(blood_prostate_cor, file = "data/blood_prostate_cor.Rdata")
load("data/blood_prostate_cor.Rdata")
res$cpg <- add_blood_cols(res$cpg, blood_prostate_cor)

# check activating cases --------------------------------------------------
tophit_act <- res$cpg %>% arrange(pdelta) %>% filter(classSimplified=="Activating", sigdelta)
i=2
plot_eqtl(res$ori, snp = tophit_act$snp[i], probe = tophit_act$probe[i], gene = tophit_act$gene[i], plot_what = "all")
#does activating cases have more 0 methylation?
t <- res$cpg %>% filter(classSimplified=="Activating", sigdelta) %>% pull(probe)
plot(density(rowMeans(mdata[t,])), main = "Chinese prostate data", xlab = "Mean methylation")
t <- res$cpg %>% filter(classSimplified=="Insulating", sigdelta) %>% pull(probe)
lines(density(rowMeans(mdata[t,])), col = 2)
legend("topleft", bty = "n", lty = 1, col=1:2, legend =c("Activating", "Insulating"), cex = 0.8)

