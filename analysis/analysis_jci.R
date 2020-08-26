jci_links <- get_jci_links(res$snp, "/Volumes/Seagate_Backup_Plus/ENCODE/chia-pet/134260-JCI-RG-RV-3_sd_416940.xlsx")
t <- res$snp %>% filter(sigdelta, classSimplified=="Activating")
t2 <- jci_links[res$snp$sigdelta & res$snp$classSimplified=="Activating",]
t <- res$snp %>% filter(sigdelta, classSimplified=="Insulating")
t2 <- jci_links[res$snp$sigdelta & res$snp$classSimplified=="Insulating",]
t2$snp <- row.names(t2)
table(is.na(t2$RWPE.1))
table(is.na(t2$LNCaP))
table(is.na(t2$VCaP))
length(unique(t2[!is.na(t2$RWPE.1) | !is.na(t2$LNCaP) |
                       !is.na(t2$VCaP) | !is.na(t2$DU145), 'snp']))
row.names(t) <- t$snp
all(row.names(t) == row.names(t2))
cells <- colnames(t2)[6:9]
flag <- c()
for (i in 1:nrow(t)){
    flag[i] <- ""
    for (cell in cells){
        if (t$gene[i] %in% unlist(strsplit(t2[[cell]][i],";"))){
            flag[i] <- paste0(flag[i],"+",cell)
        } else if (!is.na(t2$RWPE.1[i])){
            flag[i] <- paste0(flag[i],"+other")
        } else {
            flag[i] <- paste0(flag[i],"+NA")
        }
    }
}

