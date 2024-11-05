## colors
## cell types - CyCIF
uniq.cols <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(13)
names(uniq.cols) <- 1:13

uniq.cts <- c("Cancer","Fibro","Endo","CD8T","CD4T_nonTreg","CD4T_Treg","T_other","NK","B",
              "Mac_CD68","Mac_CD68_CD163","Mac_CD163","Imm_other","all_other","outOfROI")　# 'unknown' is removed

col.idx <- c(CD8T=1,NK=2,CD4T_nonTreg=3,T_other=7,Endo=6,Imm_other=5,Fibro=4,
             Cancer=8,Mac_CD163=9,Mac_CD68_CD163=10,Mac_CD68=11,CD4T_Treg=12,B=13)

ct.cols <- c(uniq.cols[col.idx[uniq.cts[1:13]]],"grey70")
names(ct.cols) <- uniq.cts[1:14]


## WEHI - scRNA-seq
used.cts <- c("Cancer","Fibroblast","Pericyte","Endothelial","Mono_Mac","Dendritic","Mast","T","NK","B","Plasma")
exp.ct.cols1 <- RColorBrewer::brewer.pal(11,"Spectral")[c(7,1,6,9,3,8,10,4,2,11,5)]
names(exp.ct.cols1) <- used.cts
if(0){
  plot(1:11,pch=19,col=exp.ct.cols1)
  text(1:11,1:11,paste0(1:11,names(exp.ct.cols1)))
}

st.uniq.cols <- c("#f3756c","#7caf42","#16bcc2")
names(st.uniq.cols) <- c("HR+","HER2+","TNBC")

### module colors
mod.uniq.cols <- c("steelblue",RColorBrewer::brewer.pal(3,"Accent")[c(1,3,2)])#f3756c","#7caf42","#16bcc2","#a780ba")
names(mod.uniq.cols) <- c("ERS","TNFa/NFkB","IFN-I","APM/TC")

## cell types - snRNA-seq
exp.ct.cols2 <- c(exp.ct.cols1[1:4],Adipo="grey50",exp.ct.cols1[5:11],Other="grey70")
names(exp.ct.cols2) <- c("Cancer","Fibro","Peri","Endo","Adipo","Mono_Mac","DC","Mast","T","NK","B","Plasma","Other")




