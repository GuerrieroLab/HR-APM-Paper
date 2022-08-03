# Load packages
library(RColorBrewer)
library(dendextend)
library(gplots)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(rstatix)
library(gtools)

#Load in CCLE data and 163 ERS, APM/TC panel genes
setwd("/Users/kenichishimada/Dropbox (HMS)/_projects_/B06 03 HR-APM/ERS_APM Manuscript/CCLE Data analysis/rda")
load("ccle_RNA-seq.rda")
load("common_signatures.rda")
comm.gns <- c(comm.esr,comm.apm.tc)
cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
hc <- function(x)hclust(x,method="complete")

#Filter for 163 ERS an APM/TC panel genes
ccle.exp <- log2(ccle.m1.breast[match(comm.gns,ccle.syms),]+1)
rownames(ccle.exp) <- comm.gns
sds <- apply(ccle.exp,1,sd)
which(sds==0) ## remove these two genes, IFI30, IFNG
# ccle.exp <- ccle.exp[sds>0,]
#sapply(ccle.exp,function(x)sum(is.na(x))) # no missing values

#set col colors
side.col.annotation<- sub.ccle %>% select("Cell.line.name","PAM50.mRNA") %>% arrange(PAM50.mRNA,Cell.line.name)
unique(side.col.annotation$PAM50.mRNA)

# order of annotations/colors are defined here
ccle.exp.graph <- ccle.exp[,match(side.col.annotation$Cell.line.name,colnames(ccle.exp))]

groups<- side.col.annotation$PAM50.mRNA
coloursSamples <- factor(groups, levels=c("Basal-like","Her2amp","Luminal A","Luminal B"))
coloursSamples <- colorRampPalette(c("navyblue","hotpink","purple","green"))(length(unique(coloursSamples)))[factor(coloursSamples)]

lcol <- brewer.pal(3,"Set1")[(comm.gns %in% comm.esr)+1]
#lcol <- lcol[!comm.gns %in% c("IFI30","IFNG")]

#Plot heatmap for CCLE PAM50
setwd("/Users/kenichishimada/Dropbox (HMS)/_projects_/B06 03 HR-APM/ERS_APM Manuscript/CCLE Data analysis/Graphs")
pdf(file="CCLE_heatmap.pdf")
heatmap_ccle<- heatmap.2(ccle.exp.graph, trace="none",dendrogram='none',Rowv=FALSE,Colv=FALSE, 
                         col=cols[],
                         RowSideColors=lcol,ColSideColors= coloursSamples,
                         hclustfun=hc,cexRow = 0.2,cexCol = 0.2, scale = "row")
legend("topright",c("Basal-like","Her2amp","Luminal A","Luminal B"),fill =c("navyblue","hotpink","purple","green"),border=NA, box.lwd = 0) 
dev.off()

##total 161 genes, removed IFI30, and IFNG for cor matrix as not expression was detected.
ccle.exp.graph.cor<- ccle.exp.graph[!comm.gns %in% c("IFI30","IFNG"),]
ccle_cormat<- cor(t(ccle.exp.graph.cor), use ="everything",method= "spearman")
#set ERS, APM/TC gene colors
lcol.cor <- brewer.pal(3,"Set1")[(colnames(ccle_cormat) %in% comm.esr)+1]

setwd("/Users/kenichishimada/Dropbox (HMS)/_projects_/B06 03 HR-APM/ERS_APM Manuscript/CCLE Data analysis/Graphs")
pdf(file="ccle_scmatrix_2021.08.02.pdf")
heatmap.2(ccle_cormat, trace="none",#dendrogram='none',Rowv=FALSE,Colv=FALSE, 
          RowSideColors=lcol.cor,ColSideColors=lcol.cor,hclustfun=hc,
          col=cols,cexRow = 0.2,cexCol = 0.2) 
dev.off()

## k-means (via GMM)
library(mclust)
bics <- mclustBIC(ccle.exp,G=1:50,modelNames=c("EII"))
plot(bics) # G=4 looks acceptable
m <- Mclust(ccle.exp,G=4,modelNames="EII",x=bics)

cl.km <- factor(m$classification) # 161 genes
cl.km <- cl.km[!comm.gns %in% c("IFI30","IFNG")]

gs <- rep(c("ESR","TC"),c(length(comm.esr),length(comm.apm.tc)))
gs <- gs[!comm.gns %in% c("IFI30","IFNG")]

table(cluster=cl.km,genesig=gs)

length(cl.km)

cl.hc <- factor(cutree(hc1$colDendrogram,k=4))
lcol2.cor <- brewer.pal(6,"Set1")[as.numeric(cl.hc)] ## cl.hc

hc1 <- heatmap.2(ccle_cormat, trace="none",#dendrogram='none',Rowv=FALSE,Colv=FALSE, 
          RowSideColors=lcol2.cor,ColSideColors=lcol2.cor,hclustfun=hc,
          col=cols,cexRow = 0.2,cexCol = 0.2) 

gns.km <- tapply(names(cl.km),cl.km,identity)
gns.hc <- tapply(names(cl.hc),cl.hc,identity)

table(cl.km,cl.hc)

diag(ccle_cormat) <- NA

cl <- cl.km
mean.cor <- sapply(levels(cl),function(x){
  sapply(levels(cl),function(y){
    mean(ccle_cormat[cl==x,cl==y],na.rm=T)
  })
})

slp.cor <- sapply(levels(cl),function(x){
  sapply(levels(cl),function(y){
    mean(ccle_cormat[cl==x,cl==y],na.rm=T)
  })
})

uniq.lcol <- brewer.pal(4,"Set1")
steps <- seq(-1,1,length=50)
d <- 50-which.min(abs(max(abs(mean.cor))-steps))
sub.cols <- cols[(d+1):(50-d)]

heatmap.2(mean.cor, trace="none",#dendrogram='none',Rowv=FALSE,Colv=FALSE, 
          RowSideColors=uniq.lcol,ColSideColors=uniq.lcol,hclustfun=hc,
          col=sub.cols,cexRow = 0.2,cexCol = 0.2) 
