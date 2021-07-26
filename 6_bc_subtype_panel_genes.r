library(RColorBrewer)
library(gplots)
library(ggplot2)

setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/RData")
x <- load(file="gene_list.rda") # gene.sig, category, hallmarks
y <- load(file="tcga_RNA-seq.rda") # tcga.ids,tcga.clinical,tcga.syms,tcga.m1,tcga.m3
z <- load(file="metabric_subset.rda") # meta.ids,meta.clinical,meta.syms,meta.m1,meta.m3

setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/RData")
load(file="common_signatures.rda") # comm.apm.tc,comm.esr,comm.esr2

comm.gns <- c(comm.esr,comm.apm.tc)
lcol <- brewer.pal(3,"Set1")[(comm.gns %in% comm.esr)+1]

## sutypes
smpl.class <- names(tcga.ids)
hc <- function(x)hclust(x,method="complete") 

## TCGA - subtypes
tcga.exp <- tcga.m3[match(comm.gns,tcga.syms),] # 778 genes, 1093 samples
rownames(tcga.exp) <- comm.gns

tcga.sp <- list()
for(cl in smpl.class){
	used.ids <- tcga.ids[[cl]]
	tcga.exp.sub <- tcga.exp[,used.ids]
	tcga.sp[[cl]] <- cor(t(tcga.exp.sub),method="spearman",use="everything")
}

tcga.exp <- tcga.m3[match(comm.gns,tcga.syms),] # 778 genes, 1093 samples
rownames(tcga.exp) <- comm.gns

##
setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/figs")
for(cl in smpl.class){
	cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
	png(paste0("tcga_",cl,"_signature.png"),width=700,height=700)
	par(cex.main=3)
	heatmap.2(tcga.sp[[cl]],#Colv=F,Rowv=F,dendrogram="none",
		RowSideColors=lcol,ColSideColors=lcol,hclustfun=hc,
		trace="none",col=cols,main=cl,key.title="")
	dev.off()
}

## METABRIC - subtypes
meta.exp <- meta.m3[match(comm.gns,meta.syms),] # 778 genes, 1440 samples
rownames(meta.exp) <- comm.gns

meta.sp <- list()
for(cl in smpl.class){
	used.ids <- meta.ids[[cl]]
	meta.exp.sub <- meta.exp[,used.ids]
	meta.sp[[cl]] <- cor(t(meta.exp.sub),method="spearman",use="everything")
}


setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/figs")
for(cl in smpl.class){
	cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
	png(paste0("meta_",cl,"_signature.png"),width=700,height=700)
	par(cex.main=3)
	heatmap.2(meta.sp[[cl]],#Colv=F,Rowv=F,dendrogram="none",
		RowSideColors=lcol,ColSideColors=lcol,hclustfun=hc,
		trace="none",col=cols,main=cl,key.title="")
	dev.off()
}


## TCGA
tcga.sp.panel <- list()
tcga.sp.panel$esr <- tcga.sp.panel$apm.tc <- tcga.sp.panel$between <- list()

for(cl in smpl.class[c(1,3,2,4)]){
	tmp <- tcga.sp[[cl]]
	g1 <- tmp[comm.esr,comm.esr] # 378
	g1 <- g1[lower.tri(g1)]
	g2 <- tmp[comm.apm.tc,comm.apm.tc] # 9045
	g2 <- g2[lower.tri(g2)]
	g3 <- as.vector(tmp[comm.esr,comm.apm.tc])# 3780

	tcga.sp.panel$esr[[cl]] <- g1
	tcga.sp.panel$apm.tc[[cl]] <- g2
	tcga.sp.panel$between[[cl]] <- g3
}

names(tcga.sp.panel$esr) <- names(tcga.sp.panel$apm.tc) <- names(tcga.sp.panel$between) <- c("HR","DP","HER2","TN")

tcga.esr <- data.frame(
	cor=unlist(tcga.sp.panel$esr),
	class=factor(rep(c("HR","DP","HER2","TN"),each=378),levels=c("HR","DP","HER2","TN")))

tcga.apm.tc <- data.frame(
	cor=unlist(tcga.sp.panel$apm.tc),
	class=factor(rep(c("HR","DP","HER2","TN"),each=9045),levels=c("HR","DP","HER2","TN")))

tcga.between <- data.frame(
	cor=unlist(tcga.sp.panel$between),
	class=factor(rep(c("HR","DP","HER2","TN"),each=3780),levels=c("HR","DP","HER2","TN")))

plot.cor <- function(x,title="TCGA",ylab){
	ggplot(x,aes(class,cor)) + 
	geom_violin(fill='lightblue') +
	# geom_jitter(height = 0, width = 0.2,size=.5) +
	labs(title=title,
        x ="", y = ylab) +
    theme_bw() +
	geom_hline(yintercept = 0) +
    theme(plot.title=element_text(hjust=0.5),text = element_text(size = 20))
}

setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/figs")
png("tcga_esr_cor.png",width=400,height=700)
plot.cor(tcga.esr,title="TCGA",ylab="Intra-cluster correlation coefficients (ESR)")
dev.off()
png("tcga_apm-tc_cor.png",width=400,height=700)
plot.cor(tcga.apm.tc,title="TCGA",ylab="Intra-cluster correlation coefficients (APM/TC)")
dev.off()
png("tcga_between_cor.png",width=400,height=700)
plot.cor(tcga.between,title="TCGA",ylab="Inter-cluster correlation coefficients (ESR vs APM/TC)")
dev.off()

p1 <- sapply(1:4,function(i)t.test(tcga.sp.panel$esr[[i]],alternative="greater")$p.value)
p2 <- sapply(1:4,function(i)t.test(tcga.sp.panel$apm.tc[[i]],alternative="greater")$p.value)
p3 <- sapply(1:4,function(i)t.test(tcga.sp.panel$between[[i]],alternative="less")$p.value)

## METABRIC
meta.sp.panel <- list()
meta.sp.panel$esr <- meta.sp.panel$apm.tc <- meta.sp.panel$between <- list()

for(cl in smpl.class[c(1,3,2,4)]){
	tmp <- meta.sp[[cl]]
	g1 <- tmp[comm.esr,comm.esr] # 378
	g1 <- g1[lower.tri(g1)]
	g2 <- tmp[comm.apm.tc,comm.apm.tc] # 9045
	g2 <- g2[lower.tri(g2)]
	g3 <- as.vector(tmp[comm.esr,comm.apm.tc])# 3780

	meta.sp.panel$esr[[cl]] <- g1
	meta.sp.panel$apm.tc[[cl]] <- g2
	meta.sp.panel$between[[cl]] <- g3
}

names(meta.sp.panel$esr) <- names(meta.sp.panel$apm.tc) <- names(meta.sp.panel$between) <- c("HR","DP","HER2","TN")

meta.esr <- data.frame(
	cor=unlist(meta.sp.panel$esr),
	class=factor(rep(c("HR","DP","HER2","TN"),each=378),levels=c("HR","DP","HER2","TN")))

meta.apm.tc <- data.frame(
	cor=unlist(meta.sp.panel$apm.tc),
	class=factor(rep(c("HR","DP","HER2","TN"),each=9045),levels=c("HR","DP","HER2","TN")))

meta.between <- data.frame(
	cor=unlist(meta.sp.panel$between),
	class=factor(rep(c("HR","DP","HER2","TN"),each=3780),levels=c("HR","DP","HER2","TN")))

setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/figs")
png("meta_esr_cor.png",width=400,height=700)
plot.cor(meta.esr,title="METABRIC",ylab="Intra-cluster correlation coefficients (ESR)")
dev.off()
png("meta_apm-tc_cor.png",width=400,height=700)
plot.cor(meta.apm.tc,title="METABRIC",ylab="Intra-cluster correlation coefficients (APM/TC)")
dev.off()
png("meta_between_cor.png",width=400,height=700)
plot.cor(meta.between,title="METABRIC",ylab="Inter-cluster correlation coefficients (ESR vs APM/TC)")
dev.off()

p1 <- sapply(1:4,function(i)t.test(meta.sp.panel$esr[[i]],alternative="greater")$p.value)
p2 <- sapply(1:4,function(i)t.test(meta.sp.panel$apm.tc[[i]],alternative="greater")$p.value)
p3 <- sapply(1:4,function(i)t.test(meta.sp.panel$between[[i]],alternative="less")$p.value)

## ESR1 expression 
## TCGA
esr.exp <- tcga.m3[match("ESR1",tcga.syms),unlist(tcga.ids[c(1,3,2,4)])]
cls <- factor(rep(c("HR","DP","HER2","TN"),sapply(tcga.ids[c(1,3,2,4)],length)),
	levels=c("HR","DP","HER2","TN"))

df <- data.frame(esr1=esr.exp,class=cls)

p <- ggplot(df,aes(class,esr1)) + 
	geom_violin(fill='grey80') +
	geom_jitter(height = 0, width = 0.2,size=.5) +
	labs(title="TCGA",
        x ="", y = "ESR1 expression (log, z-score)") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5),text = element_text(size = 20))
p	

## METABRIC
esr.exp <- meta.m3[match("ESR1",meta.syms),unlist(meta.ids[c(1,3,2,4)])]
cls <- factor(rep(c("HR","DP","HER2","TN"),sapply(meta.ids[c(1,3,2,4)],length)),
	levels=c("HR","DP","HER2","TN"))

df <- data.frame(esr1=esr.exp,class=cls)

p <- ggplot(df,aes(class,esr1)) + 
	geom_violin(fill='grey80') +
	geom_jitter(height = 0, width = 0.2,size=.5) +
	labs(title="METABRIC",
        x ="", y = "ESR1 expression (log, z-score)") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5),text = element_text(size = 20))
p	



