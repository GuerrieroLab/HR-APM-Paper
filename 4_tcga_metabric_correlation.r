library(RColorBrewer)
library(dendextend)
library(gplots)

##
setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/RData")
x <- load(file="gene_list.rda") # gene.sig, category, hallmarks
y <- load(file="tcga_RNA-seq.rda") # tcga.ids,tcga.clinical,tcga.syms,tcga.m1,tcga.m3
z <- load(file="metabric_subset.rda") # meta.ids,meta.clinical,meta.syms,meta.m1,meta.m3

sig.gns <- unique(unlist(gene.sig)) ## 988 genes in either of the three signatures
used.gns <- intersect(unlist(gene.sig),intersect(tcga.syms,meta.syms)) ## 988 -> 778
tested.gns  <- intersect(tcga.syms,meta.syms)

if(0){
	saveRDS(used.gns,file="common_sig_gns.rds")
	saveRDS(tested.gns,file="common_tested_gns.rds")
}

## overlap between genes
all.gns <- unique(c(tcga.syms,meta.syms,unlist(gene.sig)))
limma::vennDiagram(cbind(
	`TCGA`= all.gns %in% tcga.syms,
	`META`= all.gns %in% meta.syms,
	`GeneSig`= all.gns %in% unlist(gene.sig)),
main="Overlap of genes")

##
table(table(tcga.syms[tcga.syms %in% used.gns])) ## all unique
tcga.exp <- tcga.m3[match(used.gns,tcga.syms),] # 778 genes, 1093 samples
rownames(tcga.exp) <- used.gns

##
smpl.class <- names(tcga.ids)
tcga.sp <- list()
for(cl in smpl.class){
	used.ids <- tcga.ids[[cl]]
	tcga.exp.sub <- tcga.exp[,used.ids]
	tcga.sp[[cl]] <- cor(t(tcga.exp.sub),method="spearman",use="everything")
}

sapply(tcga.sp,function(x)sum(is.na(x))) ## no missing values

## 
table(table(meta.syms[meta.syms %in% used.gns])) ## 8 duplicated, leave it for now (use only the first one)

meta.exp <- meta.m3[match(used.gns,meta.syms),] # 778 genes, 1093 samples
rownames(meta.exp) <- used.gns

##
smpl.class <- names(meta.ids)
meta.sp <- list()
for(cl in smpl.class){
	used.ids <- meta.ids[[cl]]
	meta.exp.sub <- meta.exp[,used.ids]
	meta.sp[[cl]] <- cor(t(meta.exp.sub),method="spearman",use="everything")
}

sapply(meta.sp,function(x)sum(is.na(x))) ## no missing values
plot(t(meta.m1[match(c("ESR1","HLA-A"),meta.syms),meta.ids$hr]),pch=20) # no null values

##
hc <- function(x)hclust(x,method="complete") # switch from 'complete' to 'average', 'ward.D2', etc.

myheatmap <- function(x,k=10,side.cols=lcols){
	cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50) # heatmap

	Rowv  <- x %>% dist %>% hc %>% as.dendrogram %>%
		set("branches_k_color", k = k) %>% 
		set("branches_lwd", 1) %>% ladderize

	Colv  <- x %>% t %>% dist %>% hc %>% as.dendrogram %>%
		set("branches_k_color", k = k) %>%
		set("branches_lwd", 1) %>% ladderize

	st <- heatmap.2(x,
		trace="none",
		col=cols,
		Rowv=Rowv,
		Colv=Colv,
		RowSideColors=side.cols,
		ColSideColors=side.cols,
		hclustfun=hc,
		cexRow=1,
		cexCol=1
	)
	return(invisible(st))
}

##
category.tested <- category[used.gns]
lcols <- brewer.pal(6,"Set2")[as.numeric(category.tested)] ## side colors
tcga.sp1 <- tcga.sp$hr
txt <- rownames(tcga.sp1)
txt[!txt %in% c("ESR1","FOXA1","CD8A","CD8B","HLA-A","NLRC5")] <-NA
# txt.na <- txt
# txt.na[] <- NA

lcols2 <- c(NA,"red")[(!is.na(txt))+1]
colnames(tcga.sp1)  <- rownames(tcga.sp1)  <- txt

setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/figs")
png("tcga_er_v1.png",width=700,height=700)
myheatmap(tcga.sp1,k=9,side.cols=lcols)
dev.off()

png("tcga_er_v2.png",width=700,height=700)
tcga.hc <- myheatmap(tcga.sp1,k=9,side.cols=lcols2)
dev.off()

t1 <- txt[tcga.hc$rowInd]
t1[!is.na(t1)]

## cluster memberships
hc <- function(x)hclust(x,method="complete") 

hcm <- hc(dist(tcga.sp$hr))
dend <- as.dendrogram(hcm)

ks <-cutree(hcm,k=9)

new.order <- rev(unique(ks[tcga.hc$rowInd]))
ks <- factor(as.numeric(factor(ks,levels=new.order)))
names(ks) <- used.gns

tcga.clust.gns <- tapply(used.gns,ks,identity)
names(tcga.clust.gns) <- seq(tcga.clust.gns)

setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/RData")
saveRDS(tcga.clust.gns,"tcga.clust.gns.rds")

##
category.tested <- category[used.gns]
lcols <- brewer.pal(6,"Set2")[as.numeric(category.tested)] ## side colors
meta.sp1 <- meta.sp$hr
txt <- rownames(meta.sp1)
txt[!txt %in% c("ESR1","FOXA1","CD8A","CD8B","HLA-A","NLRC5")] <-NA
txt.na <- txt
txt.na[] <- NA

lcols2 <- c(NA,"red")[(!is.na(txt))+1]
colnames(meta.sp1)  <- rownames(meta.sp1)  <- txt

setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/figs")
png("meta_er_v1.png",width=700,height=700)
myheatmap(meta.sp1,k=9,side.cols=lcols)
dev.off()

png("meta_er_v2.png",width=700,height=700)
meta.hc <- myheatmap(meta.sp1,k=9,side.cols=lcols2)
dev.off()

t1 <- txt[meta.hc$rowInd]
t1[!is.na(t1)]

## cluster memberships
hc <- function(x)hclust(x,method="complete") 

hcm <- hc(dist(meta.sp$hr))
dend <- as.dendrogram(hcm)

ks <-cutree(hcm,k=9)

new.order <- rev(unique(ks[meta.hc$rowInd]))
ks <- factor(as.numeric(factor(ks,levels=new.order)))
names(ks) <- used.gns

meta.clust.gns <- tapply(used.gns,ks,identity)
names(meta.clust.gns) <- seq(meta.clust.gns)

setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/RData")
saveRDS(meta.clust.gns,"meta.clust.gns.rds")

##
plot.new()
legend("topright",c("ESR1","APM","TC","ESR1/TC","ESR1/APM","APM/TC"),
	fill=brewer.pal(6,"Set2"),ncol=2,border=NA,box.lwd=0)

