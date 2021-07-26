library(dplyr)

##
setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/RData")
load(file="gene_list.rda") # gene.sig,category,hallmarks

## expression matrix
setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/metabric")

fns <-  c("data_expression_median.txt",
	"data_mRNA_median_Zscores.txt",
	"data_mRNA_median_all_sample_Zscores.txt")

meta.mat1 <- read.table(fns[1],sep="\t",header=T,stringsAsFactors=F)
meta.mat2 <- read.table(fns[2],sep="\t",header=T,stringsAsFactors=F)
meta.mat3 <- read.table(fns[3],sep="\t",header=T,stringsAsFactors=F)

## colnames - identical
meta.m1 <- as.matrix(meta.mat1[,-(1:2)])
meta.m2 <- as.matrix(meta.mat2[,-(1:2)])
meta.m3 <- as.matrix(meta.mat3[,-(1:2)])

dim(meta.m1) # 24368 genes, 1904 samples
dim(meta.m2) # 18543 genes, 1904 samples
dim(meta.m3) # 24368 genes, 1904 samples

identical(meta.mat1$Hugo_Symbol,meta.mat3$Hugo_Symbol) # true
syms12 <- intersect(meta.mat1$Hugo_Symbol,meta.mat2$Hugo_Symbol)

i=1
plot(meta.m1[match(syms12,meta.mat1$Hugo_Symbol)[i],],meta.m2[match(syms12,meta.mat2$Hugo_Symbol)[i],],pch=".")
plot(meta.m3[match(syms12,meta.mat1$Hugo_Symbol)[i],],meta.m2[match(syms12,meta.mat2$Hugo_Symbol)[i],],pch=".")

## they are monotonously increasing relationships, so doesn't matter which samples to use.


## load clinical data
setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/metabric")
meta.clinical <- read.table("data_clinical_patient.txt",sep="\t",header=T)
table(meta.clinical$ER_IHC) ## 1817

meta.clinical2 <- read.table("data_clinical_sample.txt",sep="\t",header=T)
table(meta.clinical2$ER_STATUS) ## 1825 - not used this one

meta.clinical <- meta.clinical %>% left_join(meta.clinical2,by="PATIENT_ID")

meta.clinical$PATIENT_ID <- gsub("-",".",meta.clinical$PATIENT_ID)

##
used.ids <- colnames(meta.m3)
table(used.ids %in% meta.clinical$PATIENT_ID) ## all samples' records exist in clinical metadata

meta.clinical <- meta.clinical %>% filter(PATIENT_ID %in% used.ids)

##
limma::vennDiagram(cbind(
	`ER+`=(meta.clinical$ER_STATUS=="Positive"),
	`PR+`=(meta.clinical$PR_STATUS=="Positive"),
	`HER2+`=(meta.clinical$HER2_STATUS=="Positive")),
main="METABRIC samples")

##
meta.clin.1 <- meta.clinical %>%
	filter(ER_STATUS %in% c("Positive","Negative")) %>%
	filter(PR_STATUS %in% c("Positive","Negative")) %>%
	filter(HER2_STATUS %in% c("Positive","Negative"))

meta.er.pos.ids <- (meta.clinical %>% 
	filter(
		(ER_STATUS == "Positive" | 
		PR_STATUS == "Positive")
			& HER2_STATUS == "Negative"))$PATIENT_ID # 445

meta.her2.pos.ids <- (meta.clinical %>% 
	filter(ER_STATUS == "Negative" & 
		PR_STATUS == "Negative" &
		HER2_STATUS == "Positive"))$PATIENT_ID # 37

meta.dp.ids <- (meta.clinical %>% 
	filter((ER_STATUS == "Positive" | 
		PR_STATUS == "Positive")
			& HER2_STATUS == "Positive"))$PATIENT_ID # 126

meta.tnbc.ids <- (meta.clinical %>% 
	filter(ER_STATUS == "Negative" & 
		PR_STATUS == "Negative" &
		HER2_STATUS == "Negative"))$PATIENT_ID # 37

meta.all.ids <- meta.clinical$PATIENT_ID

meta.ids <- list(
	hr=meta.er.pos.ids,
	her2=meta.her2.pos.ids,
	dp=meta.dp.ids,
	tnbc=meta.tnbc.ids,
	all = meta.all.ids
)

sapply(meta.ids,length)

##
meta.syms <- meta.mat3$Hugo_Symbol ## genes whose expression is measured

## HLA-bimodality
library(mclust)
set.seed(12345)
hlaa <- meta.m1[match("HLA-A",meta.syms),meta.ids$hr]
hlaa.all <- meta.m1[match("HLA-A",meta.syms),]

esr1 <- meta.m1[match("ESR1",meta.syms),]
er.status <- meta.clinical$ER_STATUS[match(names(esr1),meta.clinical$PATIENT_ID)]
gnn <- Mclust(hlaa,G=2)
hlaa.th <- mean(c(max(hlaa[gnn$cl==1]),min(hlaa[gnn$cl==2])))
if(0){
	boxplot(hlaa ~ gnn$classification)
	abline(h=hlaa.th,col=2)
}

if(0){
	par(mfrow=c(2,1))
	## 1
	x <- hist(meta.m1, breaks=1000,main="all genes",xlab="log2 expression",freq=F,border=NA,
		col="lightblue")
	mod <- x$mids[which.max(x$counts)]
	abline(v=mod,col=2)

	## 2
	rng <- range(meta.m1,na.rm=T)
	x <- hist(hlaa,breaks=seq(rng[1],rng[2],length=200),plot=F)
	plot(x,col="lightblue",freq=F,border=NA,main="HLA-A",xlim=rng,xlab="log2 expression")
	abline(v=mod,col=2)
	abline(v=hlaa.th,lty=2)

	## 3
	x1 <- hist(esr1[er.status=="Positive"],breaks=seq(rng[1],rng[2],length=200),plot=F)
	x2 <- hist(esr1[er.status=="Negative"],breaks=seq(rng[1],rng[2],length=200),plot=F)

	c1 <- rgb(220,110,110,max = 255, alpha = 80, names = "lt.blue")
	c2 <- rgb(80,160,225, max = 255, alpha = 80, names = "lt.pink")

	plot(x2,freq=F,border=NA,main="ESR1",xlim=rng,col=c2)
	plot(x1,freq=F,border=NA,main="ESR1",xlim=rng,col=c1,add=T)
	legend(12,1.2,c("ER positive","ER negative"),fill=c(c1,c2))

	##
	meta.hm.exp <- data.frame(t(meta.m1[match(hallmarks,meta.syms),meta.ids$hr]))
	names(meta.hm.exp) <- hallmarks

	hi.hlaa <- hlaa > hlaa.th

	setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/figs")
	png("metabric_hallmark-scatter.png",width=700,height=700)
	plot(meta.hm.exp,pch='.',col=hi.hlaa+1,main="METABRIC")
	dev.off()

	##
	setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/RData")
	load(file="tcga_RNA-seq.rda") # tcga.ids,tcga.clinical,tcga.syms,tcga.m1,tcga.m3

	tcga.hm.exp <- data.frame(t(tcga.m3[match(hallmarks,tcga.syms),tcga.ids$hr]))
	names(tcga.hm.exp) <- hallmarks

	png("tcga_hallmark-scatter.png",width=700,height=700)
	plot(tcga.hm.exp,pch='.',main="TCGA")
	dev.off()
}

##
hi.hlaa.all <- hlaa.all > hlaa.th
hi.hlaa.ids <- names(which(hi.hlaa.all))

meta.m1 <- meta.m1[,hi.hlaa.all]
meta.m3 <- meta.m3[,hi.hlaa.all]

meta.clinical <- meta.clinical[match(hi.hlaa.ids,meta.clinical$PATIENT_ID),]

meta.ids <- lapply(meta.ids,function(x){
	x[x %in% hi.hlaa.ids]
})

sapply(meta.ids,length)

##
hlaa <- meta.m1[match("HLA-A",meta.syms),meta.ids$hr]
hist(hlaa)

##
limma::vennDiagram(cbind(
	`ER+`=(meta.clinical$ER_STATUS=="Positive"),
	`PR+`=(meta.clinical$PR_STATUS=="Positive"),
	`HER2+`=(meta.clinical$HER2_STATUS=="Positive")),
main="METABRIC samples")

##
setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/RData")
save(meta.ids,meta.clinical,meta.syms,meta.m1,meta.m3,file="metabric_subset.rda")

##
library(Rtsne)
rt <- Rtsne(t(meta.m1),is_distance=FALSE,perplexity=10)
plot(rt$Y,col=idx)
plot(m1.all.er["HLA-A",])


