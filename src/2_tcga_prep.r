library(dplyr)

## expression matrix
ma.fns <- c("data_RNA_Seq_v2_expression_median.txt",
	"data_RNA_Seq_v2_mRNA_median_Zscores.txt",
	"data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt")

setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/brca_tcga")
tcga.mat1 <- read.table(ma.fns[1],sep="\t",header=T)
tcga.mat2 <- read.table(ma.fns[2],sep="\t",header=T)
tcga.mat3 <- read.table(ma.fns[3],sep="\t",header=T)

## colnames - identical
tcga.m1 <- as.matrix(tcga.mat1[,-(1:2)])
tcga.m2 <- as.matrix(tcga.mat2[,-(1:2)])
tcga.m3 <- as.matrix(tcga.mat3[,-(1:2)])

dim(tcga.m1) # 20531 genes, 1100 samples
dim(tcga.m2) # 20440 genes, 1100 samples
dim(tcga.m3) # 20531 genes, 1100 samples

table(sub(".+\\.","",colnames(tcga.m3))) # 01 - primary, 06 - duplicated (metastasis?)

identical(tcga.mat1$Hugo_Symbol,tcga.mat3$Hugo_Symbol) # true

## convert column names
identical(colnames(tcga.mat1),colnames(tcga.mat2)) # true
identical(colnames(tcga.mat1),colnames(tcga.mat3)) # true

is.dup <- grepl("06$",colnames(tcga.m3))

tcga.m1 <- tcga.m1[,!is.dup]
tcga.m2 <- tcga.m2[,!is.dup]
tcga.m3 <- tcga.m3[,!is.dup]

colnames(tcga.m1) <- colnames(tcga.m2) <- colnames(tcga.m3) <- sub("\\.01$","",colnames(tcga.m3))

if(0){
	## they are monotonously increasing relationships, so doesn't matter which samples to use for computing spearman correlation
	i=4
	plot(tcga.m1[match(syms12,tcga.mat1$Hugo_Symbol)[i],],
		tcga.m2[match(syms12,tcga.mat2$Hugo_Symbol)[i],],
		pch=".")
	plot(tcga.m1[match(syms12,tcga.mat1$Hugo_Symbol)[i],],
		tcga.m3[match(syms12,tcga.mat2$Hugo_Symbol)[i],],
		pch=".")
	plot(tcga.m3[match(syms12,tcga.mat1$Hugo_Symbol)[i],],
		tcga.m2[match(syms12,tcga.mat2$Hugo_Symbol)[i],],
		pch=".")
}

## load clinical data
setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/brca_tcga")
cl <- readLines("data_bcr_clinical_data_patient.txt")[-(1:4)]
clt <- do.call(rbind,strsplit(cl,split="\t"))
cn <- clt[1,]
tcga.clinical <- data.frame(clt[-1,],stringsAsFactors=F)
names(tcga.clinical) <- cn
tcga.clinical$PATIENT_ID <- gsub("-",".",tcga.clinical$PATIENT_ID)

used.ids <- colnames(tcga.m3)
table(used.ids %in% tcga.clinical$PATIENT_ID) ## all samples' records exist in clinical metadata

tcga.clinical <- tcga.clinical %>% filter(PATIENT_ID %in% used.ids)

## # patients per disease
tcga.clin.1 <- tcga.clinical %>%
	filter(ER_STATUS_BY_IHC %in% c("Positive","Negative")) %>%
	filter(PR_STATUS_BY_IHC %in% c("Positive","Negative")) %>%
	filter(IHC_HER2 %in% c("Positive","Negative"))

mm <- cbind(
	`ER+`=(tcga.clin.1$ER_STATUS_BY_IHC=="Positive"),
	`PR+`=(tcga.clin.1$PR_STATUS_BY_IHC=="Positive"),
	`HER2+`=(tcga.clin.1$IHC_HER2=="Positive"))
rownames(mm) <- tcga.clin.1$PATIENT_ID

limma::vennDiagram(mm,main="TCGA samples")

tcga.er.pos.ids <- (tcga.clinical %>% 
	filter((ER_STATUS_BY_IHC == "Positive" | PR_STATUS_BY_IHC == "Positive")
			& IHC_HER2 == "Negative"))$PATIENT_ID # 445

tcga.her2.pos.ids <- (tcga.clinical %>% 
	filter(ER_STATUS_BY_IHC == "Negative" & 
		PR_STATUS_BY_IHC == "Negative" &
		IHC_HER2 == "Positive"))$PATIENT_ID # 37

tcga.dp.ids <- (tcga.clinical %>% 
	filter((ER_STATUS_BY_IHC == "Positive" | PR_STATUS_BY_IHC == "Positive")& 
		IHC_HER2 == "Positive"))$PATIENT_ID # 126

tcga.tnbc.ids <- (tcga.clinical %>% 
	filter(ER_STATUS_BY_IHC == "Negative" & 
		PR_STATUS_BY_IHC == "Negative" &
		IHC_HER2 == "Negative"))$PATIENT_ID # 37

tcga.all.ids <- tcga.clinical$PATIENT_ID

tcga.ids <- list(
	hr=tcga.er.pos.ids,
	her2=tcga.her2.pos.ids,
	dp=tcga.dp.ids,
	tnbc=tcga.tnbc.ids,
	all = tcga.all.ids
)
sapply(tcga.ids,length)

## Be careful, Hugo symbols (or Entrez IDs) are not unique.
table(table(tcga.mat1$Hugo)) # not unique
table(table(tcga.mat2$Hugo)) # not unique

tcga.syms <- tcga.mat3$Hugo_Symbol ## genes whose expression is measured

setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/RData")
save(tcga.ids,tcga.clinical,tcga.syms,tcga.m1,tcga.m3,file="tcga_RNA-seq.rda")
