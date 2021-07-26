library(openxlsx)

## gene list
setwd("~/Dropbox (HMS)/_BTIL_/progress/hr_apm_paper/data/gene_list/Big Gene List/")

## apm
apm1 <- read.xlsx("APM gene list 2021.06.16.xlsx",1,colNames=F)
title.id <- which(is.na(apm1[[2]]))
titles <- apm1[[1]][title.id]
apm.list <- list()
for(i in seq(title.id)){
	start <- title.id[i]+1	
	end <- c(title.id,nrow(apm1)+1)[i+1]-1
	apm.list[[i]] <- apm1[start:end,]
}

apm.gns <- sapply(apm.list,function(x)x[[1]])
names(apm.gns) <- titles

if(0){
	apm2 <- read.xlsx("APM gene list 2021.06.16.xlsx",2)
	uniq.apm.gns <- unique(unlist(apm.gns))
	uniq.apm.gns[which(!uniq.apm.gns %in% apm2$Gene)] # why they're removed?, H2AZ1, UBE2D2, UBE2D3
}

##  esr
esr1 <- read.xlsx("ESR1 Gene lists 2021.06.16.xlsx",1,colNames=F)
title.id <- which(is.na(esr1[[2]]))
titles <- esr1[[1]][title.id]
esr.list <- list()
for(i in seq(title.id)){
	start <- title.id[i]+1	
	end <- c(title.id,nrow(esr1)+1)[i+1]-1
	esr.list[[i]] <- esr1[start:end,]
}

esr.gns <- sapply(esr.list,function(x)x[[1]][-1])
names(esr.gns) <- titles

if(0){
	esr2 <- read.xlsx("APM gene list 2021.06.16.xlsx",2)
	uniq.esr.gns <- unique(unlist(esr.gns))
	uniq.esr.gns[which(!uniq.esr.gns %in% esr2$Gene)] # why they're removed?, H2AZ1, UBE2D2, UBE2D3
}

##
sapply(esr.gns,function(x)table(table(x)))
table(table(unlist(esr.gns)))

## tc
tc <- read.xlsx("T cell gene list 2021.06.16.xlsx",1,colNames=T)
tc[[3]] <- sub("azizi et al 2018","azizi et al. 2018",tc[[3]])
source <- strsplit(tc$Literature,", ?")
src.names <- names(table(unlist(source)))
tc[[4]][grep("BIOCARTA",tc[[3]])]

tc.gns <- lapply(src.names,function(x)tc$Genes[grepl(x,tc[[3]])])
names(tc.gns) <- src.names

## looking into overlaps between genes
apm <- unique(unlist(apm.gns))
tc <- unique(unlist(tc.gns))
esr <- unique(unlist(esr.gns))

gene.sig <- list(esr=esr,apm=apm,tc=tc)
all.gns <- unique(unlist(gene.sig)) # 988 genes

ol <- sapply(gene.sig,function(x){
	all.gns %in% x
})
colnames(ol) <- c("ESR1","APM","TC")

limma::vennDiagram(ol,main="overlap")

##
idx <- ol %*% c(1,2,4)
uniq.category <- c("ESR1","APM","ESR1/APM","TC","ESR1/TC","APM/TC")
category <- factor(uniq.category[idx],levels=uniq.category[c(1,2,4,3,5,6)])
names(category) <- all.gns

hallmarks <- c("ESR1","FOXA1","CD8A","CD8B","HLA-A","NLRC5")
category[hallmarks]

##
setwd("~/Dropbox (HMS)/_BTIL_/progress/hr_apm_paper/data/RData")
save(gene.sig,category,hallmarks,file="gene_list.rda")

##
