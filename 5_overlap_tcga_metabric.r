setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/RData")
load(file="gene_list.rda") # gene.sig,category,hallmarks
tcga.clust.gns <- readRDS("tcga.clust.gns.rds")
meta.clust.gns <- readRDS("meta.clust.gns.rds")
used.gns <- readRDS(file="common_sig_gns.rds")

tcga.esr <- tcga.clust.gns[[7]]
tcga.esr2 <- tcga.clust.gns[[4]]
tcga.apm.tc <- tcga.clust.gns[[3]]
meta.esr <- meta.clust.gns[[7]]
meta.esr2 <- meta.clust.gns[[5]]
meta.apm.tc <- unlist(meta.clust.gns[c(2,4)])

names(tcga.clust.gns) <- c(1:2,"3(apm.tc)",4:6,"7(esr)",8:9)
names(meta.clust.gns) <- c(1,"2(hlaa)",3,"4(apm.tc)",5:6,"7(esr)",8:9)

##
sapply(hallmarks,function(x)grep(x,tcga.clust.gns))
sapply(hallmarks,function(x)grep(x,meta.clust.gns))

if(0){
	nlps <- sapply(tcga.clust.gns,function(x){
		sapply(meta.clust.gns,function(y){
			tmp <- table(used.gns %in% x,used.gns %in% y)[c("TRUE","FALSE"),c("TRUE","FALSE")]
			nlp <- -log10(fisher.test(tmp,alternative="greater")$p.value)
		})
	})

	round(nlps)
}

##
esr.gns <- unique(c(tcga.esr,meta.esr))
esr2.gns <- unique(c(tcga.esr2,meta.esr2))
apm.tc.gns <- unique(c(tcga.apm.tc,meta.apm.tc))

mm1 <- sapply(list(tcga.apm.tc,meta.apm.tc),function(x)apm.tc.gns %in% x)
colnames(mm1) <- c("TCGA","METABRIC")
par(mar=c(0,0,0,0))
limma::vennDiagram(mm1,main="")
box()
mtext("APM/TC panel",side=3,line=1)
comm.apm.tc <- apm.tc.gns[which(rowSums(mm1)==2)]

mm2 <- sapply(list(tcga.esr,meta.esr),function(x)esr.gns %in% x)
colnames(mm2) <- c("TCGA","METABRIC")
limma::vennDiagram(mm2)
comm.esr <- esr.gns[which(rowSums(mm2)==2)]

mm3 <- sapply(list(tcga.esr2,meta.esr2),function(x)esr2.gns %in% x)
colnames(mm3) <- c("TCGA","METABRIC")
limma::vennDiagram(mm3)
comm.esr2 <- esr2.gns[which(rowSums(mm3)==2)]

fisher.test(matrix(c(28,41,27,778-(28+41+27)),ncol=2),alternative="greater")$p.value
fisher.test(matrix(c(112,29,39,778-(112+29+39)),ncol=2),alternative="greater")$p.value
fisher.test(matrix(c(135,19,66,778-(135+19+66)),ncol=2),alternative="greater")$p.value



if(0){
	setwd("/Users/kenichishimada/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/RData")
	save(comm.apm.tc,comm.esr,comm.esr2,file="common_signatures.rda")
}

##
if(0){
	save(tcga.esr,
		tcga.apm.tc,
		meta.esr,
		meta.apm,
		meta.tc,
		comm.esr1,
		comm.apm,
		comm.tc,
		file="tcga_metabric_signature.rda")
}else{
	setwd("/Users/kenichishimada/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/RData")
	load("tcga_metabric_signature.rda")
}

setwd("/Users/kenichishimada/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/RData")
load(file="gene_list.rda") # apm.gns,tc.gns,esr.gns

barplot(sapply(tcga.clust.gns,function(x){
	sum(x %in% unlist(esr.gns))/length(x)
}))  # 2

barplot(sapply(meta.clust.gns,function(x){
	sum(x %in% unlist(esr.gns))/length(x)
}))  # 9

esr2 <- intersect(tcga.clust.gns[[4]],meta.clust.gns[[5]])

##
load("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/RData/msigdb_v6.1.rda")
keggs <- msigdb[grep("^KEGG_",names(msigdb))]
gos <- msigdb[grep("^GO_",names(msigdb))]

setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/RData")
tested.gns <- readRDS(file="common_tested_gns.rds")

library(org.Hs.eg.db)
kg <- msigdb[sub("_.+","",names(msigdb)) %in% c("HALLMARK","GO","KEGG","REACTOME","BIOCARTA")]
barplot(sort(table(sub("_.+","",names(kg)))),las=2)

# kg <- c(keggs,gos)
kg.syms <- unlist(mget(unique(unlist(kg)),org.Hs.egSYMBOL,ifnotfound=NA))
all.genes <- intersect(kg.syms,tested.gns) # 13321
tab <- as.matrix(sort(table(sub("_.+","",names(kg))),decreasing=T))
colnames(tab)[1] <- "# gene sets"

kg.syms <- lapply(kg,function(x){
	s <- unlist(mget(x,org.Hs.egSYMBOL,ifnotfound=NA))
	s <- s[!is.na(s)]
	s[s %in% all.genes]
})

ge10 <- sapply(kg.syms,length)>=20 & sapply(kg.syms,length) <= 500
kg.syms <- kg.syms[ge10]

all.genes <- intersect(unlist(kg.syms),tested.gns) # 14787 (all pathways)

##
pvals.esr <- parallel::mclapply(kg.syms,function(x){
	tab <- table(all.genes %in% comm.esr, all.genes %in% x)[c("TRUE","FALSE"),c("TRUE","FALSE")]
	pval <- (fisher.test(tab,alternative="greater")$p.value)
},mc.cores=6)

pvals.esr2 <- parallel::mclapply(kg.syms,function(x){
	tab <- table(all.genes %in% esr2, all.genes %in% x)[c("TRUE","FALSE"),c("TRUE","FALSE")]
	pval <- (fisher.test(tab,alternative="greater")$p.value)
},mc.cores=6)

pvals.apm.tc <- parallel::mclapply(kg.syms,function(x){
	tab <- table(all.genes %in% comm.apm.tc, all.genes %in% x)[c("TRUE","FALSE"),c("TRUE","FALSE")]
	pval <- (fisher.test(tab,alternative="greater")$p.value)
},mc.cores=6)

##
library(qvalue)
nlqs.esr <- -log10(qvalue(unlist(pvals.esr))$qvalue)
nlqs.esr <- nlqs.esr[nlqs.esr > -log10(0.05)]
par(mar=c(5,30,5,1))
barplot(sort(nlqs.esr),horiz=T,col="steelblue",las=1,xlab="-log10(P-value, adusted)",
	main="Pathways enriched in ESR gene panel")
head(sort(nlqs.esr,decreasing=T),20)

nlqs.esr2 <- -log10(qvalue(unlist(pvals.esr2))$qvalue)
nlqs.esr2 <- nlqs.esr2[nlqs.esr2 > -log10(0.05)]
par(mar=c(5,30,5,1))
barplot(rev(head(sort(nlqs.esr2,decreasing=T),20)),
	horiz=T,col="steelblue",las=1,xlab="-log10(P-value, adusted)",
	main="Pathways enriched in ESR2 gene panel")
head(sort(nlqs.esr2,decreasing=T),20)

nlqs.apm.tc <- -log10(qvalue(unlist(pvals.apm.tc))$qvalue)
nlqs.apm.tc <- -log10(qvalue(unlist(pvals.apm.tc))$qvalue)
nlqs.apm.tc <- nlqs.apm.tc[nlqs.apm.tc > -log10(0.05)]
par(mar=c(5,30,5,1))
barplot(rev(head(sort(nlqs.apm.tc,decreasing=T),20)),
	horiz=T,col="steelblue",las=1,xlab="-log10(P-value, adusted)",
	main="Pathways enriched in APM/TC gene panel")
head(sort(nlqs.apm.tc,decreasing=T),20)

##
setwd("~/Dropbox (HMS)/_BTIL_/projects/hr_apm_paper/data/RData")
save.image("tmp_062621.rda")

