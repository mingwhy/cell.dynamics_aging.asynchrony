
library(gridExtra) 
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);library(grid)
library(tidyverse);
library(zellkonverter)
library(SummarizedExperiment)
library(ggpointdensity)
library(viridis);library(RColorBrewer);library(grDevices)
library(ggpubr)

###############################################################################################
## reading in gene id.mapping and dn.ds infomation
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt')  #readin_h5ad.R
filter_mouse_human=data.table::fread('../../src/mouse_human.dnds.txt')
head(filter_mouse_human)
dim(filter_mouse_human) #23089 genes
dim(id.mapping) #20449
sum(filter_mouse_human$ensembl_gene_id %in% id.mapping$ensembl_gene_id) #18311

gene.meta=merge(filter_mouse_human,id.mapping,by.x='ensembl_gene_id',by.y='ensembl_gene_id')
gene.meta$omega=gene.meta$hsapiens_homolog_dn/gene.meta$hsapiens_homolog_ds
gene.meta[is.infinite(gene.meta$omega),]$omega=NA
gene.meta[is.nan(gene.meta$omega),]$omega=NA
gene.meta=gene.meta[!is.na(gene.meta$omega),]
dim(gene.meta)
#[1] 18045    12

table(gene.meta$hsapiens_homolog_orthology_type)
#ortholog_many2many  ortholog_one2many   ortholog_one2one 
#2643               1870              13532 

gene.meta.1to1=gene.meta[gene.meta$hsapiens_homolog_orthology_type=='ortholog_one2one',]
dim(gene.meta.1to1) #13532
sum(gene.meta.1to1$omega>=1) #5

summary(gene.meta$omega) #non-zero
gene.meta2<-gene.meta %>% group_by(hsapiens_homolog_orthology_type) %>% 
  summarise(n = n(),median=median(omega))
max(gene.meta$omega)

pdf('mouse_human_orthologs.pdf',useDingbats = T,height = 7)
ggplot(gene.meta,aes(x=`hsapiens_homolog_orthology_type`,y=omega))+geom_violin()+
  geom_jitter(size=0.02)+scale_y_log10(limits=c(0.001,4))+theme_classic(base_size = 15)+
  ylab('Nonsynonymous/synonymous substitution rate (dN/dS) ratio \nof mouse-human orthologous genes')+
  xlab('Types of homolog orthologs')+
  stat_summary(fun=median,geom='point',col='red',size=2)+
  geom_text(data = gene.meta2, 
            aes(y = 2.5, label = paste0(n, ' genes\n','median = ',round(median,3))),
            hjust = 0.5, vjust = 0, 
            position = position_nudge(y = 0.08), size =5) 
dev.off()


