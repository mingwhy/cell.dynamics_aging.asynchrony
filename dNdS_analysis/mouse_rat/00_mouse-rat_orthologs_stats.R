
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
filter_mouse_rat=data.table::fread('../src/mouse_rat.dnds.txt')

head(filter_mouse_rat)
dim(filter_mouse_rat) #30033 genes
dim(id.mapping) #20449
sum(filter_mouse_rat$ensembl_gene_id %in% id.mapping$ensembl_gene_id) #20697

gene.meta=merge(filter_mouse_rat,id.mapping,by.x='ensembl_gene_id',by.y='ensembl_gene_id')
gene.meta$omega=gene.meta$rnorvegicus_homolog_dn/gene.meta$rnorvegicus_homolog_ds
gene.meta[is.infinite(gene.meta$omega),]$omega=NA
gene.meta[is.nan(gene.meta$omega),]$omega=NA
gene.meta=gene.meta[!is.na(gene.meta$omega),]
dim(gene.meta)
#[1] 19959    12

gene.meta.1to1=gene.meta[gene.meta$rnorvegicus_homolog_orthology_type=='ortholog_one2one',]
dim(gene.meta.1to1) #14185
sum(gene.meta.1to1$omega>=1) #76

table(gene.meta$rnorvegicus_homolog_orthology_type)
#ortholog_many2many  ortholog_one2many   ortholog_one2one 
#3383               2391              14185 
summary(gene.meta$omega)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.001136 0.074084 0.168511 0.256657 0.352365 3.568905

summary(gene.meta$omega) #non-zero
gene.meta2<-gene.meta %>% group_by(rnorvegicus_homolog_orthology_type) %>% 
  summarise(n = n(),median=median(omega))
max(gene.meta$omega)

pdf('mouse_rat_orthologs.pdf',useDingbats = T,height = 7)
ggplot(gene.meta,aes(x=`rnorvegicus_homolog_orthology_type`,y=omega))+geom_violin()+
  geom_jitter(size=0.02)+scale_y_log10(limits=c(0.001,4.3))+theme_classic(base_size = 12)+
  ylab('Nonsynonymous/synonymous substitution rate (dN/dS) ratio \nof mouse-rat orthologous genes')+
  xlab('Types of homolog orthologs')+
  stat_summary(fun=median,geom='point',col='red',size=2)+
  geom_text(data = gene.meta2, 
            aes(y = max(gene.meta$omega)*0.9, label = paste0(n, ' genes\n','median = ',round(median,3))),
            hjust = 0.5, vjust = 0, 
            position = position_nudge(y = 0.08), size = 4) 
dev.off()


