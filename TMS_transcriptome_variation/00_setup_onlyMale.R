
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(zellkonverter)
library(SummarizedExperiment)
library(ggpointdensity)
library(viridis)
library(patchwork); 

######################################################################
if(!file.exists('select.tc.h5ad')){
  #https://figshare.com/articles/dataset/tms_gene_data_rv1/12827615
  inp.sce<-readH5AD('TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
  colData(inp.sce)$tissue_cell.type=paste(inp.sce$tissue,inp.sce$cell_ontology_class,sep=':')
  
  sub.sce=inp.sce[,inp.sce$sex=='male']
  x=table(sub.sce$tissue_cell.type,sub.sce$age)
  x=as.data.frame(x[,c(1,4)]) #3m and 24m
  colnames(x)=c('tc','age','ncell')
  x=x %>% spread(age,ncell)
  dim(x) #202
  data.table::fwrite(x,'ncell_per.age_per.tc_raw.txt',sep='\t',quote=F)
  
  y=x[x[,2]>=50 & x[,3]>=50, ] #both age groups contain >=50 cells
  dim(y) #72 cell states
  sum(c(y$`3m`,y$`24m`)) #47898 cells in total
  
  y[grep('Brain',y$tc),]
  y$tc=as.character(y$tc)
  unique(unlist(lapply(strsplit(y$tc,'\\:'),'[',1))) #22 unique tissues
  unique(unlist(lapply(strsplit(y$tc,'\\:'),'[',2))) #52 unique cell types
  data.table::fwrite(y,'ncell_per.age_per.tc.txt',sep='\t',quote=F)
  
  
  # #only keep tc enough cells
  #sum(y$tc %in% tc.orders) 
  #pick.cell.types=y[y$tc %in% tc.orders,]$tc 
  pick.cell.types=y$tc
  sce=sub.sce[,sub.sce$tissue_cell.type %in% pick.cell.types]
  
  unique(sce$age) #'3m','18m','21m','24m'
  sce$age=droplevels(sce$age)  
  unique(sce$age) 
  sce$age=factor(sce$age,levels=c('3m','18m','24m'))
  cell.meta=colData(sce)
  #rowData(sce)
  length(table(cell.meta$tissue_cell.type)) 
  writeH5AD(sce, 'select.tc.h5ad')
}

y=data.table::fread('ncell_per.age_per.tc.txt')
y$tissue=sapply(strsplit(y$tc,':'),'[',1)
y$cell.type=sapply(strsplit(y$tc,':'),'[',2)
dim(y) #72
table(y$tissue) #22 tissues
length(table(y$cell.type)) #52 cell.type

