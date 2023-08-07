
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(zellkonverter)
library(SummarizedExperiment)
library(ggpointdensity)
library(viridis)
library(patchwork); #for plot_annotation

one.sex='male'
min.ncell=40
age.group=c(3,24)
#age.group=c(18,24)
age.group.names=paste0(age.group,'m')

######################################################################
if(!file.exists('select.tc.h5ad')){
  inp.sce<-readH5AD('~/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
  colData(inp.sce)$tissue_cell.type=paste(inp.sce$tissue,inp.sce$cell_ontology_class,sep=':')
  
  sub.sce=inp.sce[,inp.sce$sex==one.sex]
  x=table(sub.sce$tissue_cell.type,sub.sce$age) #3m 18m 21m 24m
  head(x)
  
  x=as.data.frame(x[,age.group.names])
  head(x)
  colnames(x)=c('tc','age','ncell')
  
  y=x %>% group_by(tc,age) %>% summarise(enough=ncell>=min.ncell) %>% 
    summarise(n.age.group=sum(enough))
  y=y[y$n.age.group==length(age.group.names),]
  dim(y) #58
  
  # only keep tc enough cells
  pick.cell.types=y$tc
  sce=sub.sce[,sub.sce$tissue_cell.type %in% pick.cell.types & sub.sce$age %in% age.group.names]
  
  unique(sce$age)
  sce$age=droplevels(sce$age)  
  unique(sce$age) 
  sce$age=factor(sce$age,levels=age.group.names)
    
  # integrate with cell lifespan
  cell.lifespan=data.table::fread('../src/Dataset_S1.txt')
  tc.names=unique(sce$tissue_cell.type)
  tc.names=intersect(tc.names,cell.lifespan$`cell type annotation in TMS`); 
  length(tc.names) #29tc
  
  sce=sce[,sce$tissue_cell.type %in% tc.names]
  
  unique(sce$tissue_cell.type) #29tc
  unique(sce$mouse.id)
  unique(sce$age)
  
  writeH5AD(sce, 'select.tc.h5ad')
}


