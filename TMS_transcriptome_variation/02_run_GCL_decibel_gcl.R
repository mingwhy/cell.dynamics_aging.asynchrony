
library(zellkonverter)
library(SingleCellExperiment)
library(tidyverse)
library(ggExtra) #
library(scRNAseq)
library(ggplot2);theme_set(theme_classic())
library(scran)
library(ggpubr)
library(Seurat)
library(grDevices);library(RColorBrewer)
library(SingleCellExperiment)
plotcol <- brewer.pal(8,'Dark2')
library(viridis)

mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.6)),
  colhead = list(fg_params=list(cex = 0.6)),
  rowhead = list(fg_params=list(cex = 0.6)))

source('src_cellular_lifespan.R')
sce.minCellCounts=readRDS('sce_minCellCounts.rds')
tc.names=names(sce.minCellCounts) 
tc.names=intersect(tc.names,cell.lifespan$`cell type annotation in TMS`);
sce.minCellCounts=sce.minCellCounts[tc.names]

out_folder='sce.minCellCounts_gcl'; dir.create(out_folder)
out_plot='sce.minCellCounts_gcl.pdf'
out_file='sce.minCellCounts_gcl.txt'

########################################################################
## GCL code from decibel package: https://gitlab.com/olgaibanez/decibel
## https://gitlab.com/olgaibanez/decibel/-/blob/main/module/decibel.py
## def gcl(adata, num_divisions)
## save this function python code as 'decibel_module_gcl.py' 
if(F){  #run once
  
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/Users/mingyang/anaconda3/envs/myenv/bin/python")
py_config()
sys=import('sys')
sys$path
# make sure you have 'decibel_module_gcl.py' in the working env
sys$path=c(sys$path,'./') #sys.path.append('./decibel/module/')
sys$path 
gcl_lib=import('decibel_module_gcl') 
gcl_lib$gcl

py_run_string('from pathlib import Path')
py_run_string('import sys')
py_run_string('from pathlib import Path')
py_run_string('from os import path')
# for threading:
py_run_string('import threading')
# math tools:
np=import('numpy')
py_run_string('from scipy.spatial.distance import cdist')
py_run_string('import math')
py_run_string('import random')
# gcl library import:
py_run_string('from matplotlib import pyplot as plt')

######################################################################
## run

# check #gene per tc
ngenes<-lapply(tc.names,function(tc){
  sce.tc=sce.minCellCounts[[tc]]
  ages=unique(sce.tc$age)
  
   tmp=sapply(ages,function(age){
      tmp=sce.tc[,sce.tc$age==age]
      csv_mat=assay(tmp,'counts')
      csv_mat=csv_mat[Matrix::rowSums(csv_mat>0)>=10,] #filter gene
      csv_mat=csv_mat[,Matrix::colSums(csv_mat>0)>=100] #filter cell
      nrow(csv_mat) #(3000, 202), 3000 gene by 202 cell
  })
  return(tmp)
})
df.ngenes=Reduce(`rbind`,ngenes)
rownames(df.ngenes)<-tc.names

##
start=proc.time();

tc.names=names(sce.minCellCounts)
num_divisions=as.integer(100) 

for(tc in tc.names){
  
  sce.tc=sce.minCellCounts[[tc]]
  ages=unique(sce.tc$age)
  
  (out.file=paste0(out_folder,'/',gsub(':','__',tc),'.txt'))
  if(file.exists(out.file)){cat('tc',tc,'out.file exists\n');next}
  
  res.df=data.frame();
  for(age in ages){
    tmp=sce.tc[,sce.tc$age==age]
    csv_mat=assay(tmp,'counts')
    csv_mat=csv_mat[Matrix::rowSums(csv_mat>0)>=10,] #filter gene
    csv_mat=csv_mat[,Matrix::colSums(csv_mat>0)>=100] #filter cell
    csv_mat=as.matrix(csv_mat) 
    res.values=gcl_lib$gcl(csv_mat, num_divisions)
    res=data.frame(tc=tc,gcl=unlist(res.values),age=rep(age,num_divisions))
    res.df=rbind(res.df,res)
  }
  data.table::fwrite(res.df,out.file)
}
print(proc.time()-start) 
}

########################################################################
(files=Sys.glob(paste0(out_folder,'/*txt')))
x=lapply(files,function(file){x=data.table::fread(file)})
res.df=as.data.frame(Reduce(`rbind`,x))

#res.df$age=factor(res.df$age,levels=c('3m','18m','24m'))
res.df$age=factor(res.df$age,levels=c('3m','24m'))

res.df=res.df[res.df$tc %in% cell.lifespan$`cell type annotation in TMS`,]
dim(table(res.df$tc,res.df$age)) #39

res.df$cell.type=res.df$tc
res.df$cell.dist=1-res.df$gcl

tmp=res.df %>% group_by(cell.type,age) %>% summarise(age.median=median(cell.dist))
tmp1 = tmp %>% spread(age,age.median)
sum(is.na(tmp1$`24m`))
tmp1$log2ratio=log2(tmp1$`24m`/tmp1$`3m`)


plots <- lapply(tc.names,function(i){
  df=res.df[res.df$cell.type==i,]
  ratio.value=tmp1[tmp1$cell.type==i,]$log2ratio
  ggplot(df,aes(x=age,y=cell.dist))+
    geom_violin()+geom_jitter(size=0.1)+
    stat_summary(fun.y=median, geom="point", size=2, color="red")+
    theme_classic()+ylab('1-GCL' )+
    ggtitle(paste0(i,'\nlog2(O/Y)=',round(ratio.value,5)))+
    theme(plot.title = element_text(size=11))
})
length(plots)


pdf(out_plot,useDingbats = T,width = 14,height = 18)
grid.arrange(grobs=plots,ncol=4)
dev.off()

data.table::fwrite(res.df,out_file)

