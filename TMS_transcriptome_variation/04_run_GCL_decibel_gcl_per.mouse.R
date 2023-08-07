library(zellkonverter)
library(SingleCellExperiment)
library(tidyverse)
library(ggExtra) 
library(gridExtra)
library(scRNAseq)
library(ggplot2);theme_set(theme_classic())
library(scran)
library(ggpubr)
library(Seurat)
library(grDevices);library(RColorBrewer)
library(SingleCellExperiment)
library(ggpmisc) #https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
#https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
scaleFUN <- function(x) sprintf("%.3f", x)
mycols=c("#999999", "#E69F00");
#mycols=c("#999999",  "#56B4E9","#E69F00");
plotcol <- brewer.pal(8,'Dark2')
library(viridis)

age.groups=c('3m','24m')
min.ncell_per.mouse=15;

cell.lifespan=data.table::fread('../src/Dataset_S1.txt')
cell.lifespan$Ron_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`

sce.minCellCounts=readRDS('sce_minCellCounts.rds')
tc.names=names(sce.minCellCounts) 
length(tc.names) #39

tc.names=intersect(tc.names,cell.lifespan$`cell type annotation in TMS`);
length(tc.names) #39
sce.minCellCounts=sce.minCellCounts[tc.names]

out_folder='sce.minCellCounts_gcl_per.mouse/'; dir.create(out_folder)
out_plot='sce.minCellCounts_gcl_per.mouse.pdf'
out_file='sce.minCellCounts_gcl_per.mouse.txt'

npoints<-lapply(tc.names,function(tc){
  sce.tc=sce.minCellCounts[[tc]]
  #ages=names(sce.tc)
  ages=unique(sce.tc$age)
  
  df=colData(sce.tc)
  df$mouse.id_age=paste(df$mouse.id,df$age)
  
  # sample equal number of cells
  mouse.id.names=names(which(table(df$mouse.id_age)>=min.ncell_per.mouse))
  table(df[df$mouse.id_age %in% mouse.id.names,]$mouse.id_age)
  
  res.out<-lapply(mouse.id.names,function(i){
    tmp=sce.tc[,df$mouse.id_age==i]
    if(T){ #use all expr genes
      csv_mat=assay(tmp,'counts')
      csv_mat=csv_mat[Matrix::rowSums(csv_mat>0)!=0,,drop=FALSE] #filter gene
      csv_mat=csv_mat[,Matrix::colSums(csv_mat>0)>=100,,drop=FALSE] #filter cell
      dim(csv_mat) 
    }
    if(ncol(csv_mat)<min.ncell_per.mouse){return(NULL)}
    #cat('tc',tc,'of mouse.id',i,'is done\n')
    return(i)
  })
  return(c(tc,length(res.out)))
})
as.data.frame(Reduce(`rbind`,npoints[order(tc.names)]))

########################################################################
## GCL code from decibel package: https://gitlab.com/olgaibanez/decibel
## https://gitlab.com/olgaibanez/decibel/-/blob/main/module/decibel.py
## def gcl(adata, num_divisions)
## save this function python code as 'decibel_module_gcl.py' 
if(TRUE){  #run once
    
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
  
  start=Sys.time();
  
  tc.names=names(sce.minCellCounts)
  num_divisions=as.integer(100) 
  
  for(tc in tc.names){
    
    sce.tc=sce.minCellCounts[[tc]]
    ages=unique(sce.tc$age)
    
    df=colData(sce.tc)
    df$mouse.id_age=paste(df$mouse.id,df$age)
    
    # sample equal number of cells
    mouse.id.names=names(which(table(df$mouse.id_age)>=min.ncell_per.mouse))
    table(df[df$mouse.id_age %in% mouse.id.names,]$mouse.id_age)
    ncell=min(table(df[df$mouse.id_age %in% mouse.id.names,]$mouse.id_age))
    
    (out.file=paste0(out_folder,'/',gsub(':','__',tc),'.txt'))
    if(file.exists(out.file)){cat('tc',tc,'out.file exists\n');next}
        
    res.out<-lapply(mouse.id.names,function(i){
      tmp=sce.tc[,df$mouse.id_age==i]
      if(T){ #use all expr genes
        csv_mat=assay(tmp,'counts')
        csv_mat=csv_mat[Matrix::rowSums(csv_mat>0)!=0,,drop=FALSE] #filter gene
        csv_mat=csv_mat[,Matrix::colSums(csv_mat>0)>=100,,drop=FALSE] #filter cell
        dim(csv_mat) 
      }
      if(ncol(csv_mat)<min.ncell_per.mouse){return(NULL)}
      
      csv_mat=as.matrix(csv_mat) 
      res.values=gcl_lib$gcl(csv_mat, num_divisions)
      
      res=data.frame(tc=tc,gcl=unlist(res.values),
                     age=rep(tmp$age[[1]],num_divisions),'mouse.id'=rep(tmp$mouse.id[[1]],num_divisions))
      cat('tc',tc,'of mouse.id',i,'is done\n')
      return(res)
    })
    
    res.out<-Filter(Negate(is.null), res.out)
    res.df.out=as.data.frame(Reduce(`rbind`,res.out))
    data.table::fwrite(res.df.out,out.file)
  }
  
  end=Sys.time();
  print(end-start) 

}

############################################################
(files=Sys.glob(paste0('sce.minCellCounts_gcl_per.mouse/*txt')))
x=lapply(files,function(file){x=data.table::fread(file)})
res.df=as.data.frame(Reduce(`rbind`,x))

res.df %>% group_by(tc,age) %>% summarise(nrep=n())

res.df$age=factor(res.df$age,levels=age.groups)

res.df=res.df[res.df$tc %in% cell.lifespan$`cell type annotation in TMS`,]
dim(table(res.df$tc,res.df$age)) #39

res.df$cell.type=res.df$tc
res.df$cell.dist=1-res.df$gcl
out_file="sce.minCellCounts_gcl_per.mouse.txt"
data.table::fwrite(res.df,out_file)

########################################################################
res.df=data.table::fread(out_file)
res.df$age=factor(res.df$age,levels=age.groups)

tmp0=res.df %>% group_by(cell.type,age,mouse.id) %>% summarise(age.average.per.mouse=mean(cell.dist))
tmp=tmp0 %>% group_by(cell.type,age) %>% summarise(age.average=mean(age.average.per.mouse))

tmp1 = tmp %>% spread(age,age.average)
sum(is.na(tmp1[[age.groups[[1]]]]))
sum(is.na(tmp1[[age.groups[[2]]]]))
tmp1$log2ratio=log2(tmp1[[age.groups[[2]]]]/tmp1[[age.groups[[1]]]])
sum(tmp1$log2ratio>0,na.rm=T) 

plots <- lapply(unique(tmp1$cell.type),function(i){
    df=res.df[res.df$cell.type==i,]
    df=df[order(df$age),]
    df$mouse.id=factor(df$mouse.id,unique(df$mouse.id))
    ratio.value=tmp1[tmp1$cell.type==i,]$log2ratio
    #ggplot(df,aes(x=age,y=cell.dist))+
    ggplot(df,aes(x=mouse.id,y=cell.dist,col=age,fill=age))+
      geom_violin()+#geom_jitter(size=0.1)+
      stat_summary(fun.y=median, geom="point", size=1, color="black")+
      theme_classic()+ylab('1-GCL' )+
      ggtitle(paste0(i,'\nlog2(O/Y)=',sprintf('%.3f',ratio.value)))+
      scale_color_manual(name='Month',values=mycols)+
      scale_fill_manual(name='Month',values=mycols)+
      theme(plot.title = element_text(size=9),legend.position = 'none',
            axis.text.x = element_text(size=6))+
      xlab('')+ylab('')
})
length(plots)
  
pdf(out_plot,useDingbats = T,width = 12,height = 17)
grid.arrange(grobs=plots,ncol=4)
dev.off()
