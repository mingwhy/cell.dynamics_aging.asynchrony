
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)
library(ggpointdensity)
library(viridis)
library(patchwork); #for plot_annotation
library(scMerge)
library(SeuratWrappers)


## subsample cell read count
library(countland) #https://github.com/shchurch/countland/blob/main/tutorials_and_vignettes/R_tutorials_and_vignettes/vignette-tutorial.Rmd

if(!file.exists('sce_minCellCounts.rds')){
  sce=readH5AD('select.tc.h5ad') 
  unique(sce$age)
  assayNames(sce)<-'counts'
  
  (tcs=unique(sce$tissue_cell.type))
  
  set.seed(84095) # set random seed for reproducibility
  minCellCounts<-lapply(tcs,function(tc){
    tmp=sce[,sce$tissue_cell.type==tc]
    #summary(tmp$n_counts)
    m=assay(tmp,'counts')
    
    C <- countland(m)
    C <- Subsample(C,cell_counts='min')
    return(C)
  })
  
  Matrix::colSums(minCellCounts[[1]]@subsample)
  minCellCounts[[1]]@subsample[1:3,1:3]
  names(minCellCounts)<-tcs
  
  # create single cell object
  tc.names=names(minCellCounts)
  
  ## filter gene (expr in at least 10% cells per cell type) and filter cell (expr >=100 genes)
  sce.minCellCounts<-lapply(tc.names,function(x){
    tmp=minCellCounts[[x]]@subsample
    
    # filter gene
    n.expr.cell = Matrix::rowSums(tmp>0)
    tmp=tmp[n.expr.cell>= 0.1*ncol(tmp),]
    
    # filter cell
    n.expr.gene=Matrix::colSums(tmp>0)
    mat=tmp[,n.expr.gene>=100] #expr at least 100 genes
    
    cell.info=colData(sce[,sce$tissue_cell.type==x])
    cell.info=cell.info[n.expr.gene>=100,] #discard cells
    
    cat(x);print(table(cell.info$age)) 
    
    sce.filtered <- SingleCellExperiment(list(counts=mat),colData=cell.info) 
    sce.filtered 
  })
  names(sce.minCellCounts)<-tc.names
  saveRDS(sce.minCellCounts, file='sce_minCellCounts.rds')
}


