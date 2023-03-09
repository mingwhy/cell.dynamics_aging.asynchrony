
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

#########################################################
## subsample cell number to be equal in both age gropus
## read in data (a "SingleCellExperiment" object)
if(!file.exists('select.tc_subsampled.h5ad')){
  sce=readH5AD('select.tc.h5ad') 
  sce=sce[,sce$age %in% c('3m','24m')]  
  
  (pick.cell.types=as.character(unique(sce$tissue_cell.type))) 
  
  ## sub-sample cell: min(3,18,24), if 18m contain<50cell, discard this age group
  cell.meta=colData(sce)
  as.data.frame(cell.meta) %>% group_by(tissue_cell.type,age,sex) %>% summarise(n=n())
  
  min.ncell=50
  (tcs=unique(sce$tissue_cell.type))
  sce.subs<-lapply(tcs,function(tc){
    sce.tmp=sce[,sce$tissue_cell.type==tc]
    meta.tmp=colData(sce.tmp)
    x=names(which(table(meta.tmp$age)>=min.ncell))
    sce.tmp=sce.tmp[,sce.tmp$age %in% x]
    sce.tmp$age=factor(sce.tmp$age,levels=x)
    meta.tmp=colData(sce.tmp)
    n=min(table(meta.tmp$age))
    cat(as.character(tc),x,n,'\n')
    tmp=lapply(levels(sce.tmp$age),function(age.i){
      x=sce.tmp[,sce.tmp$age==age.i]
      set.seed(1224)
      i=sample(1:ncol(x),n,replace = F)
      x[,i]
    })
    sce.sub<-Reduce(`cbind`,tmp)
    sce.sub
  })
  lapply(sce.subs,dim)
  sce<-Reduce(cbind,sce.subs)
  assayNames(sce)[1] <- "counts" 
  sce #22966 17442
  writeH5AD(sce, 'select.tc_subsampled.h5ad')
}

sce=readH5AD('select.tc_subsampled.h5ad')
(pick.cell.types=as.character(unique(sce$tissue_cell.type))) #72 tc
table(sce$age,sce$tissue_cell.type) #same #cell for each tc across age groups

#################################################################################
## remove gene which didn't expr at all 
sce.filter<-lapply(pick.cell.types,function(tc){
    tmp=sce[,sce$tissue_cell.type==tc]
    mat=assay(tmp,'counts')
    #i=which(Matrix::rowSums(mat>0)> ncol(mat)*0.1)
    i=which(Matrix::rowSums(mat>0)!=0)
    tmp[i,]
  })
names(sce.filter)<-pick.cell.types;


## subsample cell read count
library(countland) #https://github.com/shchurch/countland/blob/main/tutorials_and_vignettes/R_tutorials_and_vignettes/vignette-tutorial.Rmd

set.seed(84095) # set random seed for reproducibility
minCellCounts<-lapply(sce.filter,function(tmp){
  m=assay(tmp,'counts')
  
  C <- countland(m)
  C <- Subsample(C,cell_counts='min')
  # preview the resulting subsampled count matrix
  #C@subsample[1:10,1:10]
  return(C)
})

Matrix::colSums(minCellCounts[[1]]@subsample)
minCellCounts[[1]]@subsample[1:3,1:3]
names(minCellCounts)<-names(sce.filter)

# create single cell object
tc.names=names(minCellCounts)
#tc.names=tc.names[tc.names!="Thymus:thymocyte"] #38tc

## filter gene (expr in at least 10% cells per cell type) and filter cell (expr >=100 genes)
sce.minCellCounts<-lapply(tc.names,function(x){
  tmp=minCellCounts[[x]]@subsample
  
  # filter gene
  n.expr.cell = Matrix::rowSums(tmp>0)
  tmp=tmp[n.expr.cell>= 0.1*ncol(tmp),]
  
  # filter cell
  n.expr.gene=Matrix::colSums(tmp>0)
  #cat(x);print(summary(n.expr.gene))
  mat=tmp[,n.expr.gene>=100] #expr at least 100 genes
  
  cell.info=colData(sce.filter[[x]]) #modify respective cell.meta data
  cell.info=cell.info[n.expr.gene>=100,]
  cat(x);print(table(cell.info$age)) #all contain >=50cells per age per tc
  
  #sce <- SingleCellExperiment(list(counts=minCellCounts[[x]]@subsample),colData=colData(sce.filter[[x]]) )
  sce <- SingleCellExperiment(list(counts=mat),colData=cell.info) 
  sce 
})
names(sce.minCellCounts)<-tc.names
saveRDS(sce.minCellCounts, file='sce_minCellCounts.rds')




