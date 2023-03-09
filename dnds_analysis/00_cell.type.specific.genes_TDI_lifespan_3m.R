
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

####################################
## read in mouse turnover rate data
source('../TMS_transcriptome_variation/src_cellular_lifespan.R')
dim(cell.lifespan)

cell.lifespan$`tissue: cell.type in mouse`=cell.lifespan$`cell type annotation in TMS`
cell.lifespan$human_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`

tc.orders=cell.lifespan[order(cell.lifespan$lifespan),]$`tissue: cell.type in mouse`
tc.orders #37 tc

###############################################################################################
## reading in gene id.mapping and dn.ds infomation
id.mapping=data.table::fread('fac_20449genes_id.mapping.txt')  
filter_mouse_rat=data.table::fread('mouse_rat.dnds.txt')
head(filter_mouse_rat)
dim(filter_mouse_rat) #30033 genes
dim(id.mapping) #20449
sum(filter_mouse_rat$ensembl_gene_id %in% id.mapping$ensembl_gene_id) #20697

x=filter_mouse_rat[filter_mouse_rat$rnorvegicus_homolog_orthology_type=='ortholog_one2one',]
anyDuplicated(x$ensembl_gene_id)
#mouse_rat_dnds=filter_mouse_rat[!duplicated(filter_mouse_rat$ensembl_gene_id),] #random choose one ortholog in rat
mouse_rat_dnds=x; #16494

head(mouse_rat_dnds)
head(id.mapping)
gene.meta=merge(mouse_rat_dnds,id.mapping,by.x='ensembl_gene_id',by.y='ensembl_gene_id')
dim(gene.meta) #14821
gene.meta$omega=gene.meta$rnorvegicus_homolog_dn/gene.meta$rnorvegicus_homolog_ds
gene.meta[is.infinite(gene.meta$omega),]$omega=NA
gene.meta[is.nan(gene.meta$omega),]$omega=NA
summary(gene.meta$omega)  #max=2.24 
gene.meta[which(gene.meta$omega>=1),] #75 genes with omega>1

#remove genes with dn/ds>1 (more likely to be under positive selection)
length(which(gene.meta$omega<=1)) #14110 genes
gene.meta=gene.meta[!is.na(gene.meta$omega),]
sum(gene.meta$omega>=1) #76 genes

gene.meta=gene.meta[which(gene.meta$omega<1),]
dim(gene.meta) #14185 or 14109    12

summary(replicate(100,median(gene.meta[sample(1:nrow(gene.meta),500,replace = F),]$omega)))
summary(replicate(100,median(gene.meta[sample(1:nrow(gene.meta),3500,replace = F),]$omega)))


######################################################################
## read in TMS data

if(T){
  sce=readH5AD('../TMS_transcriptome_variation/select.tc.h5ad') # 22966 31001 
  sce_naive=sce[,sce$age=='3m']
  sce_naive #22966 10763 
  
  assayNames(sce_naive)<-'counts'
  sce_naive<-logNormCounts(sce_naive, log=FALSE, pseudo.count=1) #if log, the cell.lib.size range 2~3 orders
  assayNames(sce_naive)
  
}


## use average gene expr per cell type to select cell type-specific genes
if(!file.exists('df.expr.per.tc.rds')){
  expr.per.tc<-lapply(cell.lifespan$`tissue: cell.type in mouse`,function(tc){
    tmp=sce_naive[,sce_naive$tissue_cell.type==tc] #raw count data
    if(F){
      counts=assay(sce_naive,'X')
      libsizes <- colSums(counts)
      size.factors <- libsizes/mean(libsizes)
      logcounts=  log2(t(t(counts)/size.factors) + 1)
      Matrix::rowMeans(logcounts)
    }
    
    expr.m=assay(tmp,'normcounts')
    expr.m=expr.m[,tmp$age=='3m']
    Matrix::rowMeans(expr.m)
  })
  df.expr.per.tc=as.data.frame(Reduce(`cbind`,expr.per.tc))
  colnames(df.expr.per.tc)=cell.lifespan$`tissue: cell.type in mouse`
  sum(Matrix::rowSums(df.expr.per.tc)==0) #534, remove non-expressed genes
  df.expr.per.tc=df.expr.per.tc[Matrix::rowSums(df.expr.per.tc)!=0,]
  dim(df.expr.per.tc) 
  hist(as.numeric(df.expr.per.tc[1,]))
  saveRDS(as.data.frame(df.expr.per.tc), 'df.expr.per.tc.rds')
}

##########################################################################
## look at TDI using only tissue-specific genes 
#Understanding Tissue-Specific Gene Regulation, https://www.sciencedirect.com/science/article/pii/S2211124717314183?via%3Dihub
df.expr.per.tc=readRDS('df.expr.per.tc.rds')
df.expr.per.tc[1:3,1:3]
x=as.numeric(df.expr.per.tc[1,])
(x-median(x))/ as.numeric(quantile(x,c(0.75))-quantile(x,c(0.25)))
zscore=apply(df.expr.per.tc,1,function(x){
  x=as.numeric(x)
  (x-median(x))/ as.numeric(quantile(x,c(0.75))-quantile(x,c(0.25)))
})
dim(zscore) #cell type by tissue 

cut.offs=c(0,0.5,1,1.5, 2,2.5, 3,3.5); #larger than 3.5, too few express gene per cell
tmp=sapply(cut.offs,function(cut.off){apply(zscore,1,function(j){sum(j>cut.off,na.rm=T)})})
colnames(tmp)=as.character(cut.offs)
rownames(tmp)=colnames(df.expr.per.tc)

for(cut.off.value in cut.offs){

  output.file=paste0('mouse_male_binary_TDI_ts_',cut.off.value,'.rds')
  #output.file=paste0('mouse_male_binary_TDI_ts_',cut.off.value,'_filter.dnds1.rds')
  
  if(file.exists(output.file)){next}
  ts.genes<-apply(zscore,1,function(i){
    names(which(i>cut.off.value))
  })
  names(ts.genes)<-colnames(df.expr.per.tc)
  tmp=data.frame(cell.type=names(ts.genes),ngene=sapply(ts.genes,length))
  summary(tmp$ngene)
  
  ## calculate dnds using specific genes per cell type
  tc.names=cell.lifespan$`tissue: cell.type in mouse`
  mouse_tcs_TDI.list<-lapply(tc.names,function(tc){
    sce_naive_one=sce_naive[,sce_naive$tissue_cell.type==tc] #raw count data
    
    assayNames(sce_naive_one)
    expr.m=assay(sce_naive_one,'normcounts')
    
    overlap.genes=intersect(rownames(expr.m),gene.meta$mgi_symbol)
    overlap.genes=intersect(ts.genes[[tc]],overlap.genes)
    cat(tc,nrow(expr.m),length(overlap.genes),'\n');
    
    expr.m=expr.m[overlap.genes,]
    n.expr.gene=Matrix::colSums(expr.m>0) #control for the number of expressed genes （https://academic.oup.com/gbe/article/12/4/300/5807614?login=true
    
    i=match(overlap.genes,gene.meta$mgi_symbol)
    gene.meta.m=gene.meta[i,]
    dim(gene.meta.m);dim(expr.m)
    sum(gene.meta.m$mgi_symbol==rownames(expr.m))
    
    expr.m.binary=expr.m;
    expr.m.binary[expr.m.binary>0]=1;
    x=gene.meta.m$omega %*% as.matrix(expr.m.binary)
    index.per.cell=x/n.expr.gene;
    
    length(index.per.cell)
    cell.meta=colData(sce_naive_one)
    dim(cell.meta)
    cell.meta$TDI=index.per.cell;
    cell.meta$n_expr_gene=n.expr.gene
    cell.meta=as.data.frame(cell.meta)
    return(cell.meta)
  })
  
  mouse_tcs_TDI=as.data.frame(Reduce(`rbind`,mouse_tcs_TDI.list))
  head(mouse_tcs_TDI)
  summary(mouse_tcs_TDI$sizeFactor) #as we use subsampled UMI
  summary(mouse_tcs_TDI$n_expr_gene) #44 to 4482
  
  saveRDS(mouse_tcs_TDI,output.file)
}
 

