
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

pick.age='3m';
#pick.age='24m';
######################################################################
## read in mouse turnover rate data
cell.lifespan=data.table::fread('../../src/Dataset_S1.txt')
tc.orders=cell.lifespan[order(cell.lifespan$lifespan),]$`cell type annotation in TMS`
tc.orders #39 tc

###############################################################################################
## reading in gene id.mapping and dn.ds infomation
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt')  #readin_h5ad.R
filter_mouse_rat=data.table::fread('../../src/mouse_rat.dnds.txt')
head(filter_mouse_rat)
dim(filter_mouse_rat) #30033 genes
dim(id.mapping) #20449
sum(filter_mouse_rat$ensembl_gene_id %in% id.mapping$ensembl_gene_id) #20697

tmp=filter_mouse_rat$rnorvegicus_homolog_dn/filter_mouse_rat$rnorvegicus_homolog_ds
tmp[!is.infinite(tmp) & !is.nan(tmp)]->tmp
summary(tmp)

x=filter_mouse_rat[filter_mouse_rat$rnorvegicus_homolog_orthology_type=='ortholog_one2one',]
anyDuplicated(x$ensembl_gene_id)
mouse_rat_dnds=x; 
nrow(mouse_rat_dnds); #16494

head(mouse_rat_dnds)
head(id.mapping)
gene.meta=merge(mouse_rat_dnds,id.mapping,by.x='ensembl_gene_id',by.y='ensembl_gene_id')
dim(gene.meta) #14821
gene.meta$omega=gene.meta$rnorvegicus_homolog_dn/gene.meta$rnorvegicus_homolog_ds
gene.meta[is.infinite(gene.meta$omega),]$omega=NA
gene.meta[is.nan(gene.meta$omega),]$omega=NA
gene.meta=gene.meta[!is.na(gene.meta$omega),]
dim(gene.meta) #14185    12

summary(gene.meta$omega)  #max=2.248201
nrow(gene.meta[which(gene.meta$omega>=1),]) #76 genes with omega>=1

#remove genes with dn/ds>1 (more likely to be under positive selection)
nrow(gene.meta) #14185
length(which(gene.meta$omega<=1)) #14110 genes

sum(gene.meta$omega>=1) #76 genes
#gene.meta=gene.meta[which(gene.meta$omega<1),] #remove or not
dim(gene.meta) #14185 or 14109 genes

summary(gene.meta$omega)
summary(replicate(100,median(gene.meta[sample(1:nrow(gene.meta),500,replace = F),]$omega)))
summary(replicate(100,median(gene.meta[sample(1:nrow(gene.meta),3500,replace = F),]$omega)))
#sum(gene.meta$rnorvegicus_homolog_orthology_confidence)
#gene.meta=gene.meta[gene.meta$rnorvegicus_homolog_orthology_confidence==1,]

######################################################################
## read in TMS data

if(T){
  sce=readH5AD('../../0708_TMS_male_3m_24m/select.tc.h5ad') # 22966 31001 
  sce_naive=sce[,sce$age==pick.age]
  sce_naive 
  
  assayNames(sce_naive)<-'counts'
  sce_naive<-logNormCounts(sce_naive, log=FALSE, pseudo.count=1) #if log, the cell.lib.size range 2~3 orders
  assayNames(sce_naive)
  
  #tmp=assay(sce_naive,'normcounts')
  #unique(colSums(tmp))
}

## use average gene expr per cell type to select cell type-specific genes
ts.file=paste0('df.expr.per.tc_',pick.age,'.rds')
if(!file.exists(ts.file)){
  expr.per.tc<-lapply(cell.lifespan$`cell type annotation in TMS`,function(tc){
    tmp=sce_naive[,sce_naive$tissue_cell.type==tc] #raw count data
    expr.m=assay(tmp,'normcounts')
    expr.m=expr.m[,tmp$age==pick.age]
    Matrix::rowMeans(expr.m)
  })
  df.expr.per.tc=as.data.frame(Reduce(`cbind`,expr.per.tc))
  colnames(df.expr.per.tc)=cell.lifespan$`cell type annotation in TMS`
  sum(Matrix::rowSums(df.expr.per.tc)==0) #534, remove non-expressed genes
  df.expr.per.tc=df.expr.per.tc[Matrix::rowSums(df.expr.per.tc)!=0,]
  dim(df.expr.per.tc) 
  hist(as.numeric(df.expr.per.tc[1,]))
  saveRDS(as.data.frame(df.expr.per.tc), ts.file)
}

##########################################################################
## look at TDI using only tissue-specific genes 
#Understanding Tissue-Specific Gene Regulation, https://www.sciencedirect.com/science/article/pii/S2211124717314183?via%3Dihub
#df.expr.per.tc=readRDS(ts.file)
df.expr.per.tc=readRDS('./all.genes/df.expr.per.tc_3m.rds')
df.expr.per.tc[1:3,1:3]
x=as.numeric(df.expr.per.tc[1,])
(x-median(x))/ as.numeric(quantile(x,c(0.75))-quantile(x,c(0.25)))
zscore=apply(df.expr.per.tc,1,function(x){
  x=as.numeric(x)
  (x-median(x))/ as.numeric(quantile(x,c(0.75))-quantile(x,c(0.25)))
})
dim(zscore) #cell type by tissue 

#cut.offs=c(0,1,1.5, 2,2.5, 3,3.5, 4,4.5, 5);
cut.offs=c(0,0.5,1,1.5, 2,2.5, 3,3.5); #larger than 3.5, too few express gene per cell
tmp=sapply(cut.offs,function(cut.off){apply(zscore,1,function(j){sum(j>cut.off,na.rm=T)})})
colnames(tmp)=as.character(cut.offs)
rownames(tmp)=colnames(df.expr.per.tc)
#summary(tmp[,'2']) #942 ~ 5307 genes

for(cut.off.value in cut.offs){
  #output.file=paste0('mouse_male_binary_TDI_ts_',cut.off.value,"_",pick.age,'_filter.dnds1.rds')
  output.file=paste0('mouse_male_binary_TDI_ts_',cut.off.value,"_",pick.age,'.rds')
  if(file.exists(output.file)){next}
  ts.genes<-apply(zscore,1,function(i){
    names(which(i>cut.off.value))
  })
  names(ts.genes)<-colnames(df.expr.per.tc)
  tmp=data.frame(cell.type=names(ts.genes),ngene=sapply(ts.genes,length))
  summary(tmp$ngene)
  
  ## calculate dnds using specific genes per cell type
  tc.names=cell.lifespan$`cell type annotation in TMS`
  mouse_tcs_TDI.list<-lapply(tc.names,function(tc){
    sce_naive_one=sce_naive[,sce_naive$tissue_cell.type==tc] #raw count data    
    assayNames(sce_naive_one)
    expr.m=assay(sce_naive_one,'normcounts')
    
    overlap.genes=intersect(rownames(expr.m),gene.meta$mgi_symbol)
    overlap.genes=intersect(ts.genes[[tc]],overlap.genes)
    cat(tc,nrow(expr.m),length(overlap.genes),'\n');
    
    expr.m=expr.m[overlap.genes,]
    n.expr.gene=Matrix::colSums(expr.m>0) 

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
 
 


