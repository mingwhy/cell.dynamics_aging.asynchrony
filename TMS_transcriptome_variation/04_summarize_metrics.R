
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)

used.stat='cv2'; 
source('src_cellular_lifespan.R')
tc.orders=cell.lifespan[order(cell.lifespan$lifespan),]$`cell type annotation in TMS`
tc.orders

#############################################################
## gene-level metric
res.gene=readRDS('LMM_cv2_out.rds')
dim(res.gene) 
head(res.gene)
colnames(res.gene)[1]='cell.type'

res.gene=res.gene[res.gene$cell.type %in% tc.orders,]
res.gene=res.gene[order(res.gene$cell.type),]

head(res.gene)
summary(p.adjust(res.gene$pvalue,method='BH')) 

max(p.adjust(res.gene$pvalue,method='bonferroni')) #1
res.gene$FDR=p.adjust(res.gene$pvalue,method='BH')
sum(res.gene$FDR<0.05 & res.gene$beta>0) #35 tc
sum(res.gene$FDR<0.05 & res.gene$beta<0) #3 tc
res.gene[res.gene$FDR>=0.05,] #only 1 tc, SCAT:B cell

############################################################
## cell: res_cell.cors.rds
df=readRDS('res_cell.cors.rds')
df=df[df$age %in% c('3m','24m'),]
df$age=droplevels(df$age)
df$dist=1-df$cell.cors;

df.pval=df %>% group_by(cell.type) %>%  
  summarize(pval = wilcox.test(dist ~ age)$p.value)
df.pval$p.adj=p.adjust(df.pval$pval,method='BH')

tmp=df %>% group_by(cell.type,age) %>% summarise(age.median=median(dist))
tmp1 = tmp %>% spread(age,age.median)
tmp1=tmp1[!is.na(tmp1$`24m`),]
tmp1$log2ratio=log2(tmp1$`24m`/tmp1$`3m`)
res.spearman=tmp1;

res.spearman=res.spearman[res.spearman$cell.type %in% tc.orders,]
unique(res.spearman$cell.type) #39 tc
res.spearman=res.spearman[order(res.spearman$cell.type),]

res.spearman=merge(res.spearman,df.pval)
sum(res.spearman$log2ratio>0) #26
sum(res.spearman$log2ratio<0) #13
sum(res.spearman$log2ratio>0 & res.spearman$p.adj<0.05) #26 
summary(res.spearman[res.spearman$log2ratio>0 & res.spearman$p.adj<0.05,]$p.adj)

############################################################
## cell: GCL
df=data.table::fread('sce.minCellCounts_gcl.txt')
df$cell.type=df$tc
df$dist=1-df$gcl;
df  %>% group_by(cell.type,age) %>% summarize(n()) #100 reps per tc

df.pval=df %>% group_by(cell.type) %>%  
  summarize(pval = wilcox.test(dist ~ age)$p.value)
df.pval$p.adj=p.adjust(df.pval$pval,method='BH')

tmp=df %>% group_by(cell.type,age) %>% summarise(age.median=median(dist))
tmp1 = tmp %>% spread(age,age.median)
tmp1=tmp1[!is.na(tmp1$`24m`),]
tmp1$log2ratio=log2(tmp1$`24m`/tmp1$`3m`)
res.gcl=tmp1;

res.gcl=res.gcl[res.gcl$cell.type %in% tc.orders,]
unique(res.gcl$cell.type) 
res.gcl=res.gcl[order(res.gcl$cell.type),]

res.gcl=merge(res.gcl,df.pval)
sum(res.gcl$log2ratio>0) #34
sum(res.gcl$log2ratio<0) #5
sum(res.gcl$log2ratio>0 & res.gcl$p.adj<0.05) #34
summary(res.gcl[res.gcl$log2ratio>1 & res.gcl$p.adj<0.05,]$p.adj)

#######################################################################
## make sure cell.types are all the same order among different metrics
dim(res.spearman);sum(res.spearman$cell.type==res.gcl$cell.type)
sum(res.spearman$cell.type==res.gene$cell.type);

df=data.frame(res.spearman$log2ratio,res.gcl$log2ratio,res.gene$beta)
if(used.stat=='cv2'){
  colnames(df)=c('1-spearman','1-GCL','gene_CV2')
}

df$cell.type=res.gene$cell.type
data.table::fwrite(df,'age.change_all.metrics.txt')


