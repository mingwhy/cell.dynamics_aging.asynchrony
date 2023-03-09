
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)
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

n.cell.per.tc<-lapply(tc.names,function(tc){
  sce_naive=sce.minCellCounts[[tc]]
  assayNames(sce_naive)
  expr.m=assay(sce_naive,'counts')
  
  out=lapply(unique(sce_naive$age),function(i){
    m=expr.m[,sce_naive$age==i]
    m=m[Matrix::rowSums(m)!=0,]
    ncol(m)
  })
  return(unlist(out))
})
df.n.cell.per.tc=Reduce(`rbind`,n.cell.per.tc)
rownames(df.n.cell.per.tc)=tc.names
df.n.cell.per.tc=df.n.cell.per.tc[order(df.n.cell.per.tc[,1]),]
head(df.n.cell.per.tc)
tail(df.n.cell.per.tc)
summary(df.n.cell.per.tc[,1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#55.0    78.5   119.0   222.2   239.5  2253.0 

###########################################################
# spearman's cor coeff on the down-sampled expr.mat across all genes 
# between all pairwise cell comparisons per cell type per age group
if(!file.exists('res_cell.cors.rds')){
  DM.out<-lapply(tc.names,function(tc){
    sce_naive=sce.minCellCounts[[tc]]
    assayNames(sce_naive)
    expr.m=assay(sce_naive,'counts')
  
    out=lapply(unique(sce_naive$age),function(i){
      m=expr.m[,sce_naive$age==i]
      m=m[Matrix::rowSums(m)!=0,]
     
      # take some time to compute
      cell.cors=cor(as.matrix(m),method='spearman')
      cell.cors=cell.cors[upper.tri(cell.cors,diag=FALSE)]
      df=data.frame('cell.type'=tc,'cell.cors'=cell.cors,'age'=i)
      return(df)
    })
    out=as.data.frame(Reduce(`rbind`,out))
    return(out)
  })
  
  DM.out=as.data.frame(Reduce(`rbind`,DM.out))
  saveRDS(DM.out,'res_cell.cors.rds')
}


###########################################################
## plot per cell.type: young vs old
res_dist=readRDS('res_cell.cors.rds')

res_dist$cell.type=factor(res_dist$cell.type,
                          levels=sort(unique(res_dist$cell.type)))
res_dist$cell.dist=1-res_dist$cell.cors; 
res_dist$age=factor(res_dist$age,levels=c('3m','24m'))

res_dist_tc<-res_dist[res_dist$cell.type %in% tc.names,]

tmp=res_dist_tc %>% group_by(cell.type,age) %>% summarise(age.median=median(cell.dist))
tmp1 = tmp %>% spread(age,age.median)
sum(is.na(tmp1$`24m`))
#tmp1=tmp1[!is.na(tmp1$`24m`),]
tmp1$log2ratio=log2(tmp1$`24m`/tmp1$`3m`)

plots <- lapply(tc.names,function(i){
  df=res_dist_tc[res_dist_tc$cell.type==i,]
  ratio.value=tmp1[tmp1$cell.type==i,]$log2ratio
  ggplot(df,aes(x=age,y=cell.dist))+
    geom_violin()+#geom_jitter(size=0.1)+
    stat_summary(fun.y=median, geom="point", size=2, color="red")+
    theme_classic()+ylab('1-Spearman\'s rho' )+
    ggtitle(paste0(i,'\nlog2(O/Y)=',round(ratio.value,5)))+
    theme(plot.title = element_text(size=9))
})
length(plots)



pdf('res_1-spearman_dist.pdf',useDingbats = T,width = 12,height = 17)
grid.arrange(grobs=plots,ncol=4)
dev.off()

