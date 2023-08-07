
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)
library(ggpmisc) #https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
#https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
scaleFUN <- function(x) sprintf("%.3f", x)
mycols=c("#999999", "#E69F00");
#mycols=c("#999999",  "#56B4E9","#E69F00");

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

n.cell.per.tc<-lapply(tc.names,function(tc){
  sce_naive=sce.minCellCounts[[tc]]
  assayNames(sce_naive)
  expr.m=assay(sce_naive,'counts')
  cat(tc,'\n')
  
  df=colData(sce_naive)
  df$mouse.id_age=paste(df$mouse.id,df$age)
  mouse.id.names=names(which(table(df$mouse.id_age)>=min.ncell_per.mouse))
  table(df[df$mouse.id_age %in% mouse.id.names,]$mouse.id_age)
  ncell=min(table(df[df$mouse.id_age %in% mouse.id.names,]$mouse.id_age))
  
  #out=lapply(unique(sce_naive$age),function(i){
  out=lapply(mouse.id.names,function(i){
    m=expr.m[,df$mouse.id_age==i]
    m=m[Matrix::rowSums(m)!=0,,drop=FALSE]
    m=m[,Matrix::colSums(m>0)>=100,drop=FALSE]
    if(ncol(m)<min.ncell_per.mouse){return(NULL)}
    ncol(m)
  })
  return(unlist(out))
})
summary(unlist(n.cell.per.tc))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#15.0    31.0    48.0   104.2   119.0  1215.0

###########################################################
# spearman's cor coeff on the down-sampled expr.mat across all genes 
# between all pairwise cell comparisons per cell type per age group
if(!file.exists('res_cell.cors_per.mouse.rds')){
  
  DM.out<-lapply(tc.names,function(tc){
    sce_naive=sce.minCellCounts[[tc]]
    assayNames(sce_naive)
    
    df=colData(sce_naive)
    df$mouse.id_age=paste(df$mouse.id,df$age)
    
    # sample equal number of cells
    mouse.id.names=names(which(table(df$mouse.id_age)>=min.ncell_per.mouse))
    table(df[df$mouse.id_age %in% mouse.id.names,]$mouse.id_age)
    ncell=min(table(df[df$mouse.id_age %in% mouse.id.names,]$mouse.id_age))
    
    out=lapply(mouse.id.names,function(i){
      sce.tc=sce_naive[,df$mouse.id_age==i]
      if(ncol(sce.tc)<min.ncell_per.mouse){return(NULL)}
      
      tmp=sce.tc;
      m=assay(tmp,'counts')
      m=m[Matrix::rowSums(m)!=0,]
      m=m[,Matrix::colSums(m>0)>=100,drop=FALSE] #expr >=200 genes per cell 
      dim(m)
      if(ncol(m)<min.ncell_per.mouse){return(NULL)}
      
      # take some time to compute
      cell.cors=cor(as.matrix(m),method='spearman')
      #cell.cors=cor(as.matrix(m2),method='spearman')
      cell.cors=cell.cors[upper.tri(cell.cors,diag=FALSE)]
      df=data.frame('cell.type'=tc,'cell.cors'=cell.cors,'age'=tmp$age[1],'mouse.id'=tmp$mouse.id[1])
      return(df)
    })
    
    out2<-Filter(Negate(is.null), out)
    df.out2=as.data.frame(Reduce(`rbind`,out2))
    return(df.out2)
  })
  
  DM.out=as.data.frame(Reduce(`rbind`,DM.out))
  saveRDS(DM.out,'res_cell.cors_per.mouse.rds')
}

###########################################################
## plot per cell.type: young vs old
res_dist=readRDS('res_cell.cors_per.mouse.rds')

res_dist$cell.type=factor(res_dist$cell.type,
                          levels=sort(unique(res_dist$cell.type)))
res_dist$cell.dist=1-res_dist$cell.cors; #dist=1-spearman.cor.coeff, apply to both positive and negative values 

res_dist$age=factor(res_dist$age,levels=age.groups)

res_dist_tc<-res_dist[res_dist$cell.type %in% tc.names,]

tmp0=res_dist_tc %>% group_by(cell.type,age,mouse.id) %>% summarise(age.average.per.mouse=median(cell.dist))
tmp=tmp0 %>% group_by(cell.type,age) %>% summarise(age.average=mean(age.average.per.mouse))

tmp1 = tmp %>% spread(age,age.average)
nrow(tmp1)
sum(is.na(tmp1[[age.groups[[1]]]]))
sum(is.na( tmp1[[age.groups[[2]]]]))

tmp1$log2ratio=log2(tmp1[[age.groups[[2]]]]/ tmp1[[age.groups[[1]]]])
tmp1=tmp1[!is.na(tmp1$log2ratio),]
nrow(tmp1) 
sum(tmp1$log2ratio>0) 

plots <- lapply(unique(tmp1$cell.type),function(i){
  df=res_dist_tc[res_dist_tc$cell.type==i,]
  df=df[order(df$age),]
  df$mouse.id=factor(df$mouse.id,unique(df$mouse.id))
  ratio.value=tmp1[tmp1$cell.type==i,]$log2ratio
  ggplot(df,aes(x=mouse.id,y=cell.dist,col=age,fill=age))+
    geom_violin()+#geom_jitter(size=0.1)+
    stat_summary(fun.y=median, geom="point", size=1, color="black")+
    theme_classic()+ylab(expression(paste('1-Spearman\'s ',rho)))+
    ggtitle(paste0(i,'\nlog2(O/Y)=',sprintf('%.3f',ratio.value)))+
    scale_color_manual(name='Month',values=mycols)+
    scale_fill_manual(name='Month',values=mycols)+
    theme(plot.title = element_text(size=9),legend.position = 'none',
          axis.text.x = element_text(size=6))+
    xlab('')+ylab('')
})
length(plots)

pdf('res_1-spearman_dist_per.mouse.pdf',useDingbats = T,width = 12,height = 17)
grid.arrange(grobs=plots,ncol=4)
dev.off()

