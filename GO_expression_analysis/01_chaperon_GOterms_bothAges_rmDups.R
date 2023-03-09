
library(org.Mm.eg.db,verbose=F,quietly=T)
library(GO.db);
options(stringsAsFactors = F)
library(igraph)
library(Matrix);
library(ggplot2);library(gridExtra)
library(ggpubr)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
library(UpSetR)

####################################
if(!file.exists('GOid_GOterm.rds')){
  slim.go2fb=readRDS('Mmus_go2SYMBOL_BP.rds')
  length(slim.go2fb) #1882 or 12545 (all BP GOterms)
  all.go.id=names(slim.go2fb)
  goterms=Term(GOTERM)
  x=lapply(all.go.id,function(i) GOTERM[[i]]@Term)
  df.go=data.frame(GOid=all.go.id,GOterm=unlist(x))
  dim(df.go) #12545
  saveRDS(df.go,'GOid_GOterm.rds')
}

df.go=readRDS('GOid_GOterm.rds')
df.go.chaperone=df.go[grep('chaperone',ignore.case = t,df.go$GOterm),]
rownames(df.go.chaperone)=df.go.chaperone$GOid
dim(df.go.chaperone) #14 x 2

slim.go2fb=readRDS('Mmus_go2SYMBOL_BP.rds')
go2gene=slim.go2fb[names(slim.go2fb) %in% df.go.chaperone$GOid]
x=sort(sapply(go2gene,nrow))
#GO:0051086 GO:0090035 GO:1903644 GO:1904714 GO:1904716 GO:1904764 GO:1903645 GO:1903646 GO:1904715 GO:0061740 
#1          1          1          1          1          1          2          2          2          3 
#GO:0061684 GO:0051131 GO:0061077 GO:0051085 
#7         18         25         33 
keep.go=names(x[x>5])
df.go.chaperone=df.go.chaperone[keep.go,]
go2gene=go2gene[keep.go]

#https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
input.list=lapply(go2gene,'[[',2)
names(input.list)
input.list1=input.list
sum(names(input.list1)==df.go.chaperone$GOid)
names(input.list1)=paste0(df.go.chaperone$GOid,'\n',df.go.chaperone$GOterm)

p1=upset(fromList(input.list1),nsets = length(input.list1), order.by = "freq")

pdf('upset.plot_chaperon_GOterm.pdf',useDingbats = T)
print(p1)
dev.off()

input.list2=input.list[-4]
upset(fromList(input.list2), order.by = "freq")

####################################
library(zellkonverter)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(ggpointdensity)
library(viridis)
library(patchwork); #for plot_annotation
library(lme4)
library(variancePartition)
library(ggpubr)
options(stringsAsFactors = F)
library(org.Mm.eg.db,verbose=F,quietly=T)
library(GO.db);
####################################
## reading in gene id.mapping and extract MHC genes
id.mapping=data.table::fread('fac_20449genes_id.mapping.txt')  

####################################
## read in mouse turnover rate data
source('../TMS_transcriptome_variation/src_cellular_lifespan.R')
dim(cell.lifespan)

cell.lifespan$`tissue: cell.type in mouse`=cell.lifespan$`cell type annotation in TMS`
cell.lifespan$human_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`

tc.orders=cell.lifespan[order(cell.lifespan$lifespan),]$`tissue: cell.type in mouse`
tc.orders #37 tc

########################################################################################
## use expression proportion (remember to change `expr.m=expr.m[,tmp$age==___ ]`)
pick.age='3m';
#pick.age='24m';
output_file=paste0('mouse_male_GO_',pick.age,'_libNorm_Chaperon.rds');
output_file

if(!file.exists(output_file)){
  if(T){
    sce=readH5AD('select.tc.h5ad')
    assayNames(sce)<-'counts'
    sce_naive=sce[,sce$age %in% c('3m','24m')]    
    assayNames(sce_naive)<-'counts'
    sce_naive<-logNormCounts(sce_naive, log=FALSE, pseudo.count=1)     
  }
    
  tc.names=tc.orders
  
  mouse_tcs_TDI.list<-lapply(tc.names,function(tc){
    
    tmp=sce_naive[,sce_naive$tissue_cell.type==tc] #raw count data
    
    expr.m=assay(tmp,'normcounts')
    expr.m=expr.m[,tmp$age==pick.age]
    
    n.expr.gene=Matrix::colSums(expr.m>0) 
    expr.m=expr.m[, n.expr.gene>=100] 
       
    cell.expr=Matrix::colSums(expr.m)
    all.genes=rownames(expr.m)
    
    tmp.go=names(input.list);
    go.mean.expr<-lapply(tmp.go, function(go){
      overlap.genes=intersect(all.genes,input.list[[go]])
      if(length(overlap.genes)==0){return(c(0,NA))}
      expr.m.tmp=expr.m[overlap.genes,,drop=F]
      go.expr=Matrix::colSums(expr.m.tmp)
      return(c(length(overlap.genes),mean(go.expr/cell.expr)))  
    })
    
    return(as.data.frame(Reduce(`rbind`,go.mean.expr)))
  })
  tmp.go=names(input.list);
  mouse_tcs_TDI.list2<-lapply(mouse_tcs_TDI.list,function(i){
    colnames(i)=c('ngene','score');
    i$GO.id=tmp.go;
    i})
  mouse_tcs_TDI=as.data.frame(Reduce(`rbind`,mouse_tcs_TDI.list2))
  mouse_tcs_TDI$cell.type=rep(tc.names, sapply(mouse_tcs_TDI.list2,nrow) )
  
  saveRDS(mouse_tcs_TDI,output_file)
}


