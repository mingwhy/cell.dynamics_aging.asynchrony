
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
  slim.go2fb=readRDS('../src/Mmus_go2SYMBOL_BP.rds')
  length(slim.go2fb) #1882 or 12545 (all BP GOterms)
  all.go.id=names(slim.go2fb)
  goterms=Term(GOTERM)
  x=lapply(all.go.id,function(i) GOTERM[[i]]@Term)
  df.go=data.frame(GOid=all.go.id,GOterm=unlist(x))
  dim(df.go) #12545
  saveRDS(df.go,'GOid_GOterm.rds')
}

df.go=readRDS('GOid_GOterm.rds')
df.go.pick=df.go[grep('DNA repair',ignore.case = t,df.go$GOterm),]
rownames(df.go.pick)=df.go.pick$GOid
dim(df.go.pick) #11 x 2

slim.go2fb=readRDS('../src/Mmus_go2SYMBOL_BP.rds')
go2gene=slim.go2fb[names(slim.go2fb) %in% df.go.pick$GOid]
x=sort(sapply(go2gene,nrow)) #11
x

keep.go=names(x[x>5]) # at least 6 gene members
length(keep.go)
df.go.pick=df.go.pick[keep.go,]
go2gene=go2gene[keep.go]
length(go2gene) #6

#https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
input.list=lapply(go2gene,'[[',2)
names(input.list)
input.list1=input.list
sum(names(input.list1)==df.go.pick$GOid)
names(input.list1)=paste0(df.go.pick$GOid,'\n',df.go.pick$GOterm)

p1=upset(fromList(input.list1),nsets = length(input.list1), order.by = "freq")

pdf('upset.plot_DNA.repair_GOterm.pdf',useDingbats = T)
print(p1)
dev.off()


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
## reading in gene id.mapping 
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt') 
gene.meta=id.mapping
####################################
## read in mouse turnover rate data
cell.lifespan=data.table::fread('../src/Dataset_S1.txt')
cell.lifespan$Ron_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`
tc.orders=cell.lifespan[order(cell.lifespan$lifespan),]$`cell type annotation in TMS`
tc.orders #39 tc

####################################
## use expression proportion
#pick.age='3m';
pick.age='24m';
if(!dir.exists('select.tc.h5ad_result_goterms')){dir.create('select.tc.h5ad_result_goterms')}
output_file=paste0('select.tc.h5ad_result_goterms/mouse_male_GO_',pick.age,'_libNorm2_DNArepair.rds');
output_file

if(!file.exists(output_file)){
  if(T){
    sce=readH5AD('../0708_TMS_male_3m_24m/select.tc.h5ad') 
    unique(sce$age) #3m 18m 24m
    assayNames(sce)<-'counts'
    sce_naive=sce[,sce$age %in% c('3m','24m')]
    table(sce_naive$age)
    assayNames(sce_naive)<-'counts'
    sce_naive<-logNormCounts(sce_naive, log=TRUE, pseudo.count=1.1) 
    assayNames(sce_naive)
  }
  
  
  tc.names=tc.orders
  
  mouse_tcs_TDI.list<-lapply(tc.names,function(tc){
    
    tmp=sce_naive[,sce_naive$tissue_cell.type==tc] #raw count data
    
    expr.m=assay(tmp,'logcounts')
    expr.m=expr.m[,tmp$age==pick.age]
    n.expr.gene=Matrix::colSums(expr.m>0) 
    expr.m=expr.m[, n.expr.gene>=100] #keep cells which expr>=100 genes
    
    cell.expr=Matrix::colSums(expr.m)
    all.genes=rownames(expr.m)
    
    tmp.go=names(input.list);
    go.median.expr<-lapply(tmp.go, function(go){
      overlap.genes=intersect(all.genes,input.list[[go]])
      if(length(overlap.genes)==0){return(c(0,NA))}
      
      expr.m.tmp=expr.m[overlap.genes,,drop=F]
      go.expr=Matrix::colSums(expr.m.tmp)
      return(c(length(overlap.genes),mean(go.expr/cell.expr)))   
    })
    
    return(as.data.frame(Reduce(`rbind`,go.median.expr)))
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

####################################
### see correlation
mouse_tcs_TDI=readRDS(output_file)
summary(mouse_tcs_TDI$ngene) #1~3006
length(unique(mouse_tcs_TDI$GO.id))

df.go.score<-mouse_tcs_TDI[,-1] %>% spread(GO.id,score)
dim(df.go.score); #cell type by go term
rownames(df.go.score)<-df.go.score$cell.type
df.go.score=df.go.score[,-1]
all.go=colnames(df.go.score)
dim(df.go.score); #39 x 6

overlap.tc=intersect(rownames(df.go.score),cell.lifespan$`cell type annotation in TMS`)

df1=df.go.score[overlap.tc,]
df2=cell.lifespan[match(overlap.tc,cell.lifespan$`cell type annotation in TMS`),]

df2$duplicate='tissue-specific estimate';
df2[df2$Ron_tc %in% df2$Ron_tc[duplicated(df2$Ron_tc)],]$duplicate='non tissue-specific estimate'
table(df2$duplicate)

na.number=apply(df1,2,function(i) sum(is.na(i)))
na.number

sum(rownames(df1)==df2$`cell type annotation in TMS`) #39
df1$`cell type annotation in TMS`=df2$`cell type annotation in TMS`

df1.long=reshape2::melt(df1)
colnames(df1.long)[2]='GOid'
intersect(colnames(df1.long),colnames(df2))
df3=merge(df1.long,df2)

#######################################################################################
## use average for non-tissue specific ones
x=df3 %>% group_by(Ron_tc,GOid) %>% summarise(average=mean(value))
df4=x %>% spread(GOid,average)
dim(df4) 

# at least 10 cell types terms express this GO term's gene members
tc.size=unlist(lapply(2:ncol(df4),function(i) sum(df4[,i]!=0) ))


intersect(colnames(df4),colnames(df2)) #Ron_tc
df5=merge(df4,df2[!duplicated(df2$Ron_tc),],all.x=TRUE)
dim(df5) 

Matrix=df4[,-1];
dim(Matrix) 

df5$turnover=1-exp(-1/df5$lifespan)
SampleAge=df5$lifespan

cor.coeffs.list=t(apply(Matrix,2,function(x){
  i=cor.test(x,SampleAge,method='kendall',use='pairwise.complete.obs')
  c(i$estimate,i$p.value)
}))
cor.coeffs=as.data.frame(cor.coeffs.list)
colnames(cor.coeffs)=c('Spearman.rho','Pval')
cor.coeffs$GOid=rownames(cor.coeffs)
cor.coeffs=cor.coeffs[!is.na(cor.coeffs$Pval),]
cor.coeffs$FDR=p.adjust(cor.coeffs$Pval,method='BH')
cor.coeffs$tc.size=tc.size
  
dim(cor.coeffs) 
cor.coeffs=merge(cor.coeffs,df.go.pick)
cor.coeffs=cor.coeffs[order(cor.coeffs$Pval),]
cor.coeffs

plot_go_id=cor.coeffs$GOid


plots<-lapply(1:nrow(cor.coeffs),function(i){
  goid=cor.coeffs$GOid[i]
  goterm=cor.coeffs[cor.coeffs$GOid==goid,]$GOterm
  if(cor.coeffs[cor.coeffs$GOid==goid,]$tc.size>10){
    (cor0=round(cor.coeffs[cor.coeffs$GOid==goid,]$Spearman.rho,3))
    #(pval0=round(cor.coeffs[cor.coeffs$GOid==goid,]$Pval,6))
    (pval0=round(cor.coeffs[cor.coeffs$GOid==goid,]$FDR,6))
  }else{
    cor0=NA
    #(pval0=round(cor.coeffs[cor.coeffs$GOid==goid,]$Pval,6))
    pval0=NA
  }
  
  ggplot(df5,aes(x=lifespan,y=df5[,paste(goid)]))+
    geom_point(size=3,shape=16)+
    scale_x_log10()+ylab('GO term activity score')+
    ggtitle(paste0('Age ',pick.age,'\n',goid,', ',goterm,'\nKendall\'s tau = ',cor0,', P.adj value=',pval0))+    
    xlab('Cell lifespan (day)')+
    scale_shape_discrete(name='Cell lifespan estimate')
})

plots2=lapply(plots,function(i) 
  i+theme_classic(base_size = 17)+theme(legend.position = 'none',
                                        plot.margin = unit(c(1,1,1,1), "cm"),
                                        #plot.title = element_text(size = 12, face = "bold")))
                                        plot.title = element_text(size = 15)))

pdf(paste0('GO_DNArepair_cell.lifespan_',pick.age,'_rmDups_onePage.pdf'),useDingbats = T,
    width = 6,height = 5,pointsize=12)
for(i in plots2){print(i)}
dev.off()
