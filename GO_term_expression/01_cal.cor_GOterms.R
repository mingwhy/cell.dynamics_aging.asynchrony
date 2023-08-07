
####################################
#library(zellkonverter)
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
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt') 
gene.meta=id.mapping
####################################
## read in mouse turnover rate data
cell.lifespan=data.table::fread('../src/Dataset_S1.txt')
cell.lifespan$Ron_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`
tc.orders=cell.lifespan[order(cell.lifespan$lifespan),]$`cell type annotation in TMS`
tc.orders #39 tc

########################################################################################
(files=Sys.glob('select.tc.h5ad_result_goterms/mouse_male_GO_*_libNorm2_*.rds'))

cor.coeffs.all=list()
df5.all=list();

all.plots=list();
for(output_file in files){
  pick.age=gsub('GO_|_','',stringr::str_extract(output_file,'GO_.+?_'))
  
  
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
  dim(df4) #21 x 5
  
  # at least 10 cell types terms express this GO term's gene members
  tc.size=unlist(lapply(2:ncol(df4),function(i) sum(df4[,i]!=0) ))
  
  
  intersect(colnames(df4),colnames(df2)) #Ron_tc
  df5=merge(df4,df2[!duplicated(df2$Ron_tc),],all.x=TRUE)
  dim(df5) #21
  
  Matrix=df4[,-1];
  dim(Matrix) # ncell.type x 1061 GO
  
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
  cor.coeffs$age=pick.age
  cor.coeffs$tc.size=tc.size
  cor.coeffs$filename=output_file;
  
  dim(cor.coeffs) #6 GO term x _
  x=lapply(cor.coeffs$GOid,function(i) GOTERM[[i]]@Term)
  cor.coeffs$GOterm=unlist(x)
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
                                          plot.title = element_text(size = 14)))  
  
  pdf(gsub('.rds','.pdf',basename(output_file)),useDingbats = T,
      width = 7,height = 5,pointsize=12)
  for(i in plots2){print(i)}
  dev.off()
  
  df5$age=pick.age
  cor.coeffs.all[[basename(output_file)]]=cor.coeffs
  df5.all[[basename(output_file)]]=df5
}

cor.coeffs.df=Reduce(`rbind`,cor.coeffs.all)
length(unique(cor.coeffs.df$GOid)) #17
dim(cor.coeffs.df) #34
sum(cor.coeffs.df$FDR<0.05) #10
#cor.coeffs.df[cor.coeffs.df$FDR<0.05,]
saveRDS(cor.coeffs.df,'cor.coeffs.df.rds')

length(df5.all) #6
saveRDS(df5.all,'GOscore_tc.rds')

