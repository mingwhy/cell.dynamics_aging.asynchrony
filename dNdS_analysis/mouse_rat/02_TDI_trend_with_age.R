
library(gridExtra) 
#library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);library(grid)
library(tidyverse);
library(zellkonverter)
library(SummarizedExperiment)
library(ggpointdensity)
library(viridis);library(RColorBrewer);library(grDevices)
library(ggpubr)

my.cols=c( "#56B4E9","#E69F00" ); # blue (young),yellow (old)

######################################################################
## read in mouse turnover rate data
cell.lifespan=data.table::fread('../../src/Dataset_S1.txt')
cell.lifespan$Ron_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`

tc.orders=cell.lifespan[order(cell.lifespan$lifespan),]$`cell type annotation in TMS`
tc.orders #39 tc

######################################################################
## read in TDI
#cut.off.value=2;

df.res.out.young=readRDS('./filtered.dnds1/TDI_binary_ts_all.cutoff.values_3m.rds')
#df.res.out.young=readRDS('./all.genes/TDI_binary_ts_all.cutoff.values_3m.rds')
df.res.out.young$cut.off.value
df.res.out.young$age
summary(df.res.out.young[df.res.out.young$cut.off.value==2,]$ngene)
#116-971
summary(df.res.out.young[df.res.out.young$cut.off.value==2,]$ncell)
#56~2251

df.res.out.old=readRDS('./filtered.dnds1/TDI_binary_ts_all.cutoff.values_24m.rds')
#df.res.out.old=readRDS('./all.genes/TDI_binary_ts_all.cutoff.values_24m.rds')
df.res.out.old$cut.off.value
df.res.out.old$age

################################################################
pdf('TDI_lifespan_binary_regression_2ages_two_filter.dnds1.pdf',useDingbats = T,height = 6,width = 10)
#pdf('TDI_lifespan_binary_regression_2ages_two.pdf',useDingbats = T,height = 6,width = 10)
for(cut.off.value in sort(unique(df.res.out.young$cut.off.value))){
 
  df.res.out.3m=df.res.out.young[df.res.out.young$cut.off.value==cut.off.value,]
  
  df.res.out.24m=df.res.out.old[df.res.out.old$cut.off.value==cut.off.value,]
  
  my.list=list(df.res.out.3m[!duplicated(df.res.out.3m$Ron_tc),],
                   df.res.out.24m[!duplicated(df.res.out.24m$Ron_tc),])
  
  plots=lapply(1:2,function(i){
    tmp=my.list[[i]]
    #tmp=df.res.out.3m
    age=tmp$age[1]
    test.out=cor.test(tmp$average,log10(tmp$lifespan),method='pearson',use='pairwise.complete.obs')
    (cor1=test.out$estimate)
    (pval1=test.out$p.value)
    ggplot(tmp,aes(x=lifespan,y=average))+
      #geom_point(aes(col=Ron_tc,shape=duplicate),size=3)+
      geom_point(size=4,shape=16)+
      scale_x_log10()+#ylab(paste0('Median dN/dS per cell type at ',age))+
      ylab(paste0('Mean dN/dS per cell type at ',age))+
      theme_classic(base_size = 15)+
      ggtitle(paste0('Age ',age,'\nOnly cell type-specific expressed genes\ngene specificity score > ',cut.off.value,
                     #'\nSpearman\'s rho = ',round(cor0,3),
                     #', P value =',round(pval0,6),
                     '\nPearson\'s r = ',round(cor1,3),
                     ', P value =',round(pval1,6)))+
      #'\nKendall\'s tau =',round(cor2,3),
      #', P value =',round(pval2,6)))+
      xlab('Cell lifespan (day)')+
      geom_smooth(method=lm , color=my.cols[i], fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),
                  se=TRUE,level=0.95) + 
      scale_shape_discrete(name='Cell lifespan estimate')+
      theme(axis.title = element_text(size=17),
            plot.title = element_text(size=10))
  })
  #plots[[1]]
  
  grid.arrange(grobs=plots,ncol=2)
}
dev.off()

######################################################################
## select cutoff=2
cut.off.value=2
df.res.out.3m=df.res.out.young[df.res.out.young$cut.off.value==cut.off.value,]
df.res.out.24m=df.res.out.old[df.res.out.old$cut.off.value==cut.off.value,]
my.list=list(df.res.out.3m[!duplicated(df.res.out.3m$Ron_tc),],
              df.res.out.24m[!duplicated(df.res.out.24m$Ron_tc),])
plots.cutoff2=lapply(1:2,function(i){
    tmp=my.list[[i]]
    age=tmp$age[1]
    test.out=cor.test(tmp$average,log10(tmp$lifespan),method='pearson',use='pairwise.complete.obs')
    (cor1=test.out$estimate)
    (pval1=test.out$p.value)
    ggplot(tmp,aes(x=lifespan,y=average))+
      #geom_point(aes(col=Ron_tc,shape=duplicate),size=3)+
      geom_point(size=4,shape=16)+
      scale_x_log10()+#ylab(paste0('Median dN/dS per cell type at ',age))+
      ylab(paste0('Mean dN/dS per cell type at ',age))+
      theme_classic(base_size = 15)+
      ggtitle(paste0('Age ',age,'\nOnly cell type-specific expressed genes\ngene specificity score > ',cut.off.value,
                     '\nPearson\'s r = ',round(cor1,3),
                     ', P value =',round(pval1,6)))+
      xlab('Cell lifespan (day)')+
      geom_smooth(method=lm , color=my.cols[i], fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),
                  se=TRUE,level=0.95) + 
      scale_shape_discrete(name='Cell lifespan estimate')+
      theme(axis.title = element_text(size=17),
            plot.title = element_text(size=10))
})
grid.arrange(grobs=plots.cutoff2,ncol=2)

####################################################################################################################
## double axis plot for two ages
#https://github.com/mingwhy/fly.brain.core_coexpr.net/blob/main/comparison/01_plot_rank.aggre_all.four_fdr.R

df.res.out.young=readRDS('./filtered.dnds1/TDI_binary_ts_all.cutoff.values_3m.rds')
df.res.out.old=readRDS('./filtered.dnds1/TDI_binary_ts_all.cutoff.values_24m.rds')
#df.res.out.young=readRDS('./all.genes/TDI_binary_ts_all.cutoff.values_3m.rds')
#df.res.out.old=readRDS('./all.genes/TDI_binary_ts_all.cutoff.values_24m.rds')


## handle young
df.res.out=df.res.out.young
cuts=sort(as.numeric(as.character(unique(df.res.out$cut.off.value))))
out=sapply(cuts,function(i){
  tmp=df.res.out[df.res.out$cut.off.value==i,] %>% group_by(Ron_tc,age) %>% mutate(average=mean(mean.TDI))
  tmp=tmp[!duplicated(tmp$Ron_tc),]
  dim(tmp) #21
  
  test.out=cor.test(tmp$average,log10(tmp$lifespan),method='pearson',use='pairwise.complete.obs')
  (cor1=test.out$estimate)
  (pval1=test.out$p.value)
  return(c(i,cor1,pval1))
})
df.out=as.data.frame(t(out))
colnames(df.out)=c('Gene.specificity.score.cutoff','pearson.r','pval')
df.out$abs.r=abs(df.out$pearson.r)
#ggplot(df.out,aes(x=Gene.specificity.score.cutoff,y=abs.r))+geom_point(size=5)+theme_classic(base_size = 15)+
#  xlab('Gene specificity score cutoff value') + ylab('absolute Pearson\'s correlation coefficient')

## handle old
df.res.out=df.res.out.old
cuts=sort(as.numeric(as.character(unique(df.res.out$cut.off.value))))
out=sapply(cuts,function(i){
  tmp=df.res.out[df.res.out$cut.off.value==i,] %>% group_by(Ron_tc,age) %>% mutate(average=mean(mean.TDI))
  tmp=tmp[!duplicated(tmp$Ron_tc),]
  dim(tmp) #21
  
  test.out=cor.test(tmp$average,log10(tmp$lifespan),method='pearson',use='pairwise.complete.obs')
  (cor1=test.out$estimate)
  (pval1=test.out$p.value)
  return(c(i,cor1,pval1))
})
df.out.old=as.data.frame(t(out))
colnames(df.out.old)=c('Gene.specificity.score.cutoff','pearson.r','pval')
df.out.old$abs.r=abs(df.out.old$pearson.r)


df.res.out$cut.off.value=as.numeric(as.character(df.res.out$cut.off.value))

min(df.res.out$ngene)/min(df.out$abs.r)
scaleFactor=4000
p=ggplot()+
  geom_point(aes(x=df.out$Gene.specificity.score.cutoff,y=df.out$abs.r),size=3,col=my.cols[1])+
  geom_line(aes(x=df.out$Gene.specificity.score.cutoff,y=df.out$abs.r),col=my.cols[1],lwd=0.6)+
  
  geom_point(aes(x=df.out.old$Gene.specificity.score.cutoff,y=df.out.old$abs.r),size=3,col=my.cols[2])+
  geom_line(aes(x=df.out.old$Gene.specificity.score.cutoff,y=df.out.old$abs.r),col=my.cols[2],lwd=0.6)+
  
  geom_jitter(aes(x=df.res.out$cut.off.value,y=df.res.out$ngene/scaleFactor),color='grey',width=0.1)+
  #theme_classic()+scale_y_continuous(limits=c(100,max(df.res.out$ngene)+10))+
  
  #stat_summary(aes(x=tmp1$cutoff,y=tmp1$logP),fun.y=median, geom="point", size=1, color="black")+
  xlab('Gene specificity score cutoff value')+
  scale_y_continuous(
    name = "Absolute Pearson\'s correlation coefficient", # Features of the first axis
    sec.axis = sec_axis(trans=~.*scaleFactor,
                        name="Number of expressed genes"))+  # Add a second axis and specify its features
  theme_classic(base_size = 15)+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        legend.position = 'none')
p

length(plots.cutoff2)
plots.cutoff2[[3]]<-p
library(egg)

pdf('supp_ts_double_axis_filtered.dnds1.pdf',useDingbats = T,width = 14,height = 5.5)
ggarrange(plots.cutoff2[[1]],plots.cutoff2[[2]],plots.cutoff2[[3]],ncol=3,
          labels = c("A", "B","C"))
#grid.arrange(grobs=plots.cutoff2,ncol=3)
dev.off()



