
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

######################################################################
## read in mouse turnover rate data
cell.lifespan=data.table::fread('../../src/Dataset_S1.txt')
cell.lifespan$Ron_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`

tc.orders=cell.lifespan[order(cell.lifespan$lifespan),]$`cell type annotation in TMS`
tc.orders #39 tc

##########################################################################
## plot
#need to change!!!
pick.age='3m' 
#pick.age='24m'
#(files=Sys.glob(paste0('mouse_male_binary_TDI_ts_*',pick.age,'_filter.dnds1.rds')))
(files=Sys.glob(paste0('mouse_male_binary_TDI_ts_*',pick.age,'.rds')))

plots1=list()
plots2=list()
res.out=list()
for(file in files){ #max 3.5
  
  cut.off.value=as.numeric(strsplit(stringr::str_extract(file,'ts_.+_\\d+m'),'_')[[1]][[2]])
  
  mouse_tcs_TDI=readRDS(file)
  summary(mouse_tcs_TDI$n_expr_gene)
  mouse_tcs_TDI$tissue_cell.type=droplevels(mouse_tcs_TDI$tissue_cell.type)
  #df=as.data.frame(Reduce(`rbind`,mouse_tcs_TDI))
  #df=mouse_tcs_TDI;
  df=mouse_tcs_TDI[mouse_tcs_TDI$n_expr_gene>=100,]
  table(df$tissue_cell.type)
  #head(df)
  #colnames(df)
  
  # order tc by median TDI at 3m
  one.age=subset(df,age==pick.age);
  x=one.age %>% group_by(tissue_cell.type) %>% summarise(median.TDI=median(TDI),mean.TDI=mean(TDI))
  x=x[order(x$median.TDI),]
  x
  
  cell.meta=df
  cell.meta$tissue_cell.type=factor(cell.meta$tissue_cell.type,levels=as.character(x$tissue_cell.type))

  
  p2=ggplot(subset(cell.meta,age==pick.age),aes(x=tissue_cell.type,y=TDI))+
    #facet_wrap(.~age,ncol=1)+
    geom_jitter(size=0.2,col='grey60')+theme_classic()+
    ylab('Median dN/dS of expressed genes per cell')+
    theme(axis.text.x = element_text(size=8,hjust=1,vjust=0.5,angle=90),
          legend.position = 'none')+xlab('')+
    ggtitle(paste0('Cell type-specific expressed genes, gene specificity score>',cut.off.value))+
    #stat_summary(fun=median, geom="point", shape=16, size=1.2, color="red", fill="red") 
    stat_summary(fun=mean, geom="point", shape=16, size=1.2, color="red", fill="red") 
  #p2  
  
  ## combine lifespan and TDI
  sum(unique(df$tissue_cell.type) %in% cell.lifespan$`cell type annotation in TMS`) #39
  
  #ages=c('3m','18m','24m')
  ages=c('3m','24m')
  out=lapply(ages,function(age){
    one.age=df[df$age==age,]
    one.age=one.age[!is.na(one.age$TDI),]
    #one.age %>% group_by(tissue_cell.type) %>% summarise(median.TDI=median(TDI))
    #https://www.tidyverse.org/blog/2020/03/dplyr-1-0-0-summarise/
    one.age %>% group_by(tissue_cell.type) %>% 
      dplyr::summarise(ncell=n(),ngene=median(n_expr_gene), mean.TDI=mean(TDI),
              #x = quantile(TDI, c(0.25, 0.5, 0.75)), q = c(0.25, 0.5, 0.75)) %>%
              x = quantile(TDI, c( 0.5 )), q = c( 0.5)) %>%
      spread(q,x)  %>% 
      mutate(age=age)
    
  })
  df.out=as.data.frame(Reduce(`rbind`,out))
  head(df.out)
  df2=merge(df.out,cell.lifespan,by.x='tissue_cell.type',by.y='cell type annotation in TMS')
  
  df2$turnover=1-exp(-1/df2$lifespan);
  summary(df2$turnover)
  
  head(df2)
  df2$duplicate='tissue-specific estimate';
  tmp=df2[df2$age==pick.age,]
  df2[df2$Ron_tc %in% tmp$Ron_tc[duplicated(tmp$Ron_tc)],]$duplicate='non tissue-specific estimate'
  dim(df2) 
  table(df2$duplicate)
  
  ## use average for non-tissue specific ones
  #df3=df2 %>% group_by(Ron_tc,age) %>% mutate(average=mean(`0.5`))
  df3=df2 %>% group_by(Ron_tc,age) %>% mutate(average=mean(mean.TDI))
  
  ages=pick.age
  plots<-lapply(ages,function(age){
    tmp=df3[df3$age==age,]
    tmp=tmp[!duplicated(tmp$Ron_tc),]
    dim(tmp) #21
    cat(age,' cor(turnover,n_gene_exp) is ',cor(tmp$lifespan,tmp$ngene,method = 'spearman'),'\n')
    cat(age,' cor(TDI,n_gene_exp) is ',cor(tmp$average,tmp$ngene,method = 'spearman'),'\n')
    cat(age,' cor(TDI,turnover) is ',cor(tmp$average,tmp$lifespan,method = 'spearman'),'\n')
    
    test.out=cor.test(tmp$average,log10(tmp$lifespan),method='spearman',use='pairwise.complete.obs')
    (cor0=test.out$estimate)
    (pval0=test.out$p.value)
    
    test.out=cor.test(tmp$average,log10(tmp$lifespan),method='pearson',use='pairwise.complete.obs')
    (cor1=test.out$estimate)
    (pval1=test.out$p.value)
    
    test.out=cor.test(tmp$average,log10(tmp$lifespan),method='kendall',use='pairwise.complete.obs')
    (cor2=test.out$estimate)
    (pval2=test.out$p.value)
    
    ggplot(tmp,aes(x=lifespan,y=average))+
      #geom_point(aes(col=Ron_tc,shape=duplicate),size=3)+
      geom_point(size=4,shape=16)+
      scale_x_log10()+#ylab(paste0('Median dN/dS per cell type at ',age))+
      ylab(paste0('Mean dN/dS per cell type at ',age))+
      theme_classic(base_size = 15)+
      ggtitle(paste0('Only cell type-specific expressed genes\ngene specificity score > ',cut.off.value,
                     #'\nSpearman\'s rho = ',round(cor0,3),
                     #', P value =',round(pval0,6),
                     '\nPearson\'s r = ',round(cor1,3),
                     ', P value =',round(pval1,6)))+
                     #'\nKendall\'s tau =',round(cor2,3),
                     #', P value =',round(pval2,6)))+
      xlab('Cell lifespan (day)')+
      geom_smooth(method=lm , color="black", fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),
                  se=TRUE,level=0.95) + 
      scale_shape_discrete(name='Cell lifespan estimate')
    #scale_color_viridis(name='Cell type annotation from Sender and Milo (2021)',option='turbo',discrete=T)
    
  })      
  
  #pdf(paste0('TDI_lifespan_binary_ts',cut.off.value,'.pdf'),useDingbats = T,height = 14,width = 12)  
  #grid.arrange(grobs=plots,ncol=2)
  #grid.arrange(p2+ggtitle(paste0('Cell type-specific expressed genes, gene specificity score > ',cut.off.value)),
  #             plots[[1]],ncol=1)
  #dev.off()
  
  plots1[[as.character(cut.off.value)]]=p2+
    ggtitle(paste0('Cell type-specific expressed genes, gene specificity score > ',cut.off.value));
  plots2[[as.character(cut.off.value)]]=plots[[1]];
  
  df3$cut.off.value=cut.off.value
  res.out[[as.character(cut.off.value)]]=df3
  
}

df.res.out=as.data.frame(Reduce(`rbind`,res.out))

#saveRDS(df.res.out,paste0('TDI_binary_ts_all.cutoff.values_',pick.age,'_filter.dnds1.rds'))
saveRDS(df.res.out,paste0('TDI_binary_ts_all.cutoff.values_',pick.age,'.rds'))


df.res.out$cut.off.value #numerical values
df.res.out$cut.off.value=factor(df.res.out$cut.off.value,levels=sort(unique(df.res.out$cut.off.value)))
min(df.res.out$ngene);

#pdf('ts.cutoff_ngene.pdf',useDingbats = T)
#ggplot(df.res.out,aes(x=cut.off.value,y=ngene))+geom_jitter(color='grey',width=0.1)+
#  theme_classic()+scale_y_continuous(limits=c(100,max(df.res.out$ngene)+10))+
#  xlab('Gene specificity score cutoff value')+ylab('Number of expressed genes')
#dev.off()

plots1[[1]] #boxplot
plots2[[1]] #scattor plot
length(plots2) 


plots2.rm.legends=plots2[order(as.numeric(names(plots2)))]

#pdf(paste0('TDI_lifespan_binary_ts_all.cutoff.values_',pick.age,'_filter.dnds1.pdf'),useDingbats = T,height = 14,width = 11)  
pdf(paste0('TDI_lifespan_binary_ts_all.cutoff.values_',pick.age,'.pdf'),useDingbats = T,height = 14,width = 11)  
#grid.arrange(plots2.rm.legends[[1]],plots2.rm.legends[[2]],plots2.rm.legends[[3]],legend,ncol=2)
grid.arrange(grobs=plots2.rm.legends,ncol=2)
dev.off()

########## double axis plot 
#https://github.com/mingwhy/fly.brain.core_coexpr.net/blob/main/comparison/01_plot_rank.aggre_all.four_fdr.R
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
ggplot(df.out,aes(x=Gene.specificity.score.cutoff,y=abs.r))+geom_point(size=5)+theme_classic(base_size = 15)+
  xlab('Gene specificity score cutoff value') + ylab('absolute Pearson\'s correlation coefficient')

df.res.out$cut.off.value=as.numeric(as.character(df.res.out$cut.off.value))

min(df.res.out$ngene)/min(df.out$abs.r)
scaleFactor=4000
p=ggplot()+
  geom_point(aes(x=df.out$Gene.specificity.score.cutoff,y=df.out$abs.r),size=3,col='black')+
  geom_jitter(aes(x=df.res.out$cut.off.value,y=df.res.out$ngene/scaleFactor),color='grey',width=0.1)+
  #theme_classic()+scale_y_continuous(limits=c(100,max(df.res.out$ngene)+10))+
  geom_line(aes(x=df.out$Gene.specificity.score.cutoff,y=df.out$abs.r),col='black',lwd=0.6)+
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

names(plots2.rm.legends)
#pdf(paste0('supp_ts_double_axis_filter_',pick.age,'_dnds1.pdf'),useDingbats = T,width = 6,height = 5)
pdf(paste0('supp_ts_double_axis_filter_',pick.age,'.pdf'),useDingbats = T,width = 6,height = 5)
#print(plots2.rm.legends[['2']]+ggtitle(''))
print(p)
dev.off()

