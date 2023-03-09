
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

##########################################################################
## plot
(files=Sys.glob('mouse_male_binary_TDI_ts_*.rds'))

plots.barplot=list()
plots1=list()
plots2=list()
res.out=list()

for(file in files){ #max 3.5
  (cut.off.value=as.numeric(gsub('ts_|.rds','',stringr::str_extract(file,'ts_.+rds'))))
  
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
  one.age=subset(df,age=='3m');
  x=one.age %>% group_by(tissue_cell.type) %>% summarise(median.TDI=median(TDI),mean.TDI=mean(TDI))
  x=x[order(x$median.TDI),]
  x
  
  cell.meta=df
  cell.meta$tissue_cell.type=factor(cell.meta$tissue_cell.type,levels=as.character(x$tissue_cell.type))

  
  
    p2=ggplot(subset(cell.meta,age=='3m'),aes(x=tissue_cell.type,y=TDI))+
      #facet_wrap(.~age,ncol=1)+
      geom_jitter(size=0.2,col='grey60')+theme_classic()+
      ylab('Median dN/dS of expressed genes per cell at 3m')+
      ggtitle(paste0('Cell type-specific expressed genes, gene specificity score>',cut.off.value))+
      #stat_summary(fun=median, geom="point", shape=16, size=1.2, color="red", fill="red") 
      stat_summary(fun=mean, geom="point", shape=16, size=1.2, color="red", fill="red")+
      ylab('')+xlab('')+
      theme(#axis.text.x = element_text(size=8,hjust=1,vjust=0.5,angle=90),
            axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=10),
            plot.margin = unit(c(1,1,1,2.5), "cm"),
            legend.position = 'none')
    p2
   
    plots.barplot[[as.character(cut.off.value)]]<-p2
  
  
  ## combine lifespan and TDI
  sum(unique(df$tissue_cell.type) %in% cell.lifespan$`tissue: cell.type in mouse`) #39
  
  #ages=c('3m','18m','24m')
  ages=c('3m')
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
  df2=merge(df.out,cell.lifespan,by.x='tissue_cell.type',by.y='tissue: cell.type in mouse')
  
  df2$turnover=1-exp(-1/df2$lifespan);
  summary(df2$turnover)
  
  head(df2)
  df2$duplicate='tissue-specific estimate';
  tmp=df2[df2$age=='3m',]
  df2[df2$human_tc %in% tmp$human_tc[duplicated(tmp$human_tc)],]$duplicate='non tissue-specific estimate'
  dim(df2) 
  table(df2$duplicate)
  
  ## use average for non-tissue specific ones
  #df3=df2 %>% group_by(human_tc,age) %>% mutate(average=mean(`0.5`))
  df3=df2 %>% group_by(human_tc,age) %>% mutate(average=mean(mean.TDI))
  
  ages='3m'
  plots<-lapply(ages,function(age){
    tmp=df3[df3$age==age,]
    tmp=tmp[!duplicated(tmp$human_tc),]
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
      #geom_point(aes(col=human_tc,shape=duplicate),size=3)+
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
    
  pdf(paste0('TDI_lifespan_binary_ts',cut.off.value,'.pdf'),useDingbats = T,height = 14,width = 12)
  grid.arrange(p2+ggtitle(paste0('Cell type-specific expressed genes, gene specificity score > ',cut.off.value)),
               plots[[1]],ncol=1)
  dev.off()
  
  plots1[[as.character(cut.off.value)]]=p2+
    ggtitle(paste0('Cell type-specific expressed genes, gene specificity score > ',cut.off.value));
  plots2[[as.character(cut.off.value)]]=plots[[1]];
  
  df3$cut.off.value=cut.off.value
  res.out[[as.character(cut.off.value)]]=df3
  
}

length(plots.barplot)
plots.barplot<-plots.barplot[order(as.numeric(names(plots.barplot)))]

jpeg(file="TDI_by.cell.type_binary_boxplot.jpeg",width=2400, height=3000,quality = 20000,res= 120)
grid.arrange(grobs=plots.barplot,ncol=2)
dev.off()

df.res.out=as.data.frame(Reduce(`rbind`,res.out))
saveRDS(df.res.out,'TDI_binary_ts_all.cutoff.values.rds')

