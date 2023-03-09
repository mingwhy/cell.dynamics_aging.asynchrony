
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)

source('src_cellular_lifespan.R')
dim(cell.lifespan)

cell.lifespan$`tissue: cell.type in mouse`=cell.lifespan$`cell type annotation in TMS`
cell.lifespan$human_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`

tissues=sapply(strsplit(cell.lifespan$`tissue: cell.type in mouse`,'\\:'),'[[',1)
tcs=sapply(strsplit(cell.lifespan$human_tc,'\\:'),'[[',2)

x=names(which(table(cell.lifespan$human_tc)>1))
cell.lifespan$duplicate='tissue-specific estimate';
cell.lifespan[cell.lifespan$human_tc %in% x,]$duplicate='non tissue-specific estimate';
table(cell.lifespan$duplicate)

tcs[grep('endothelial',tcs,ignore.case = T)]='endothelial cells'
cell.lifespan$cell.identity=tcs
cell.lifespan$tissues=tissues;

######################################################################
## read in mean.diff and sd.diff in theoretical model
para_df=readRDS('../cell_demography_model/plug-in_Rcal_mean_var_N0F.rds')

para_df2<-para_df %>% group_by(human_tc) %>% 
  summarise(mean.ratio= mean[month==24] / mean[month==3],
            mean.diff= mean[month==24] - mean[month==3],
            sd.ratio= var[month==24]^0.5 / var[month==3]^0.5,
            sd.diff= var[month==24]^0.5 - var[month==3]^0.5);

intersect(colnames(cell.lifespan),colnames(para_df2)) #human_tc

df.model=merge(para_df2,cell.lifespan)

###################################################
# read in aging magnitude 
df=data.table::fread('age.change_all.metrics.txt')

######################################################################
## integrate expr.change and age.change 
sum(df$cell.type %in% df.model$`tissue: cell.type in mouse`) #39

df.model$cell.type=df.model$`tissue: cell.type in mouse`;
intersect(colnames(df.model),colnames(df)) #cell.type

df2=merge(df.model[,c('cell.type','cell.identity',"duplicate" ,'human_tc',
                      'tissues','lifespan',
                      'mean.diff','mean.ratio','sd.diff','sd.ratio')],
          df)
# change lifespan into turnover
df2$turnover=1-exp(-1/df2$lifespan)

colnames(df2)
df2=df2[order(df2$lifespan),]

################################################################v
## calculate correlation and plot 1-spearman
table(df2$duplicate)
df3=df2 %>% group_by(human_tc) %>% mutate('change.in.expr.var'=mean(`1-spearman`))
df4=df3[!duplicated(df3$human_tc),]
dim(df4) #21

test.out=cor.test(df4$mean.diff,df4$change.in.expr.var,method='spearman') 
(spearman.cor.value=test.out$estimate)
(pval=test.out$p.value)

test.out=cor.test(df4$mean.diff,df4$change.in.expr.var,method='pearson') 
(pearson.cor.value=test.out$estimate)
(pearson.pval=test.out$p.value)


p.rho=ggplot(df4,aes(x=mean.diff,y=change.in.expr.var))+
  geom_jitter(size=5,shape=16,width = 6)+
  theme_classic(base_size = 15)+
  geom_smooth(method=lm , color="black", fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),
              se=TRUE,level=0.95) + 
  ylab('Aging magnitude')+xlab('mean(24m)-mean(3m)')
  ggtitle(paste0('Pearson\'s r = ',round(pearson.cor.value,3),
                 ', P value=',round(pearson.pval,6)))
p.rho<-p.rho+ylab(expression(paste('Change in transcriptome variability (1-Spearman\'s ',rho,')')))+#xlab('mean(24m)-mean(3m)')+
  xlab('Change in cell age (day)')

################################################################v
## calculate correlation and plot 1-gcl
table(df2$duplicate)
df3=df2 %>% group_by(human_tc) %>% mutate('change.in.expr.var'=mean(`1-GCL`))
df4=df3[!duplicated(df3$human_tc),]
dim(df4) #21

test.out=cor.test(df4$mean.diff,df4$change.in.expr.var,method='spearman') 
(spearman.cor.value=test.out$estimate)
(pval=test.out$p.value)

test.out=cor.test(df4$mean.diff,df4$change.in.expr.var,method='pearson') 
(pearson.cor.value=test.out$estimate)
(pearson.pval=test.out$p.value)

p.gcl=ggplot(df4,aes(x=mean.diff,y=change.in.expr.var))+
  geom_jitter(size=5,shape=16,width = 6)+
  theme_classic(base_size = 15)+
  geom_smooth(method=lm , color="black", fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),
              se=TRUE,level=0.95) + 
  scale_shape_manual(name='Cell turnover estimate',values=c(1,2))+
  ylab('Aging magnitude')+xlab('mean(24m)-mean(3m)')+
  ggtitle(paste0(#'cell-to-cell heterogeneity (1-GCL )',
                 #'\nSpearman.cor.coeff=',round(spearman.cor.value,3),
                 #', P value=',round(pval,6),
                 'Pearson\'s r = ',round(pearson.cor.value,3),
                 ', P value=',round(pearson.pval,6)))
p.gcl<-p.gcl+ylab('Change in transcriptome variability (1-GCL)')+#xlab('mean(24m)-mean(3m)')+
  xlab('Change in cell age (day)')


pdf('alternative.metric_log2ratio_change.in.cell.age.pdf',useDingbats = T,width = 8,height = 11)
grid.arrange(p.rho,p.gcl,ncol=1)
dev.off()

