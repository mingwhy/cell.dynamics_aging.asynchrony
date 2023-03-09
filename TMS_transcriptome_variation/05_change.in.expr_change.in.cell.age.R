
library(tidyverse);
library(ggplot2);library(gridExtra)
library(viridis)
library(ggplot2) 
library(grid)
library(gridExtra) 

source('src_cellular_lifespan.R')
dim(cell.lifespan)

tissues=sapply(strsplit(cell.lifespan$`cell type annotation in TMS`,'\\:'),'[[',1)
tcs=sapply(strsplit(cell.lifespan$`cell type annotation in Sender and Milo 2021`,'\\:'),'[[',2)

x=names(which(table(cell.lifespan$`cell type annotation in Sender and Milo 2021`)>1))
cell.lifespan$duplicate='tissue-specific estimate';
cell.lifespan[cell.lifespan$`cell type annotation in Sender and Milo 2021` %in% x,]$duplicate='non tissue-specific estimate';
table(cell.lifespan$duplicate)

tcs[grep('endothelial',tcs,ignore.case = T)]='endothelial cells'
cell.lifespan$cell.identity=tcs
cell.lifespan$tissues=tissues;

cell.lifespan$human_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`
######################################################################
## read in mean.diff and sd.diff in theoretical model
para_df=readRDS('../cell_demography_model/plug-in_Rcal_mean_var_N0F.rds')

para_df2<-para_df %>% group_by(human_tc) %>% 
  summarise(mean.ratio= mean[month==24] / mean[month==3],
            mean.diff= mean[month==24] - mean[month==3],
            sd.ratio= var[month==24]^0.5 / var[month==3]^0.5,
            sd.diff= var[month==24]^0.5 - var[month==3]^0.5);

intersect(colnames(cell.lifespan),colnames(para_df2)) #`cell type annotation in Sender and Milo 2021`

df.model=merge(para_df2,cell.lifespan)

###################################################
# read in aging magnitude 
df=data.table::fread('age.change_all.metrics.txt')
dim(df) #39
df$change.in.expr=df$gene_CV2

######################################################################
## integrate expr.change and age.change 
sum(df$cell.type %in% df.model$`cell type annotation in TMS`) #39

df.model$cell.type=df.model$`cell type annotation in TMS`;
intersect(colnames(df.model),colnames(df)) #cell.type

df2=merge(df.model[,c('cell.type','cell.identity',"duplicate" ,'cell type annotation in Sender and Milo 2021',
                      'tissues','lifespan',
               'mean.diff','mean.ratio','sd.diff','sd.ratio')],
      df[,c('cell.type','change.in.expr')]) 
# change lifespan into turnover
df2$turnover=1-exp(-1/df2$lifespan)

colnames(df2)
df2=df2[order(df2$lifespan),]

##################################################
## calculate correlation
table(df2$duplicate)
df3=df2 %>% group_by(`cell type annotation in Sender and Milo 2021`) %>% mutate('change.in.expr.var'=mean(change.in.expr))
df4=df3[!duplicated(df3$`cell type annotation in Sender and Milo 2021`),]
dim(df4) #21

test.out=cor.test(df4$mean.diff,df4$change.in.expr.var,method='spearman')
(cor=test.out$estimate)
(pval=test.out$p.value)
test.out=cor.test(df4$mean.diff,df4$change.in.expr.var,method='pearson')
(cor.pearson=test.out$estimate) 
(pval.pearson=test.out$p.value)


p2=ggplot(df4,aes(x=mean.diff,y=change.in.expr.var))+
  geom_jitter(size=5,shape=16,width = 6)+
  theme_classic(base_size = 25)+
    ggtitle(paste0('Pearson\'s r = ',round(cor.pearson,3),
                  ', P value=',round(pval.pearson,8)))+
  geom_smooth(method=lm , color="black", fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),
                                                                   se=TRUE,level=0.95) + 
  scale_shape_manual(name='Cell turnover estimate',values=c(1,2))
 

pdf('change.in.transcriptome_change.in.age.pdf',useDingbats = T,width = 9,height = 8)
p2+ylab(expression(paste('Change in transcriptome variability (',beta[age],')')))+#xlab('mean(24m)-mean(3m)')+
  xlab('Change in cell age (day)')+
  #theme(legend.position = 'none',plot.title = element_text(size = 15, face = "bold"))
  #theme(plot.title = element_text(size = 15, face = "bold"))
  theme(plot.title = element_text(size = 0, face = "bold"))
p2+ylab(expression(paste('Change in transcriptome variability (',beta[age],')')))+#xlab('mean(24m)-mean(3m)')+
  xlab('Change in cell age (day)')+
  theme(plot.title = element_text(size = 10, face = "bold"))
dev.off()




