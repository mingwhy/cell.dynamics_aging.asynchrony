
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)
#https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
scaleFUN <- function(x) sprintf("%.3f", x)

age.groups=c('3m','24m')
output.file1=paste0('age.change_all.metrics_per.mouse_',age.groups[1],'_vs_',age.groups[2],'.txt')
output.file2=paste0('s2_',age.groups[1],'_vs_',age.groups[2],'.txt')
out.fig=paste0('all.metric_per.mouse_',age.groups[1],'_vs_',age.groups[2],'.pdf')

cell.lifespan=data.table::fread('../src/Dataset_S1.txt')
cell.lifespan$Ron_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`
tc.orders=cell.lifespan[order(cell.lifespan$lifespan),]$`cell type annotation in TMS`
tc.orders

tissues=sapply(strsplit(cell.lifespan$`cell type annotation in TMS`,'\\:'),'[[',1)
tcs=sapply(strsplit(cell.lifespan$Ron_tc,'\\:'),'[[',2)

x=names(which(table(cell.lifespan$Ron_tc)>1))
cell.lifespan$duplicate='tissue-specific estimate';
cell.lifespan[cell.lifespan$Ron_tc %in% x,]$duplicate='non tissue-specific estimate';
table(cell.lifespan$duplicate)

tcs[grep('endothelial',tcs,ignore.case = T)]='endothelial cells'
cell.lifespan$cell.identity=tcs
cell.lifespan$tissues=tissues;


if(!file.exists(output.file1)){
  
  #############################################################
  ## gene-level metric
  res.gene=readRDS('LMM_cv2_out.rds')
  dim(res.gene) 
  head(res.gene)
  colnames(res.gene)[1]='cell.type'
  
  res.gene=res.gene[res.gene$cell.type %in% tc.orders,]
  res.gene=res.gene[order(res.gene$cell.type),]
  
  head(res.gene)
  summary(p.adjust(res.gene$pvalue,method='BH')) 
  #Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  #0.0004441 0.0040813 0.0211100 0.1438579 0.1678324 0.9962423
  res.gene$FDR=p.adjust(res.gene$pvalue,method='BH')
  sum(res.gene$FDR<0.05 & res.gene$beta>0) #22 tc
  sum(res.gene$FDR<0.05 & res.gene$beta<0) #1 tc
  dim(res.gene[res.gene$FDR>=0.05,]) #16 tc
  
  ############################################################
  ## cell: res_cell.cors.rds
  df=readRDS('res_cell.cors_per.mouse.rds')
  df=df[df$age %in% age.groups,]
  df$age=droplevels(df$age)
  df$dist=1-df$cell.cors;
  
  df.pval=df %>% group_by(cell.type) %>%  
    summarize(pval = wilcox.test(dist ~ age)$p.value)
  df.pval$p.adj=p.adjust(df.pval$pval,method='BH')
  
  tmp0=df %>% group_by(cell.type,age,mouse.id) %>% summarise(age.average.per.mouse=median(dist))
  tmp=tmp0 %>% group_by(cell.type,age) %>% summarise(age.average=mean(age.average.per.mouse))
 
  tmp1 = tmp %>% spread(age,age.average)
  tmp1=tmp1[!is.na(tmp1[[age.groups[[2]]]]),]
  tmp1$log2ratio=log2(tmp1[[age.groups[[2]]]]/tmp1[[age.groups[[1]]]])
  
  res.spearman=tmp1;
  res.spearman=res.spearman[res.spearman$cell.type %in% tc.orders,]
  unique(res.spearman$cell.type) #39 tc
  res.spearman=res.spearman[order(res.spearman$cell.type),]
  
  res.spearman=merge(res.spearman,df.pval)
  sum(res.spearman$log2ratio>0) #29
  sum(res.spearman$log2ratio<0) #10
  sum(res.spearman$log2ratio>0 & res.spearman$p.adj<0.05) #25
  summary(res.spearman[res.spearman$log2ratio>0 & res.spearman$p.adj<0.05,]$p.adj)
  
  ############################################################
  ## cell: GCL
  df=data.table::fread('sce.minCellCounts_gcl_per.mouse.txt')
  df$cell.type=df$tc
  df$dist=1-df$gcl;
  df  %>% group_by(cell.type,age) %>% summarize(n()) #100 reps per tc
  
  df.pval=df %>% group_by(cell.type) %>%  
    summarize(pval = wilcox.test(dist ~ age)$p.value)
  df.pval$p.adj=p.adjust(df.pval$pval,method='BH')
  
  tmp0=df %>% group_by(cell.type,age,mouse.id) %>% summarise(age.average.per.mouse=mean(dist))
  tmp=tmp0 %>% group_by(cell.type,age) %>% summarise(age.average=mean(age.average.per.mouse))
  tmp1 = tmp %>% spread(age,age.average)
  sum(is.na(tmp1[[age.groups[[1]]]]))
  sum(is.na(tmp1[[age.groups[[2]]]]))
  tmp1$log2ratio=log2(tmp1[[age.groups[[2]]]]/tmp1[[age.groups[[1]]]])
  
  res.gcl=tmp1;
  res.gcl=res.gcl[res.gcl$cell.type %in% tc.orders,]
  unique(res.gcl$cell.type) #39tc
  res.gcl=res.gcl[order(res.gcl$cell.type),]
  
  res.gcl=merge(res.gcl,df.pval)
  sum(res.gcl$log2ratio>0) #24
  sum(res.gcl$log2ratio<0) #15
  sum(res.gcl$log2ratio>0 & res.gcl$p.adj<0.05) #17
  summary(res.gcl[res.gcl$log2ratio>1 & res.gcl$p.adj<0.05,]$p.adj)
  #Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
  #0.000e+00 0.000e+00 0.000e+00 4.649e-15 4.649e-15 1.860e-14 
  
  #######################################################################
  ## make sure cell.types are all the same order among different metrics
  dim(res.spearman);sum(res.spearman$cell.type==res.gcl$cell.type)
  sum(res.spearman$cell.type==res.gene$cell.type);
  
  df=data.frame(res.spearman$log2ratio,res.gcl$log2ratio,res.gene$beta)
  colnames(df)=c('1-spearman','1-GCL','gene_CV2')
  
  df$cell.type=res.gene$cell.type

  data.table::fwrite(df,output.file1)
}

###################################################
# read in aging magnitude 
df=data.table::fread(output.file1)

######################################################################
## read in mean.diff and sd.diff in theoretical model
age.groups.numeric=as.numeric(gsub('m','',age.groups))
para_df=readRDS(paste0('../src/change.in.cell.age_',age.groups.numeric[1],
                       '_vs_',age.groups.numeric[2],'.rds'));

para_df2<-para_df %>% group_by(Ron_tc) %>% 
  summarise(mean.ratio= mean[month==age.groups.numeric[2]] / mean[month==age.groups.numeric[1]],
            mean.diff= mean[month==age.groups.numeric[2]] - mean[month==age.groups.numeric[1]],
            sd.ratio= var[month==age.groups.numeric[2]]^0.5 / var[month==age.groups.numeric[1]]^0.5,
            sd.diff= var[month==age.groups.numeric[2]]^0.5 - var[month==age.groups.numeric[1]]^0.5);

intersect(colnames(cell.lifespan),colnames(para_df2)) #Ron_tc

df.model=merge(para_df2,cell.lifespan)
colnames(df.model)
head(df.model)
dim(df.model) #39

######################################################################
## integrate expr.change and age.change 
sum(df$cell.type %in% df.model$`cell type annotation in TMS`) #39

df.model$cell.type=df.model$`cell type annotation in TMS`;
intersect(colnames(df.model),colnames(df)) #cell.type

df2=merge(df.model[,c('cell.type','cell.identity',"duplicate" ,'Ron_tc',
                      'tissues','lifespan',
                      'mean.diff','mean.ratio','sd.diff','sd.ratio')],
          df)
# change lifespan into turnover
df2$turnover=1-exp(-1/df2$lifespan)

colnames(df2)
df2=df2[order(df2$lifespan),]

# save a table (supp.Table2)
colnames(df2)[c(1,4,3,6,7,13,11,12)]
df3=df2[,c(1,4,3,6,7,13,11,12)]

tmp=df3;
x1=para_df[para_df$month==3,c('Ron_tc','mean')]
x2=para_df[para_df$month==24,c('Ron_tc','mean')]
tmp=merge(df3,x1)
colnames(tmp)[match('mean',colnames(tmp))]='mean cell age at 3-month (in days)'

tmp1=merge(tmp,x2)
colnames(tmp1)[match('mean',colnames(tmp1))]='mean cell age at 24-month (in days)'
tmp1=tmp1[order(tmp1$lifespan),]
colnames(tmp1)[c(2,1,3,4,9,10,5:8)]
#tmp1=tmp1[,c(2,1,3:10)]

tmp1=tmp1[,c(2,1,3,4,9,10,5:8)]
data.table::fwrite(tmp1,output.file2,sep='\t')


##################################################
## calculate correlation
table(df2$duplicate)
df3=df2 %>% group_by(Ron_tc) %>% mutate('change.in.expr.var'=mean(gene_CV2))
df4=df3[!duplicated(df3$Ron_tc),]
dim(df3); #38
dim(df4); #20
table(df4$duplicate) 

test.out=cor.test(df4$mean.diff,df4$change.in.expr.var,method='spearman')
(cor2=test.out$estimate)
(pval2=test.out$p.value)
test.out=cor.test(df4$mean.diff,df4$change.in.expr.var,method='pearson')
(cor2.pearson=test.out$estimate) 
(pval2.pearson=test.out$p.value)

p.cv=ggplot(df4,aes(x=mean.diff,y=change.in.expr.var))+
  geom_jitter(aes(fill=duplicate),pch=21,size=5,width=6)+
  theme_classic(base_size = 18)+
  ggtitle(paste0('Pearson\'s r = ',round(cor2.pearson,3),
                 ', P value=',format.pval(pval2.pearson,digits=2)))+
  geom_smooth(method=lm , color="black", fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),
              se=TRUE,level=0.95) + 
  scale_fill_manual(name='Cell lifespan estimate',values=c('NA','black'))+
  scale_y_continuous(limits = c(-0.015,0.030),labels=scaleFUN)
p.cv

################################################################v
## calculate correlation and plot 1-spearman
table(df2$duplicate)
df3=df2 %>% group_by(Ron_tc) %>% mutate('change.in.expr.var'=mean(`1-spearman`))
df4=df3[!duplicated(df3$Ron_tc),]
dim(df4) #22

test.out=cor.test(df4$mean.diff,df4$change.in.expr.var,method='spearman') 
(spearman.cor.value=test.out$estimate)
(pval=test.out$p.value)

test.out=cor.test(df4$mean.diff,df4$change.in.expr.var,method='pearson') 
(pearson.cor.value=test.out$estimate)
(pearson.pval=test.out$p.value)

p.rho=ggplot(df4,aes(x=mean.diff,y=change.in.expr.var))+
  geom_jitter(aes(fill=duplicate),pch=21,size=5,width=6)+
  theme_classic(base_size = 18)+
  geom_smooth(method=lm , color="black", fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),
              se=TRUE,level=0.95) + 
  scale_fill_manual(name='Cell turnover estimate',values=c('white','black'))+
  ylab('Aging magnitude')+xlab('mean(24m)-mean(3m)')+#ylab( 'log2(ratio)')+
  ggtitle(paste0('Pearson\'s r = ',round(pearson.cor.value,3),
                 ', P value=',format.pval(pearson.pval,digits=2) ))+
  scale_y_continuous(limits=c(-0.155,0.175),labels=scaleFUN)
p.rho<-p.rho+ylab(expression(paste('Change in transcriptome variability (1-Spearman\'s ',rho,')')))+#xlab('mean(24m)-mean(3m)')+
  xlab(expression(paste(Delta,' age (day)'))) #xlab('Change in cell age (day)')
p.rho

################################################################
## calculate correlation and plot 1-gcl
table(df2$duplicate)
df3=df2 %>% group_by(Ron_tc) %>% mutate('change.in.expr.var'=mean(`1-GCL`))
df4=df3[!duplicated(df3$Ron_tc),]
dim(df4) #22

test.out=cor.test(df4$mean.diff,df4$change.in.expr.var,method='spearman') 
(spearman.cor.value=test.out$estimate)
(pval=test.out$p.value)

test.out=cor.test(df4$mean.diff,df4$change.in.expr.var,method='pearson') 
(pearson.cor.value=test.out$estimate)
(pearson.pval=test.out$p.value)

p.gcl=ggplot(df4,aes(x=mean.diff,y=change.in.expr.var))+ 
  theme_classic(base_size = 18)+
  geom_smooth(method=lm , color="black", fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),
              se=TRUE,level=0.95) + 
  geom_jitter(aes(fill=duplicate),pch=21,size=5,width=6)+
  scale_fill_manual(name='Cell turnover estimate',values=c('white','black'))+ 
  ylab('Aging magnitude')+xlab('mean(24m)-mean(3m)')+
  ggtitle(paste0('Pearson\'s r = ',round(pearson.cor.value,3),
                 ', P value=',format.pval(pearson.pval,digits=2) ))+
  scale_y_continuous(limits=c(-1.7,1.6),labels=scaleFUN)
p.gcl<-p.gcl+ylab('Change in transcriptome variability (1-GCL)')+#xlab('mean(24m)-mean(3m)')+
  xlab(expression(paste(Delta,' age (day)'))) #xlab('Change in cell age (day)')

p.cv<-p.cv+ylab(expression(paste('Change in transcriptome variability (',beta[age],')')))+
  xlab(expression(paste(Delta,' age (day)')))+
  theme(legend.position="none")

pdf(out.fig,useDingbats = T,width = 7,height = 18)
grid.arrange(p.cv+xlab(''),p.rho+theme(legend.position = "none")+xlab(''),
             p.gcl+theme(legend.position = "none")+xlab(''),ncol=1)
dev.off()

