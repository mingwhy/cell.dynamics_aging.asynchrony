
library(tidyverse);
library(ggplot2);library(gridExtra)
library(viridis)
library(ggplot2) 
library(grid)
library(gridExtra) 
library(ggpmisc) #https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
#https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
scaleFUN <- function(x) sprintf("%.3f", x)
mycols=c("#999999", "#E69F00");
mycols=c("#999999",  "#56B4E9","#E69F00");

age.groups=c(3,18,24)

cell.lifespan=data.table::fread('../src/Dataset_S1.txt')

cell.lifespan$Ron_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`
dim(cell.lifespan)

tissues=sapply(strsplit(cell.lifespan$`cell type annotation in TMS`,'\\:'),'[[',1)
tcs=sapply(strsplit(cell.lifespan$Ron_tc,'\\:'),'[[',2)

x=names(which(table(cell.lifespan$Ron_tc)>1))
length(x) 
cell.lifespan$duplicate='tissue-specific estimate';
cell.lifespan[cell.lifespan$Ron_tc %in% x,]$duplicate='non tissue-specific estimate';
table(cell.lifespan$duplicate) #24 dups + 15 ts
length(unique(cell.lifespan[cell.lifespan$duplicate=="non tissue-specific estimate" ,]$Ron_tc)) #24dups => 7 unique tc

tcs[grep('endothelial',tcs,ignore.case = T)]
tcs[grep('endothelial',tcs,ignore.case = T)]='endothelial cells'
cell.lifespan$cell.identity=tcs
cell.lifespan$tissues=tissues;

######################################################################
## read in mean.diff and sd.diff in theoretical model
para_df1=readRDS(paste0('../src/change.in.cell.age_',age.groups[1],'_vs_',age.groups[2],'.rds'));
para_df2=readRDS(paste0('../src/change.in.cell.age_',age.groups[1],'_vs_',age.groups[3],'.rds'));

para_df=as.data.frame(rbind(para_df1[para_df1$month==3,],para_df1[para_df1$month==18,],para_df2[para_df2$month==24,]))
dim(para_df)

intersect(colnames(cell.lifespan),colnames(para_df)) #`cell type annotation in Sender and Milo 2021`

df.model=merge(para_df,cell.lifespan,all.x=T)

colnames(df.model)
head(df.model)
dim(df.model) #39 * 3

##############################################################
## median CV per age group, 3 vs 18
DM.out=readRDS('../0708_TMS_male_3m_18m/MAST_nonDEgene_stats.rds')
tc.names=names(DM.out)
length(tc.names);

keep.cell.types=intersect(tc.names,cell.lifespan$`cell type annotation in TMS`)
length(keep.cell.types)
DM.out=DM.out[keep.cell.types]

res<-lapply(keep.cell.types,function(tc){
  all.ages=DM.out[[tc]]
  all.ages=all.ages[!is.na(all.ages$mean_scran),]
  
  df=all.ages[all.ages$age %in% paste0(age.groups,'m'),]
  if(length(unique(df$age))==1){return(NULL)}
  
  ngene=length(unique(df$gene))
  df$numeric.age=as.numeric(gsub('m','',df$age))
  
  df$log10_cv2_scran=log10(df$cv2_scran)
  x=df %>% group_by(numeric.age,mouse.id) %>% 
    summarise(CV.per.id=mean(log10_cv2_scran)) %>% 
    summarise(CV=mean(CV.per.id))
  x$tc=tc
  x
})
df.res1=as.data.frame(Reduce(`rbind`,res))


##############################################################
## median CV per age group, 3 vs 24
DM.out=readRDS('../0708_TMS_male_3m_24m/MAST_nonDEgene_stats.rds')
tc.names=names(DM.out)
length(tc.names);

keep.cell.types=intersect(tc.names,cell.lifespan$`cell type annotation in TMS`)
length(keep.cell.types)
DM.out=DM.out[keep.cell.types]

res<-lapply(keep.cell.types,function(tc){
  all.ages=DM.out[[tc]]
  all.ages=all.ages[!is.na(all.ages$mean_scran),]
  
  df=all.ages[all.ages$age %in% paste0(age.groups,'m'),]
  if(length(unique(df$age))==1){return(NULL)}
  
  #df=df[df$gene %in% names(which(table(df$gene)==length(age.groups))),] #expr in all age groups
  ngene=length(unique(df$gene))
  df$numeric.age=as.numeric(gsub('m','',df$age))
  
  df$log10_cv2_scran=log10(df$cv2_scran)
  x=df %>% group_by(numeric.age,mouse.id) %>% 
    summarise(CV.per.id=mean(log10_cv2_scran)) %>% 
    summarise(CV=mean(CV.per.id))
  x$tc=tc
  x
})
df.res2=as.data.frame(Reduce(`rbind`,res))
####################################################################################
## combine 3 age groups
head(df.res1)
head(df.res2)
df.res0=as.data.frame(rbind(df.res1,df.res2))
df.res=df.res0 %>% group_by(tc,numeric.age) %>% summarise(CV=mean(CV))
table(df.res$numeric.age)

df.res$tc_age=paste(df.res$tc,df.res$numeric.age,sep='-')
df.model$tc_age=paste(df.model$`cell type annotation in TMS`,df.model$month,sep='-')

df=merge(df.res,df.model[,c('Ron_tc','month','mean','duplicate','lifespan','tc_age')])

df$month=factor(df$month,levels=age.groups)

p1<-ggplot(df,aes(x=mean,y=CV,group=month,col=month))+
  geom_smooth(method=lm,level=0.99,fill = "grey",alpha=0.1)+
  theme_classic(base_size = 12)+
  xlab('Mean cell age')+ylab('Average log10(CV)')+#scale_x_log10()+
  geom_point()+scale_color_manual(name='Month',values=mycols)
p1


## took average for the same cell type
df.rmdup=df %>% group_by(numeric.age,Ron_tc) %>% mutate(CV.rm.dup=mean(CV))
df.rmdup2=df.rmdup[!duplicated(paste(df.rmdup$Ron_tc,df.rmdup$numeric.age)),]
dim(df.rmdup2)
table(df.rmdup2$numeric.age)

p3<-ggplot(df.rmdup2,aes(x=mean,y=CV.rm.dup,group=month,col=month))+
  geom_smooth(method=lm,level=0.99,fill = "grey",alpha=0.1)+
  theme_classic(base_size = 16)+
  xlab('Mean cell age (day)')+ylab('Average log10(CV)')+
  geom_point()+scale_color_manual(name='Month',values=mycols)+
  stat_poly_eq(use_label(c("R2","p")))
p3


pdf('mean.cell.age_obs.CV_male.pdf',useDingbats = T,width = 6,height = 5)
print(p3);
dev.off()

sort(unique(df.rmdup2$month))
fitted_models<-lapply(sort(unique(df.rmdup2$month)),function(i){
  tmp=df.rmdup2[df.rmdup2$month==i,]
  lm(tmp$CV.rm.dup~tmp$mean)
})
lapply(fitted_models,summary)

