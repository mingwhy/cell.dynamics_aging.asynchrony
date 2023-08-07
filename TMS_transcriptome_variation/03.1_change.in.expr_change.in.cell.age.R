
library(tidyverse);
library(ggplot2);library(gridExtra)
library(viridis)
library(ggplot2) 
library(grid)
library(gridExtra) 
#https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
scaleFUN <- function(x) sprintf("%.3f", x)

age.groups=c(3,24)

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
para_df=readRDS(paste0('../src/change.in.cell.age_',age.groups[1],'_vs_',age.groups[2],'.rds'));
#para_df=readRDS('../src/change.in.cell.age_3_vs_24.rds')
unique(para_df$month)

age1=age.groups[1]; age2=age.groups[2];
para_df2<-para_df %>% group_by(Ron_tc) %>% 
  summarise(mean.ratio= mean[month== age2 ] / mean[month== age1],
            mean.diff= mean[month==age2 ] - mean[month==age1],
            sd.ratio= var[month==age2]^0.5 / var[month==age1]^0.5,
            sd.diff= var[month==age2]^0.5 - var[month==age1]^0.5);

intersect(colnames(cell.lifespan),colnames(para_df2)) #`cell type annotation in Sender and Milo 2021`

df.model=merge(para_df2,cell.lifespan)
colnames(df.model)
head(df.model)
dim(df.model) #39

##############################################################
## gene-level metric
res.gene=readRDS('LMM_cv2_out.rds')

dim(res.gene) 
head(res.gene)
colnames(res.gene)[1]='cell.type'

res.gene=res.gene[res.gene$cell.type %in% df.model$`cell type annotation in TMS`,]
res.gene=res.gene[order(res.gene$cell.type),]

head(res.gene)
summary(p.adjust(res.gene$pvalue,method='BH')) 

df=res.gene
df$change.in.expr=df$beta
dim(df)

######################################################################
## integrate expr.change and age.change 
sum(df$cell.type %in% df.model$`cell type annotation in TMS`) #

df.model$cell.type=df.model$`cell type annotation in TMS`;
intersect(colnames(df.model),colnames(df)) #cell.type

df2=merge(df.model[,c('cell.type','cell.identity',"duplicate" ,'Ron_tc',
                      'tissues','lifespan',
                      'mean.diff','mean.ratio','sd.diff','sd.ratio')],
          df[,c('cell.type','change.in.expr')]) 
# change lifespan into turnover
df2$turnover=1-exp(-1/df2$lifespan)

colnames(df2)
df2=df2[order(df2$lifespan),]

# save a table
colnames(df2)
df3=df2[,c(1,2,3,4, 6,7,11)]
data.table::fwrite(df3,'change.in.expression_change.in.age.txt')

##################################################
## calculate correlation
table(df2$duplicate)
df3=df2 %>% group_by(Ron_tc) %>% mutate('change.in.expr.var'=mean(change.in.expr))
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

p2=ggplot(df4,aes(x=mean.diff,y=change.in.expr.var))+
  geom_jitter(aes(fill=duplicate),pch=21,size=5,width=6)+
  theme_classic(base_size = 22)+
  ggtitle(paste0('Unique cell types\n',
            'Spearman\'s rho = ',round(cor2,3),
                 ', P value=',format.pval(pval2,digits=2),  #', P value=',sprintf('%.4f',pval1)))+
                 '\nPearson\'s r = ',round(cor2.pearson,3),
                 ', P value=',format.pval(pval2.pearson,digits=2)))+
  geom_smooth(method=lm , color="black", fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),
              se=TRUE,level=0.95) + 
  scale_fill_manual(name='Cell lifespan estimate',values=c('NA','black'))+
  scale_y_continuous(labels=scaleFUN)
p2


####################################################################################
## include duplicates
test.out=cor.test(df3$mean.diff,df3$change.in.expr.var,method='spearman')
(cor2=test.out$estimate)
(pval2=test.out$p.value)
test.out=cor.test(df3$mean.diff,df3$change.in.expr.var,method='pearson')
(cor2.pearson=test.out$estimate) 
(pval2.pearson=test.out$p.value)

p1=ggplot(df3,aes(x=mean.diff,y=change.in.expr))+  
  geom_point(aes(shape=duplicate),size=3)+
  theme_classic(base_size = 22)+ggtitle(paste0('All tissue:cell type combinations\n',
                                               'Spearman\'s rho = ',round(cor2,3),
                                               ', P value=',format.pval(pval2,digits=2),  #', P value=',sprintf('%.4f',pval1)))+
                                               '\nPearson\'s r = ',round(cor2.pearson,3),
                                               ', P value=',format.pval(pval2.pearson,digits=2)))+
  
  geom_smooth(data=df3,aes(x=mean.diff,y=change.in.expr.var),method=lm , color="black", 
              fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),se=TRUE,level=0.95)+ 
  scale_shape_manual(name='Cell lifespan estimate',values=c(1,2))+
  scale_y_continuous(labels=scaleFUN)


pdf('MainFig2_change.in.transcriptome_change.in.age_ncell15.pdf',useDingbats = T,width = 11,height = 7)

print(  p2+ylab(expression(paste('Change in transcriptome variability (',beta[age],')')))+
    xlab(expression(paste(Delta,' age (day)')))+
    theme(plot.title = element_text(size = 0, face = "bold"), legend.position="right",
          legend.title=element_text(size=15),
          legend.text=element_text(size=14))
)
dev.off()
  

pdf('change.in.transcriptome_change.in.age_ncell15.pdf',useDingbats = T,width = 12,height = 14)
grid.arrange(
  p2+ylab(expression(paste('Change in transcriptome variability (',beta[age],')')))+
    xlab(expression(paste(Delta,' age (day)')))+
    theme(plot.title = element_text(size = 15, face = "bold"), legend.position="right",
          legend.title=element_text(size=15),
          legend.text=element_text(size=14)),
  p1+ylab(expression(paste('Change in transcriptome variability (',beta[age],')')))+
    xlab(expression(paste(Delta,' age (day)')))+ 
    theme(plot.title = element_text(size = 15, face = "bold"), legend.position="right",
          legend.title=element_text(size=15),
          legend.text=element_text(size=14)),
  ncol=1);
dev.off()



