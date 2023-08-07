
library(tidyverse);
library(ggplot2);library(gridExtra)
library(viridis)
library(ggplot2) 
library(grid)
library(gridExtra) 

## plot all tcs
df2=data.table::fread('change.in.expression_change.in.age.txt');
table(df2$cell.identity)
table(df2$duplicate)
df3=df2 %>% group_by(Ron_tc) %>% mutate('change.in.expr.var'=mean(change.in.expr))
df4=df3[!duplicated(df3$Ron_tc),] 
dim(df4) #22

# include duplciates
test.out=cor.test(df3$mean.diff,df3$change.in.expr.var,method='spearman')
(cor2=test.out$estimate)
(pval2=test.out$p.value)
test.out=cor.test(df3$mean.diff,df3$change.in.expr.var,method='pearson')
(cor2.pearson=test.out$estimate) 
(pval2.pearson=test.out$p.value)

p1=ggplot(df3,aes(x=mean.diff,y=change.in.expr))+
  geom_point(aes(shape=duplicate),size=3)+
  theme_classic(base_size = 15)+ggtitle(paste0('All tissue:cell type combinations\n',
                                               'Pearson\'s r = ',round(cor2.pearson,3),
                                               ', P value=',format.pval(pval2.pearson,digits=2)))+
  
  geom_smooth(data=df3,aes(x=mean.diff,y=change.in.expr.var),method=lm , color="black", 
              fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),se=TRUE,level=0.95)+ 
  scale_shape_manual(name='Cell lifespan estimate',values=c(1,2))
p1

## endothelial cells only
df3=df2[df2$cell.identity=='endothelial cells',]
dim(df3) #12
table(df3$duplicate)


test.out=cor.test(df3$change.in.expr,df3$mean.diff,method='spearman')
(cor2=test.out$estimate) 
(pval2=test.out$p.value) 
test.out=cor.test(df3$change.in.expr,df3$mean.diff,method='pearson')
(cor2.pearson=test.out$estimate) 
(pval2.pearson=test.out$p.value) 


p2=ggplot(df3,aes(x=mean.diff,y=change.in.expr))+  
  geom_point(aes(shape=duplicate), size=3)+
  theme_classic(base_size = 15)+ggtitle(paste0('Endothelial cells\n',
                                               'Pearson\'s r = ',round(cor2.pearson,3),
                                               ', P value=',format.pval(pval2.pearson,digits=2)))+
  geom_smooth(data=df4,aes(x=mean.diff,y=change.in.expr.var),method=lm , color="black", 
              fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),se=TRUE,level=0.95)+ 
  scale_shape_manual(name='Cell lifespan estimate',values=c(1,2))
p2

library(egg)
pdf('endothelial_cells_cell.age.difference_shape.pdf',useDingbats = T,width = 10,height = 12)
ggarrange(
  p1+ylab(expression(paste('Change in transcriptome variability (',beta[age],')')))+
    #xlab(expression(paste(Delta,' age (day)')))+ 
    xlab('Change in mean cell age (day)')+
    theme(plot.title = element_text(size = 15, face = "bold")),
  p2+ylab(expression(paste('Change in transcriptome variability (',beta[age],')')))+
    ylim(c(-0.01,0.03))+#xlim(c(0,420))+
    #xlab(expression(paste(Delta,' age (day)')))+
    xlab('Change in mean cell age (day)')+
  theme(plot.title = element_text(size = 15, face = "bold"))+
    annotate('text',x=280,y=0.025,label='Regression line from Figure 2B in the maintext.'),
  ncol=1,labels = c("A", "B"));    
dev.off()
