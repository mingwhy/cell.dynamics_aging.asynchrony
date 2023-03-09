
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
df3=df2 %>% group_by(human_tc) %>% mutate('change.in.expr.var'=mean(change.in.expr))
df4=df3[!duplicated(df3$human_tc),]
dim(df4) #21

# include duplciates
test.out=cor.test(df3$mean.diff,df3$change.in.expr.var,method='spearman')
(cor2=test.out$estimate)
(pval2=test.out$p.value)
test.out=cor.test(df3$mean.diff,df3$change.in.expr.var,method='pearson')
(cor2.pearson=test.out$estimate) 
(pval2.pearson=test.out$p.value)

p1=ggplot(df3,aes(x=mean.diff,y=change.in.expr))+
  geom_point(aes(shape=duplicate),size=3)+
  theme_classic(base_size = 15)+ggtitle(paste0('All cell types\n',
                                               #'Spearman\'s rho = ',round(cor2,3),
                                               #', P value=',round(pval2,8),  #', P value=',sprintf('%.4f',pval1)))+
                                               'Pearson\'s r = ',round(cor2.pearson,3),
                                               ', P value=',round(pval2.pearson,8)))+
  geom_smooth(data=df3,aes(x=mean.diff,y=change.in.expr.var),method=lm , color="black", 
              fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),se=TRUE,level=0.90)+ 
  scale_shape_manual(name='Cell lifespan estimate',values=c(1,2))
p1

## endothelial cells only
df3=df2[df2$cell.identity=='endothelial cells',]
dim(df3) #12
table(df3$duplicate)

test.out=cor.test(df3$change.in.expr,df3$mean.diff,method='spearman')
(cor2=test.out$estimate) #0.6074271  
(pval2=test.out$p.value) #0.03618097
test.out=cor.test(df3$change.in.expr,df3$mean.diff,method='pearson')
(cor2.pearson=test.out$estimate) #0.6541237
(pval2.pearson=test.out$p.value) #0.02102181


p2=ggplot(df3,aes(x=mean.diff,y=change.in.expr))+
  geom_point(aes(shape=duplicate), size=3)+
  theme_classic(base_size = 15)+ggtitle(paste0('Endothelial cells\n',
                                               #'Spearman\'s rho = ',round(cor2,3),
                                               #', P value=',round(pval2,6),  #', P value=',sprintf('%.4f',pval1)))+
                                               'Pearson\'s r = ',round(cor2.pearson,3),
                                               ', P value=',round(pval2.pearson,6)))+
  geom_smooth(data=df4,aes(x=mean.diff,y=change.in.expr.var),method=lm , color="black", 
              fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),se=TRUE,level=0.90)+ 
  scale_shape_manual(name='Cell lifespan estimate',values=c(1,2))
p2

pdf('endothelial_cells_cell.age.difference_shape.pdf',useDingbats = T,width = 10,height = 12)
grid.arrange(
  p1+ylab(expression(paste('Change in transcriptome variability (',beta[age],')')))+
    xlab('Change in cell age (day)')+theme(plot.title = element_text(size = 15, face = "bold")),
  p2+ylab(expression(paste('Change in transcriptome variability (',beta[age],')')))+
    ylim(c(-0.01,0.03))+#xlim(c(0,420))+
    xlab('Change in cell age (day)')+#theme(legend.position = 'none',plot.title = element_text(size = 15, face = "bold"))
  theme(plot.title = element_text(size = 15, face = "bold"))+
    annotate('text',x=280,y=0.025,label='Regression line from Figure 2 in the maintext.')
    )
dev.off()
