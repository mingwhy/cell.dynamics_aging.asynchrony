#install.packages('egg', dependencies = TRUE)
library(egg)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(ggpointdensity)
library(viridis)
library(patchwork); #for plot_annotation
options(stringsAsFactors = F)
library(org.Mm.eg.db,verbose=F,quietly=T)
library(GO.db);
mycols=c("#999999", "#E69F00");

cor.coeffs.df=readRDS('cor.coeffs.df.rds')
dim(cor.coeffs.df) 
cor.coeffs.df$tc.size

plot.sig=cor.coeffs.df[cor.coeffs.df$FDR<0.05 & cor.coeffs.df$tc.size>10,]
sort(plot.sig$FDR) #<0.05
sort(abs(plot.sig$Spearman.rho)) #>0.34
nrow(plot.sig) #10

df.all=readRDS('GOscore_tc.rds')
df<-lapply(df.all,function(x){
  x1=x[,grep('Ron_tc|GO:|age',colnames(x))]
  x2=reshape2::melt(x1)
  colnames(x2)=c('Ron_tc','age','GOid','score')
  x3=merge(x2,x[,c('Ron_tc','lifespan')])
  x3
})
df=as.data.frame(Reduce(`rbind`,df))
head(df)

df=df[df$GOid %in% plot.sig$GOid,]
df$age=factor(df$age,levels=c('3m','24m'))

length(unique(df$GOid))
sum(unique(df$GOid) %in% cor.coeffs.df$GOid)
df$GOid_age=paste(df$GOid,df$age)
cor.coeffs.df$GOid_age=paste(cor.coeffs.df$GOid,cor.coeffs.df$age)

df2=merge(df,cor.coeffs.df[,c(1,2,4,6,8,9)])
dim(df2)
dim(df)
df2$GOname=paste(df2$GOid,df2$GOterm,sep='\n')
df2$sig.or.not='not'
df2[df2$FDR<0.05,]$sig.or.not='yes'

summary(abs(df2[df2$sig.or.not=='yes',]$Spearman.rho))

df2$GOid=droplevels(df2$GOid)
unique(df2$GOid)
df2$GOid=factor(df2$GOid,levels=c('GO:0051131','GO:0061077','GO:0061684','GO:0006914','GO:0000422',
                                  'GO:0006281','GO:0043504'))
df2$GOname=factor(df2$GOname,unique(df2[order(df2$GOid),]$GOname))

p=ggplot(df2,aes(x=lifespan,y=score,col=age,fill=sig.or.not))+
  facet_grid(GOname~age,scales = "free_y")+
  geom_point(size=2,shape=21,col='black')+theme_classic(base_size = 12)+
  scale_x_log10()+scale_color_manual(values=mycols)+
  scale_fill_manual(values=c('white','black'))+
  #scale_fill_manual(values=mycols)+
  ylab('GO term expression')+xlab('Cell lifespan (day)')+
  theme(#panel.spacing = unit(0, "lines"),
        panel.spacing.x = unit(4, "mm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
        axis.text.y = element_text(size=8),
        strip.text.x = element_text(size = 12, colour = "black", angle = 0),
        strip.text.y = element_text(size = 10, colour = "black", angle = 0),
        legend.position = 'none')
p

tmp=df2 %>% arrange(.,GOid,age)
#df3=df2[!duplicated(df2$GOid_age),]
df3=tmp[!duplicated(tmp$GOid_age),]
df3$my_tag=paste0('FDR = ',sprintf('%.3f',df3$FDR))
df3$my_tag
df3$score

p1=p+geom_text(x = 1.2, y = Inf, size=2.8,col='black',aes(label = my_tag), 
               #data = df3,hjust=-0.4,vjust=2)
               data = df3,hjust=-1,vjust=3)
p1

pdf('MainFig3_GO_results.pdf',useDingbats = T,height = 8.5,width = 7.5)
print(p1)
dev.off()

