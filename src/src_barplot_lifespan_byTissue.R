
library(ggplot2);library(gridExtra)
library(tidyverse);
library(viridis)
######################################################################
## read in mouse turnover rate data
cell.lifespan=data.table::fread('../src/Dataset_S1.txt')
cell.lifespan$Ron_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`

dim(cell.lifespan) #39
tmp=cell.lifespan[!duplicated(cell.lifespan$Ron_tc),]
dim(tmp); #unique 22 tc
table(tmp$species) #1 human and 20 rodents

tissue=unlist(lapply(strsplit(cell.lifespan$`cell type annotation in TMS`,':'),'[',1))
cell.type=unlist(lapply(strsplit(cell.lifespan$`cell type annotation in TMS`,':'),'[',2))

cell.lifespan$tissue=tissue;
cell.lifespan$cell.type=cell.type

cell.lifespan[grep('Brain',cell.lifespan$tissue,ignore.case = T),]$tissue='Brain'
length(unique(cell.lifespan$tissue)) #17 tissue


x=names(which(table(cell.lifespan$Ron_tc)>1))
length(x) 
cell.lifespan$duplicate='tissue-specific estimate';
cell.lifespan[cell.lifespan$Ron_tc %in% x,]$duplicate='non tissue-specific estimate';
table(cell.lifespan$duplicate) #24 dups + 15 ts
length(unique(cell.lifespan[cell.lifespan$duplicate=="non tissue-specific estimate" ,]$Ron_tc)) #24dups => 7 unique tc

cell.type[grep('endothelial',cell.type,ignore.case = T)]
cell.type[grep('endothelial',cell.type,ignore.case = T)]='endothelial cells'
cell.lifespan$cell.identity=cell.type

table(cell.lifespan$duplicate)
####################
## plot by tissue
colnames(cell.lifespan)
df=cell.lifespan %>% arrange(tissue,lifespan)


hj=rep(-0.25,nrow(cell.lifespan))
hj[7]=-0.09
hj[15]=-0.6

df$lifespan2=df$lifespan;
df[df$lifespan>10,]$lifespan2=floor(df[df$lifespan>10,]$lifespan)
plot.grid<-ggplot(df,aes(x=cell.type,y=lifespan))+
    #geom_bar(aes(fill=Ron_tc),width=0.2,stat='identity')+
    facet_grid(tissue~.,scales = "free_y",space = "free")+
    geom_bar(width=0.1,stat='identity')+
    geom_point(aes(fill=duplicate),shape=21,size=2.5)+scale_fill_manual(values=c('white','black'))+
    geom_text(aes(label=round(df$lifespan2,1)), vjust=0.5,hjust=hj,size=4.5) +
    coord_flip()+
    scale_y_continuous(limits =c(0,2e6),trans = scales::pseudo_log_trans(0.1, 1e6),
                       breaks=c(0.1,1,1e1,1e2,1e3,1e4,1e5,1e6,1e6))+
    theme_classic(base_size = 15)+
    theme(legend.position = 'None')+
    ylab('Lifespan (day)')+xlab('Tissue:cell type combination in TMS')
print(plot.grid)

pdf('barplot_lifespan_byTissue.pdf',useDingbats = T,width = 10,height = 10)  
print(plot.grid+xlab('')+ylab('')+
        theme(strip.text.y = element_text(size = 10, colour = "black", angle = 0))
      )
dev.off()

