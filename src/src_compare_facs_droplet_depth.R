library(tidyverse);
library(ggplot2);
library(gridExtra)

######################################################################
## save cell.meta file for droplet
if(!file.exists('droplet_ncell_per.age_per.tc_per.sex_raw.txt')){
  
  library(zellkonverter)
  library(SingleCellExperiment)
  library(scater);library(scran)
  library(ggplot2);library(gridExtra)
  library(SummarizedExperiment)
  library(ggpointdensity)
  library(viridis)
  library(patchwork); #for plot_annotation
       
  inp.sce<-readH5AD('~/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-droplet-official-raw-obj.h5ad');  
  
  colData(inp.sce)$tissue_cell.type=paste(inp.sce$tissue,inp.sce$cell_ontology_class,sep=':')
  df=as.data.frame(colData(inp.sce))

  summary(df$n_genes) #a lot of NA
  assayNames(inp.sce)
  m=assay(inp.sce,'X')
  n_gene=Matrix::colSums(m>0)
  summary(n_gene) #min: 250
  df$n_genes=n_gene
  
  n_counts=Matrix::colSums(m)
  summary(n_counts) #min: 2500
  df$n_counts=n_counts
  
  data.table::fwrite(df,'droplet_ncell_per.age_per.tc_per.sex_raw.txt')
}

######################################################################
## save cell.meta file for facs
if(!file.exists('facs_ncell_per.age_per.tc_per.sex_raw.txt')){
  
  library(zellkonverter)
  library(SingleCellExperiment)
  library(scater);library(scran)
  library(ggplot2);library(gridExtra)
  library(SummarizedExperiment)
  library(ggpointdensity)
  library(viridis)
  library(patchwork); #for plot_annotation
  
  inp.sce<-readH5AD('~/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');      
  
  colData(inp.sce)$tissue_cell.type=paste(inp.sce$tissue,inp.sce$cell_ontology_class,sep=':')
  df=as.data.frame(colData(inp.sce))
  summary(df$n_genes)  #min 499
  summary(df$n_counts) #min 5002
  
  df.facs=df;
  df.facs=df.facs[df.facs$age!='21m',]
  df.facs %>% group_by(sex,age) %>% summarise(n=length(unique(mouse.id)))  

  data.table::fwrite(df,'facs_ncell_per.age_per.tc_per.sex_raw.txt')
}

######################################################################
## make plot
df.facs=data.table::fread('facs_ncell_per.age_per.tc_per.sex_raw.txt')
df.facs$platform='FACS'
colnames(df.facs)
df.facs=df.facs[df.facs$age!='21m',]
df.facs %>% group_by(sex,age) %>% summarise(n=length(unique(mouse.id)))
df.facs=df.facs[-grep('\\/',df.facs$mouse.id),]

df.droplet=data.table::fread('droplet_ncell_per.age_per.tc_per.sex_raw.txt')
colnames(df.droplet)
df.droplet$platform='Droplet'
df.droplet %>% group_by(sex,age) %>% summarise(n=length(unique(mouse.id)))
df.droplet %>% group_by(mouse.id,age) %>% summarise(ncell=n())

x=intersect(colnames(df.facs),colnames(df.droplet))
df=rbind(df.facs[,..x],df.droplet[,..x])
df$platform=factor(df$platform,levels=c('FACS','Droplet'))
df$sex=factor(df$sex,levels=c('male','female'))


df[df$platform=='FACS',] %>% group_by(sex,age) %>% summarise(n=n())
ages=paste0(sort(as.numeric(gsub('m','',unique(df$age)))),'m')
df$age=factor(df$age,levels=ages)
p1=ggplot(df,aes(x=age,y=n_counts,col=platform))+
  facet_wrap(sex~.,scales = "free_x")+geom_violin()+
  scale_y_log10()+theme_classic(base_size = 12)+#geom_boxplot(width=.1)+
  stat_summary(fun = median, geom = "point",
               position = position_dodge(0.9))+
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  ylab("Number of detected transcripts per cell")

x=df %>% group_by(sex,age,platform) %>% summarise(n.mouse=length(unique(mouse.id)))
x$platform=factor(x$platform,levels=levels(df$platform))

p2=ggplot(x,aes(x=age,y=n.mouse,fill=platform))+
  facet_wrap(sex~.,scales = "free_x")+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  ylab('Number of mice per condition')+
  geom_bar(stat='identity',position='dodge')+theme_classic(base_size = 12)

df1=df[df$n_genes>=500 & df$n_counts>=5000,]
x=df1  %>% group_by(platform,sex,age,tissue_cell.type) %>%
  summarise(ncell=n()) %>% filter(ncell>=40)
x1=x  %>% group_by(platform,sex,age) %>%  summarise(ntc=length(unique(tissue_cell.type))) 

p3=ggplot(x1,aes(x=age,y=ntc,fill=platform))+
  facet_wrap(sex~.,scales = "free_x")+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  ylab('Number of tissue:cell type combinations')+
  geom_bar(stat='identity',position='dodge')+theme_classic(base_size = 12)

library(egg)

pdf('compare_facs_droplet.pdf',useDingbats = T,height = 10)
ggarrange(p2,p1,p3,ncol=1,
             labels = c("A", "B","C"));
dev.off()


