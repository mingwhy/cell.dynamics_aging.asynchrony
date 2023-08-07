
library(org.Mm.eg.db,verbose=F,quietly=T)
library(GO.db);
options(stringsAsFactors = F)
library(igraph)
library(Matrix);
library(ggplot2);library(gridExtra)
library(ggpubr)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
library(UpSetR)


####################################
if(!file.exists('GOid_GOterm.rds')){
  slim.go2fb=readRDS('../src/Mmus_go2SYMBOL_BP.rds')
  length(slim.go2fb) #1882 or 12545 (all BP GOterms)
  all.go.id=names(slim.go2fb)
  goterms=Term(GOTERM)
  x=lapply(all.go.id,function(i) GOTERM[[i]]@Term)
  df.go=data.frame(GOid=all.go.id,GOterm=unlist(x))
  dim(df.go) #12545
  saveRDS(df.go,'GOid_GOterm.rds')
}

df.go=readRDS('GOid_GOterm.rds')
dim(df.go) #12545
df.go.pick=df.go[grep('chaperone|DNA repair|autophagy of mitochondrion',ignore.case = t,df.go$GOterm),]
df.go.pick=df.go[grep('chaperone|DNA repair|\\bautophagy\\b',ignore.case = t,df.go$GOterm),]
dim(df.go.pick) #41

rownames(df.go.pick)=df.go.pick$GOid
dim(df.go.pick) #41 x 2


slim.go2fb=readRDS('../src/Mmus_go2SYMBOL_BP.rds')
length(slim.go2fb) #12545
ngene.per.go=sapply(slim.go2fb,nrow)
sum(ngene.per.go>=6) #4516

go2gene=slim.go2fb[names(slim.go2fb) %in% df.go.pick$GOid]
length(go2gene)
x=sort(sapply(go2gene,nrow)) 

#keep.go=names(x[x>=5]) #13
keep.go=names(x[x>=6]) 
length(keep.go) #17

df.go.pick=df.go.pick[keep.go,]
go2gene=go2gene[keep.go]
length(go2gene) #17
#https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
input.list=lapply(go2gene,'[[',2)
names(input.list)
input.list1=input.list
sum(names(input.list1)==df.go.pick$GOid)
names(input.list1)=paste0(df.go.pick$GOid,'\n',df.go.pick$GOterm)

m = make_comb_mat(input.list1)
ss = set_size(m)
cs = comb_size(m)
#https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
p1=UpSet(m, 
      top_annotation = upset_top_annotation(m, add_numbers = TRUE),
         #set_order = order(ss),
         #comb_order = order(comb_degree(m), -cs),
      right_annotation = upset_right_annotation(m, add_numbers = TRUE),
      left_annotation = rowAnnotation(
        set_name = anno_text(set_name(m), 
                             location = 0.5, 
                             just = "center",
                             width = max_text_width(set_name(m)) + unit(1, "mm"))
       ), 
     
      #right_annotation = NULL,
      show_row_names = FALSE)
p1

pdf('upset.plot_GOterms.pdf',useDingbats = T,width = 12,height = 9)
print(p1)
dev.off()

