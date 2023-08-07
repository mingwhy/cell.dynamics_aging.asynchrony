
# GO.db vignettes, basic GO usage
# https://www.bioconductor.org/packages/release/bioc/vignettes/annotate/inst/doc/GOusage.pdf

options(stringsAsFactors = F)
library(org.Mm.eg.db,verbose=F,quietly=T)
library(GO.db);
packageVersion('org.Mm.eg.db') #"3.14.0"
packageVersion('GO.db') #"3.14.0"

## extract genes for each GO term
all.genes<-keys(org.Mm.eg.db,"ENSEMBL")
length(all.genes) #32941 genes
length(unique(all.genes)) #32941 genes
fb2go=select(org.Mm.eg.db,keys=all.genes,keytype = 'ENSEMBL',
             columns = c('GO'))
head(fb2go) #pay attention to EVIDENCE column
table(fb2go[!duplicated(fb2go$GO),]$ONTOLOGY)
#BP    CC    MF 
#12545  1751  4417  
###############################################
## select BP
sub.fb2go=fb2go[fb2go$ONTOLOGY=='BP',]
dim(sub.fb2go) #185351     4
x<-split(sub.fb2go$ENSEMBL,sub.fb2go$GO)
length(x) #12545 BP terms
go2fb<-sapply(x,function(x){ unique(x)}) #assign genes to each GO term

# map ENSEMBL to SYMBOL
length(go2fb) #12545 BP GO terms
go2fb[[1]]
# change FBgn to symbol
go2symbol=lapply(go2fb,function(i){
  AnnotationDbi::select(org.Mm.eg.db,keys=i,keytype='ENSEMBL',c('SYMBOL'))
})
names(go2symbol)
saveRDS(go2symbol,'Mmus_go2SYMBOL_BP.rds')
###############################################
## select CC
sub.fb2go=fb2go[fb2go$ONTOLOGY=='CC',]
dim(sub.fb2go) #185351     4
x<-split(sub.fb2go$ENSEMBL,sub.fb2go$GO)
length(x) #12545 BP terms
go2fb<-sapply(x,function(x){ unique(x)}) #assign genes to each GO term

# map ENSEMBL to SYMBOL
length(go2fb) #12545 BP GO terms
go2fb[[1]]
# change FBgn to symbol
go2symbol=lapply(go2fb,function(i){
  AnnotationDbi::select(org.Mm.eg.db,keys=i,keytype='ENSEMBL',c('SYMBOL'))
})
names(go2symbol)
saveRDS(go2symbol,'Mmus_go2SYMBOL_CC.rds')
###############################################
## select MF
sub.fb2go=fb2go[fb2go$ONTOLOGY=='MF',]
dim(sub.fb2go) #185351     4
x<-split(sub.fb2go$ENSEMBL,sub.fb2go$GO)
length(x) #12545 BP terms
go2fb<-sapply(x,function(x){ unique(x)}) #assign genes to each GO term

# map ENSEMBL to SYMBOL
length(go2fb) #12545 BP GO terms
go2fb[[1]]
# change FBgn to symbol
go2symbol=lapply(go2fb,function(i){
  AnnotationDbi::select(org.Mm.eg.db,keys=i,keytype='ENSEMBL',c('SYMBOL'))
})
names(go2symbol)
saveRDS(go2symbol,'Mmus_go2SYMBOL_MF.rds')

