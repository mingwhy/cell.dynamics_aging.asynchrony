# TMS Atlas info: https://github.com/czbiohub/tabula-muris-senis

# data source: https://figshare.com/projects/Tabula_Muris_Senis/64982
# folder specifications: 
# TMS_processedFiles: Processed files (to use with scanpy)
# TMS.gene.data_final: The data for the TMS gene analysis accompanying the paper "Mouse Aging Cell Atlas Analysis Reveals Global and Cell Type Specific Aging Signatures", eLife, 2021.
# TMS.Data.Objects: Tabula Muris Senis Data Objects

## read in h5ad file
library(zellkonverter)
library(SummarizedExperiment)
sce.facs<-readH5AD('~/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad'); 
sce.droplet<-readH5AD('~/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-droplet-official-raw-obj.h5ad');
sce.facs; #gene by cell: 22966 110824 
sce.droplet; #gene by cell: 20138 245389 

## sample meta information
cell.meta.facs=colData(sce.facs)
cell.meta.droplet=colData(sce.droplet)
table(cell.meta.facs$age,cell.meta.facs$sex)
table(cell.meta.droplet$age,cell.meta.droplet$sex)

## acquire gene id use biomart (https://www.biostars.org/p/301116/)
if(F){
  library(biomaRt)
  mouse= useEnsembl(version = 99, #this archived contains dn ds values
                    biomart = 'ENSEMBL_MART_ENSEMBL', 
                    dataset = 'mmusculus_gene_ensembl')
  
  listAttributes(mouse) 
  x=listAttributes(mouse) 
  x[grep('symbol',x$name),]
  mouse_genes <- getBM(attributes = c('ensembl_gene_id', 
                                      'mgi_symbol','description',
                                      'chromosome_name',  
                                      'start_position', 
                                      'end_position'),
                       mart = mouse)
  tail(mouse_genes)
  saveRDS(mouse_genes,'~/Documents/Data_mouse_aging_atlas/mmusculus_gene_ensemblv99.rds')
}
mouse_genes=readRDS('~/Documents/Data_mouse_aging_atlas/mmusculus_gene_ensemblv99.rds')
sum('0610009B22Rik' %in% mouse_genes$mgi_symbol )

## generate gene id mapping table for facs data
if(F){
  mouse_genes=readRDS('~/Documents/Data_mouse_aging_atlas/mmusculus_gene_ensemblv99.rds')
  
  genes=rownames(sce.facs)
  length(genes) # 22966
  sum(genes %in% mouse_genes$mgi_symbol) #fac:20449 
  
  mouse_genes2=mouse_genes[mouse_genes$mgi_symbol %in% genes,]
  tmp=mouse_genes2[duplicated(mouse_genes2$mgi_symbol),]
  
  ok1=mouse_genes2[!mouse_genes2$mgi_symbol %in% tmp$mgi_symbol,]
  length(unique(ok1$mgi_symbol));dim(ok1) #20158 genes
  dups=mouse_genes2[mouse_genes2$mgi_symbol %in% tmp$mgi_symbol,]
  dups=dups[order(dups$mgi_symbol),]
  head(dups)
  length(unique(dups$mgi_symbol)) #291 genes
  
  dups[-grep('^\\d+$|X',dups$chromosome_name),]
  dups=dups[grep('^\\d+$|X',dups$chromosome_name),]
  length(unique(dups$mgi_symbol)) #291 genes
  
  tmp=dups[duplicated(dups$mgi_symbol),]
  ok2=dups[!dups$mgi_symbol %in% tmp$mgi_symbol,]
  length(unique(ok2$mgi_symbol));dim(ok2) #258 genes are ok now.
  dups=dups[dups$mgi_symbol %in% tmp$mgi_symbol,]
  View(dups); #choose one with longer gene description
  dim(dups)
  dups[dups$mgi_symbol=='Ndor1',]
  ok3=c()
  for(gene in unique(dups$mgi_symbol)){
    x=dups[dups$mgi_symbol==gene,]
    ok3=rbind(ok3,x[order(x$ensembl_gene_id)[1],])
  }
  dim(ok3) #33 genes
  ok=rbind(ok1,ok2,ok3)
  dim(ok) #20449 genes
  id.mapping=ok;
  sum(duplicated(id.mapping$mgi_symbol)); #0
  data.table::fwrite(id.mapping,'fac_20449genes_id.mapping.txt')
}
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt')
dim(id.mapping) #20449 genes

## generate gene id mapping table for droplet data
if(F){
  mouse_genes=readRDS('~/Documents/Data_mouse_aging_atlas/mmusculus_gene_ensemblv99.rds')
  
  genes=rownames(sce.droplet)
  length(genes) # 20138
  sum(genes %in% mouse_genes$mgi_symbol) # droplet: 17864
  
  mouse_genes2=mouse_genes[mouse_genes$mgi_symbol %in% genes,]
  tmp=mouse_genes2[duplicated(mouse_genes2$mgi_symbol),]
  
  ok1=mouse_genes2[!mouse_genes2$mgi_symbol %in% tmp$mgi_symbol,]
  length(unique(ok1$mgi_symbol));dim(ok1) #17608 genes
  dups=mouse_genes2[mouse_genes2$mgi_symbol %in% tmp$mgi_symbol,]
  dups=dups[order(dups$mgi_symbol),]
  head(dups)
  length(unique(dups$mgi_symbol)) #256 genes
  
  dups[-grep('^\\d+$|X',dups$chromosome_name),]
  dups=dups[grep('^\\d+$|X',dups$chromosome_name),]
  length(unique(dups$mgi_symbol)) #256 genes
  
  tmp=dups[duplicated(dups$mgi_symbol),]
  ok2=dups[!dups$mgi_symbol %in% tmp$mgi_symbol,]
  length(unique(ok2$mgi_symbol));dim(ok2) #228 genes are ok now.
  dups=dups[dups$mgi_symbol %in% tmp$mgi_symbol,]
  View(dups); #choose one with longer gene description
  dim(dups)
  dups[dups$mgi_symbol=='Ndor1',]
  ok3=c()
  for(gene in unique(dups$mgi_symbol)){
    x=dups[dups$mgi_symbol==gene,]
    ok3=rbind(ok3,x[order(x$ensembl_gene_id)[1],])
  }
  dim(ok3) #28 genes
  ok=rbind(ok1,ok2,ok3)
  dim(ok) #17864 genes
  id.mapping=ok;
  sum(duplicated(id.mapping$mgi_symbol)); #0
  data.table::fwrite(id.mapping,'~/Documents/Data_mouse_aging_atlas/droplet_17864genes_id.mapping.txt')
}
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/droplet_17864genes_id.mapping.txt')
dim(id.mapping) #17864 genes

