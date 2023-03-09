# atlas info: https://github.com/czbiohub/tabula-muris-senis

# data source: https://figshare.com/projects/Tabula_Muris_Senis/64982
# TMS_processedFiles: Processed files (to use with scanpy)
# TMS.gene.data_final: TMS gene data (final) (code for paper: 2021-Mouse aging cell atlas analysis reveals global and cell type-specific aging signatures)
# TMS.Data.Objects: Tabula Muris Senis Data Objects


## read in h5ad file
#https://githubhot.com/repo/theislab/zellkonverter/issues/28
library(zellkonverter)
library(SummarizedExperiment)

#https://figshare.com/articles/dataset/tms_gene_data_rv1/12827615
sce<-readH5AD('TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');        
df.expr=assay(sce)
genes=rownames(df.expr)

## map gene id use biomart (https://www.biostars.org/p/301116/)
if(T){
  library(biomaRt)
  mouse= useEnsembl(version = 99, #this archived still contain dn ds values
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
  saveRDS(mouse_genes,'mmusculus_gene_ensemblv99.rds')
}
mouse_genes=readRDS('mmusculus_gene_ensemblv99.rds')
sum('0610009B22Rik' %in% mouse_genes$mgi_symbol )

## write id mapping table
if(T){
  length(genes) #fac: 22966, droplet: 22899
  sum(genes %in% mouse_genes$mgi_symbol) #fac:20449 , droplet: 20404
  
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

