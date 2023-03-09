
library(biomaRt)
mouse= useEnsembl(version = 99, #this archived still contain dn ds values
                             biomart = 'ENSEMBL_MART_ENSEMBL', 
                             dataset = 'mmusculus_gene_ensembl')

listAttributes(mouse) 
searchAttributes(mouse,'rnorvegicus_homolog')
searchAttributes(mouse, 'rnorvegicus_homolog_dn') #yes, there is xxx_dn
searchAttributes(mouse, 'rnorvegicus_homolog_ds') #yes, there is xxx_dn

mouse_rat <- getBM(attributes = c('ensembl_gene_id', 
                                     'rnorvegicus_homolog_ensembl_gene',  
                                     'rnorvegicus_homolog_dn', 
                                     'rnorvegicus_homolog_ds',
                                     'rnorvegicus_homolog_orthology_type',
                                  'rnorvegicus_homolog_orthology_confidence'), 
                      mart = mouse)
head(mouse_rat)
dim(mouse_rat) #70177     5

filter_mouse_rat=mouse_rat[!is.na(mouse_rat$rnorvegicus_homolog_dn),]
dim(filter_mouse_rat) #30033     5
head(filter_mouse_rat)

data.table::fwrite(filter_mouse_rat,'mouse_rat.dnds.txt')

