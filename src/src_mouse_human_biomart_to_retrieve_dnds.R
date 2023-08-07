
library(biomaRt)
mouse= useEnsembl(version = 99, #this archived still contain dn ds values
                  biomart = 'ENSEMBL_MART_ENSEMBL', 
                  dataset = 'mmusculus_gene_ensembl')

listAttributes(mouse) 
grep('hsapiens',listAttributes(mouse)$name )
listAttributes(mouse)$name[grep('hsapiens',listAttributes(mouse)$name )]

searchAttributes(mouse,'hsapiens_homolog')
searchAttributes(mouse, 'hsapiens_homolog_dn') #yes, there is xxx_dn
searchAttributes(mouse, 'hsapiens_homolog_ds"') #yes, there is xxx_dn


## mouse and Rat (https://support.bioconductor.org/p/9141035/)
mouse_human <- getBM(attributes = c('ensembl_gene_id', 
                                  'hsapiens_homolog_ensembl_gene',  
                                  'hsapiens_homolog_dn', 
                                  'hsapiens_homolog_ds',
                                  "hsapiens_homolog_orthology_type" ,
                                  "hsapiens_homolog_orthology_confidence" ), 
                   mart = mouse)
head(mouse_human)
dim(mouse_human) #61688     5

filter_mouse_human=mouse_human[!is.na(mouse_human$hsapiens_homolog_dn),]
dim(filter_mouse_human) #23089     5
head(filter_mouse_human)

data.table::fwrite(filter_mouse_human,'mouse_human.dnds.txt')

