
######################################################################
## read in mouse turnover rate data
cell.lifespan=readxl::read_xlsx('Dataset_S1.xlsx');
dim(cell.lifespan) 
colnames(cell.lifespan);
head(cell.lifespan)

cell.lifespan=cell.lifespan[cell.lifespan$`included in this anallysis (1) or not (0)`==1,] #filter based on Kim's annotation
dim(cell.lifespan) 
cell.lifespan$lifespan=cell.lifespan$`lifespan.used.in.this.study (day)`

cell.lifespan[grep('-',cell.lifespan$lifespan),]
#https://onlinelibrary.wiley.com/doi/10.1111/imr.12693
cell.lifespan[cell.lifespan$lifespan=='40-70' & cell.lifespan$`cell type annotation in Sender and Milo 2021`=='Blood:mature T cells',]$lifespan='70'; #T cell 
cell.lifespan[cell.lifespan$lifespan=='40-70' & cell.lifespan$`cell type annotation in Sender and Milo 2021`=='Blood:mature B cells',]$lifespan='40'; #B cell 
#https://jamanetwork.com/journals/jamainternalmedicine/article-abstract/565579
cell.lifespan[cell.lifespan$lifespan=='200-400',]$lifespan='400'; #hepatocyte
cell.lifespan$lifespan=as.numeric(cell.lifespan$lifespan)
dim(cell.lifespan)  

tmp=cell.lifespan[!duplicated(cell.lifespan$`cell type annotation in Sender and Milo 2021`),]
table(tmp$species.of.lifespan.used.in.this.study) #1 human and 20 rodents

tc.orders.in.mouse=cell.lifespan[order(cell.lifespan$lifespan),]$`cell type annotation in TMS`
tc.orders.in.human=cell.lifespan[order(cell.lifespan$lifespan),]$`cell type annotation in Sender and Milo 2021`
tc.orders.in.human=unique(tc.orders.in.human)
