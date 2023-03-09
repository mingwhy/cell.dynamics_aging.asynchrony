
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)
library(MAST)
library(ggpointdensity)
library(viridis)
library(patchwork); #for plot_annotation
plot_params <- list(
  geom_pointdensity(size = 0.6),
  scale_colour_viridis(name = "Density"),
  theme(  text = element_text(size = rel(3)),
          legend.position = "bottom",
          legend.text = element_text(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
          legend.key.size = unit(0.018, "npc")
  ) )

source('src_cellular_lifespan.R')
sce.minCellCounts=readRDS('sce_minCellCounts.rds')
tc.names=names(sce.minCellCounts) 
tc.names=intersect(tc.names,cell.lifespan$`cell type annotation in TMS`); #39 tc
sce.minCellCounts=sce.minCellCounts[tc.names]
 
###########################################################
## use MAST test for age effect for each gene
## `DGE_analysis.R` from https://github.com/czbiohub/tabula-muris-senis/tree/master/2_aging_signature/job.DGE_analysis
if(!file.exists('MAST_DEtest.rds')){
  start_time =proc.time()
  de_res.out<-lapply(tc.names,function(tc){
    sce_naive=sce.minCellCounts[[tc]]
    assayNames(sce_naive)
    
    # Prepare sca object
    sca <- SceToSingleCellAssay(sce_naive, class = "SingleCellAssay",check_sanity = FALSE)
    colData(sca)$age_num = as.numeric(gsub('m','',colData(sca)$age))
    colData(sca)$scaled_n_genes = scale(colData(sca)$n_genes) # n_gene (CDR)
    sca_filt = sca[rowSums(assay(sca)) != 0, ]
    
    # DGE testing (only male, so no `sex` as covariate)
    covariate = ''
    covariate = paste0(covariate, " + scaled_n_genes");
    print(paste0('covariate: ', covariate))  
    
    zlmCond <- zlm(formula = as.formula(paste0("~age_num", covariate)), sca=sca_filt)
    summaryCond <- summary(zlmCond, doLRT="age_num")
    
    # Summarize results 
    summaryDt <- summaryCond$datatable
    dt1 = summaryDt[contrast=="age_num" & component=="H", .(primerid, `Pr(>Chisq)`)]
    dt2 = summaryDt[contrast=="age_num" & component=="logFC", .(primerid, coef, z)]
    de_res = merge(dt1, dt2, by="primerid")
    colnames(de_res) <- c("gene", "age.H_p", "age.logFC", 'age.logFC_z')
    de_res$age.H_fdr <- p.adjust(de_res$age.H_p, "fdr")
    
    de_res$tc=tc
    return(de_res)
  })
  
  df.de_res.out=as.data.frame(Reduce(`rbind`,de_res.out))
  saveRDS(df.de_res.out,'MAST_DEtest.rds')
  
  print('Finished')
  print(proc.time() - start_time) #2~3hr
}

###########################################################
df.de_res.out=readRDS('MAST_DEtest.rds')
head(df.de_res.out)
tc.names=sort(unique(df.de_res.out$tc)) #39 tc

tmp=df.de_res.out %>% group_by(tc) %>% summarise(sum(age.H_p>0.1)/n())
tmp=tmp[order(tmp$`sum(age.H_p > 0.1)/n()`),]
summary(tmp$`sum(age.H_p > 0.1)/n()`)
summary(tmp[tmp$tc %in% cell.lifespan$`tissue: cell.type in mouse`,]$`sum(age.H_p > 0.1)/n()`)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1676  0.5233  0.6406  0.6010  0.7130  0.8329 
data.table::fwrite(tmp,'nonDEgene.prop.txt')

de_res=df.de_res.out[df.de_res.out$age.H_p>0.1,] # mean-invar genes

## compute mean, var, CV2, DM per gene for mean-stable genes
DM.out<-lapply(tc.names,function(tc){
  sce_naive=sce.minCellCounts[[tc]]
  assayNames(sce_naive)
  #table(Matrix::colSums(assay(sce_naive,'counts'))) #all the same library size (already down-sampled)
  
  tmp=de_res[de_res$tc==tc,]
  #summary(tmp$age.H_fdr)
  #summary(tmp$age.logFC)
  sce_naive=sce_naive[rownames(sce_naive) %in% tmp$gene,]
  expr.m=assay(sce_naive,'counts')
  
  df=lapply(levels(sce_naive$age),function(i){
    m=expr.m[,sce_naive$age==i]
    parameter_df=data.frame(gene=rownames(sce_naive))
    parameter_df$age=i
    parameter_df$mean_scran <- rowMeans(m)
    parameter_df$var_scran <- rowVars(as.matrix(m))
    parameter_df$cv2_scran <-  parameter_df$var_scran/parameter_df$mean_scran^2 
    parameter_df$DM <- DM(
      mean = parameter_df$mean_scran,
      cv2 = parameter_df$cv2_scran
    )
    parameter_df=parameter_df[!is.nan(parameter_df$cv2_scran),] #when mean=0, cv2_scran=NaN
    #summary(parameter_df$cv2_scran)
    return(parameter_df)
  })
  df=as.data.frame(Reduce(`rbind`,df))
  return(df)
})

length(DM.out)
names(DM.out)=tc.names
sapply(DM.out,dim)
saveRDS(DM.out,'MAST_nonDEgene_stats.rds')

