
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)
library(MAST)
library(ggpointdensity)
library(viridis)

min.ncell_per.mouse=15;

cell.lifespan=data.table::fread('../src/Dataset_S1.txt')

sce.minCellCounts=readRDS('sce_minCellCounts.rds')
names(sce.minCellCounts)

###########################################################
## use MAST test for age effect for each gene
## `DGE_analysis.R` from https://github.com/czbiohub/tabula-muris-senis/tree/master/2_aging_signature/job.DGE_analysis
if(!file.exists('MAST_DEtest.rds')){
  start_time =Sys.time()
  
  de_res.out<-lapply(names(sce.minCellCounts),function(tc){
    sce_naive=sce.minCellCounts[[tc]]
    assayNames(sce_naive)='counts'
    
    # Prepare sca object
    sca <- SceToSingleCellAssay(sce_naive, class = "SingleCellAssay",check_sanity = FALSE)
    colData(sca)$age_num = as.numeric(gsub('m','',colData(sca)$age))
    colData(sca)$scaled_n_genes = scale(colData(sca)$n_genes) # n_gene (CDR)
    sca_filt = sca[rowSums(assay(sca)) != 0, ]
    
    # DGE testing 
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
  end_time =Sys.time()
  end_time-start_time; 
  
}

###########################################################
## generate gene stats and add mouse.id as one column
df.de_res.out=readRDS('MAST_DEtest.rds')
head(df.de_res.out)
tc.names=sort(unique(df.de_res.out$tc)) 
length(unique(tc.names))

cutoff=0.1;
tmp=df.de_res.out %>% group_by(tc) %>% summarise(n.total.gene=n(),n.invar.gene=sum(age.H_p>cutoff))
tmp$prop=tmp$n.invar.gene/tmp$n.total.gene
tmp=tmp[order(tmp$prop),]
summary(tmp$n.total.gene)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#366    2616    3738    3562    4292    5937 

summary(tmp$n.invar.gene)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#286    1280    2328    2051    2744    3943 

summary(tmp$prop)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1570  0.4822  0.6228  0.5685  0.7034  0.8257 

tmp2=merge(tmp,cell.lifespan[,c(1,2,7)],by.x='tc',by.y='cell type annotation in TMS')
cor.test.out=cor.test(tmp2$prop,log(tmp2$lifespan))
plot(tmp2$prop,log(tmp2$lifespan),xlab='%invar.gene',ylab='log(cell.lifespan)',pch=16,
     main=paste0('Pearson r=',round(cor.test.out$estimate,3),', P = ',round(cor.test.out$p.value,3)))

de_res=df.de_res.out[df.de_res.out$age.H_p>cutoff,] # mean-invar genes

## compute mean, var, CV2, DM per gene for mean-stable genes

DM.out<-lapply(tc.names,function(tc){
  #cat(tc,"\n")
  sce_naive=sce.minCellCounts[[tc]]
  assayNames(sce_naive)<-'counts'
  
  tmp=de_res[de_res$tc==tc,]
  sce_naive=sce_naive[rownames(sce_naive) %in% tmp$gene,]
  expr.m=assay(sce_naive,'counts')
  
  # each ind contain >=20 cells per age group per tc
  ids<-unique(sce_naive$mouse.id)
  
  df.one.tc<-lapply(ids,function(id){
    sce_naive2=sce_naive[,sce_naive$mouse.id==id]
    expr.m=assay(sce_naive2,'counts')
    mouse.age<-names(which(table(sce_naive2$age)>=min.ncell_per.mouse)) # "Heart:atrial myocyte", 13 and 17 cells
    if(length(mouse.age)==0){return(NULL)}
    
    df=lapply(mouse.age,function(i){
      m=expr.m[,sce_naive2$age==i,drop=F]
      parameter_df=data.frame(gene=rownames(sce_naive2))
      parameter_df$age=i
      parameter_df$mean_scran <- rowMeans(m)
      parameter_df$var_scran <- rowVars(as.matrix(m))
      parameter_df$cv2_scran <-  parameter_df$var_scran/parameter_df$mean_scran^2 
      parameter_df$DM <- scran::DM(
        mean = parameter_df$mean_scran,
        cv2 = parameter_df$cv2_scran
      )
      parameter_df$tc=tc;
      parameter_df=parameter_df[!is.nan(parameter_df$cv2_scran),] #when mean=0, cv2_scran=NaN
      return(parameter_df)
    })
    df=as.data.frame(Reduce(`rbind`,df))
    df$mouse.id=as.character(id)
    return(df)
  })
  sapply(df.one.tc,length)
  df.one.tc2<-Filter(Negate(is.null), df.one.tc)
  
  if(length(df.one.tc2)==0){return(NULL)
  }else{
    return(as.data.frame(Reduce(`rbind`,df.one.tc2)))
  }
})

length(DM.out)
sapply(DM.out,length)
DM.out2<-Filter(Negate(is.null), DM.out)
head(DM.out2[[1]])
sapply(DM.out2,dim)
length(unique(unlist(lapply(DM.out2,'[[','tc'))))

names(DM.out2)=unique(unlist(lapply(DM.out2,'[[','tc')))
table(DM.out2[[2]]$mouse.id)
saveRDS(DM.out2,'MAST_nonDEgene_stats.rds')




