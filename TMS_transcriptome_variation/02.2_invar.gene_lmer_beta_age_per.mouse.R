
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)
library(lme4);
library(lmerTest); #contain p value
library(variancePartition)
library(viridis)

cell.lifespan=data.table::fread('../src/Dataset_S1.txt')
DM.out=readRDS('MAST_nonDEgene_stats.rds')
tc.names=names(DM.out)
length(tc.names);
length(DM.out)

used.stat='cv2'; 
keep.cell.types=intersect(tc.names,cell.lifespan$`cell type annotation in TMS`)
age.groups=c('3m','24m')

##########################################################
## gene-level metric: beta from LMM
if(TRUE){ #run once
  betas<-lapply(keep.cell.types,function(tc){
    cat(tc,'\n')
    all.ages=DM.out[[tc]]
    all.ages=all.ages[!is.na(all.ages$mean_scran),]
    
    df=all.ages[all.ages$age %in% age.groups,]
    if(length(unique(df$age))==1){return(NULL)}
    
    ngene=length(unique(df$gene))
    df$age=factor(df$age,levels=age.groups)
    
    df$numeric.age=as.numeric(gsub('m','',df$age))
    lmer_fit1<-lmer(log10(cv2_scran)~numeric.age+(1|mouse.id)+(1|gene),data=df)

    r<-tryCatch(
      {var.out=calcVarPart(lmer_fit1)},   
      error=function(c)'error')
    if(r=='error'){PVE=NA
    }else{
      var.out=calcVarPart(lmer_fit1)
      PVE=var.out[3] #gene, mouse.id, age, resids
    }
    
    x=summary(lmer_fit1)
    # get pvalue
    #https://www.r-bloggers.com/2014/02/three-ways-to-get-parameter-specific-p-values-from-lmer/
    coefs<-data.frame(coef(summary(lmer_fit1)))
    
    m.semTest<-lmer(log10(cv2_scran)~numeric.age+(1|mouse.id)+(1|gene),data=df,REML = FALSE) #only random intercept
    # get Satterthwaite-approximated degrees of freedom
    coefs$df.Satt <- coef(summary(m.semTest))[, 3]
    # get approximate p-values
    coefs$p.Satt <- coef(summary(m.semTest))[, 5]
    coefs
    
    x$coefficients
    b0=x$coefficients[1,1]
    b1=x$coefficients[2,1] #extract regression.coefficient of age
    std=x$coefficients[2,2]
    
    x=var.test(x=log10(df[df$age==age.groups[1],]$var_scran),y=log10(df[df$age==age.groups[2],]$var_scran),
               alternative='greater')#only apply to two groups
    x$estimate
    
    c(tc,PVE,b0,b1,std,x$estimate,coefs$p.Satt[2])
  })
  keep.cell.types[which(sapply(betas,is.null))] # "Heart:atrial myocyte", 11
  
  betas=Filter(Negate(is.null), betas)
  sapply(betas,length)
  
  beta.out=as.data.frame(Reduce(`rbind`,betas))
  colnames(beta.out)=c('cell.type','PVE','intercept','beta','std','var.ratio','p.Staa')
  df.var.out=data.frame(tc=beta.out$cell.type,intercept=as.numeric(beta.out$intercept),
                        beta=as.numeric(beta.out$beta),
                        beta.std=as.numeric(beta.out$std),
                        PVE=as.numeric(beta.out$PVE),
                        var.ratio=as.numeric(beta.out$var.ratio),
                        pvalue=as.numeric(beta.out$p.Staa))
  dim(df.var.out) #
  summary(df.var.out$pvalue); 
  df.var.out[which(is.na(df.var.out$pvalue)),]
  sum(df.var.out$beta>0) #17 tc, increase with age
  
  head(df.var.out)
  saveRDS(df.var.out,'LMM_cv2_out.rds')
}


#######################################################################
## plot age var changing trend with age
if.plot=FALSE

if(if.plot){
  res.LMM.gene=readRDS('LMM_cv2_out.rds')
  res.LMM.gene$FDR=p.adjust(res.LMM.gene$pvalue,method='BH')
  
  summary(res.LMM.gene$FDR)
  sum(res.LMM.gene$beta>0) #36
  sum(res.LMM.gene$FDR<0.05) #23
  sum(res.LMM.gene$beta>0 & res.LMM.gene$FDR<0.05) #22
  
  plots1=lapply(keep.cell.types,function(tc){
    all.ages=DM.out[[tc]]
    all.ages=all.ages[!is.na(all.ages$mean_scran),]
    
    df=all.ages[all.ages$age %in% age.groups,]
    if(length(unique(df$age))==1){return(NULL)}
    
    
    ngene=length(unique(df$gene))
    df$age=factor(df$age,levels=age.groups)
    
    beta=res.LMM.gene[res.LMM.gene$tc==tc,]$beta
    FDR=res.LMM.gene[res.LMM.gene$tc==tc,]$FDR
    if(abs(FDR-0)<0.05){FDR.sentence='FDR < 0.05'
    }else{FDR.sentence=paste0('FDR = ',format.pval(FDR,digits = 3))}  
    
    g1 <- ggplot(df,aes(x=age,y=log10(cv2_scran),group=gene,col=age) )+
      theme_bw()+geom_jitter(size=0.005)+#scale_color_manual(values=c("#999999", "#E69F00"))+
      geom_line(lwd=0.03,col=grDevices::adjustcolor('grey10',alpha.f = 0.2))+     
      ggtitle(paste0(tc,'\nbeta_age=',round(beta,5),'\n',FDR.sentence))+
      ylab("cv2 (log10)")+scale_color_manual(values=c("#999999", "#E69F00"))+
      theme(plot.title = element_text(size = 16),axis.title = element_text(size=14),
            legend.position = 'none')
    
    g1
  })
  plots1[[2]]
  plots1=Filter(Negate(is.null), plots1)
  length(plots1) 

  jpeg(file="LMM_geneCV2.jpeg",width=2600, height=3600,quality = 20000,res= 120)
  grid.arrange(grobs=plots1,ncol=4)
  dev.off()
}



