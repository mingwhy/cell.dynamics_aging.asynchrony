
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)
library(lme4);
library(lmerTest); #contain p value
library(variancePartition)
library(viridis)
# for grid.table (https://stackoverflow.com/questions/31776557/how-to-adjust-the-font-size-of-tablegrob)
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.6)),
  colhead = list(fg_params=list(cex = 0.6)),
  rowhead = list(fg_params=list(cex = 0.6)))

source('src_cellular_lifespan.R')
DM.out=readRDS('MAST_nonDEgene_stats.rds')
tc.names=names(DM.out)

used.stat='cv2'; 
keep.cell.types=intersect(tc.names,cell.lifespan$`cell type annotation in TMS`)

##########################################################
## gene-level metric: beta from LMM
if(TRUE){ #run once
  betas<-lapply(keep.cell.types,function(tc){
    all.ages=DM.out[[tc]]
    all.ages=all.ages[!is.na(all.ages$mean_scran),]
    
    df=all.ages[all.ages$age %in% c('3m','24m'),]
    if(length(unique(df$age))==1){return(NULL)}
    
    df=df[df$gene %in% names(which(table(df$gene)==2)),] #expr in both ages
    ngene=length(unique(df$gene))
    df$age=factor(df$age,levels=c('3m','24m'))
    
    df$numeric.age=as.numeric(gsub('m','',df$age))
    if(used.stat=='cv2'){
      lmer_fit1<-lmer(log10(cv2_scran)~numeric.age+(1|gene),data=df) #only random intercept
    }else{
      lmer_fit1<-lmer(log10(var_scran)~numeric.age+(1|gene),data=df) #only random intercept
    }
    
    var.out=calcVarPart(lmer_fit1)
    s1=paste(names(var.out),collapse = ' | ')
    s2=paste(round(var.out,4),collapse = ' | ')
    
    var.out=calcVarPart(lmer_fit1)
    PVE=var.out[2]
    
    x=summary(lmer_fit1)
    # get pvalue
    #https://www.r-bloggers.com/2014/02/three-ways-to-get-parameter-specific-p-values-from-lmer/
    coefs<-data.frame(coef(summary(lmer_fit1)))
    m.semTest<-lmer(log10(cv2_scran)~numeric.age+(1|gene),data=df,REML = FALSE) #only random intercept
    # get Satterthwaite-approximated degrees of freedom
    coefs$df.Satt <- coef(summary(m.semTest))[, 3]
    # get approximate p-values
    coefs$p.Satt <- coef(summary(m.semTest))[, 5]
    coefs
    
    x$coefficients
    b0=x$coefficients[1,1]
    b1=x$coefficients[2,1] #extract regression.coefficient of age
    std=x$coefficients[2,2]
    
    x=var.test(x=log10(df[df$age=='3m',]$var_scran),y=log10(df[df$age=='24m',]$var_scran),
               alternative='greater')#only apply to two groups
    
    x$estimate
    c(tc,PVE,b0,b1,std,x$estimate,coefs$p.Satt[2])
  })
  betas=Filter(Negate(is.null), betas)
  
  beta.out=as.data.frame(Reduce(`rbind`,betas))
  colnames(beta.out)=c('cell.type','PVE','intercept','beta','std','var.ratio','p.Staa')
  df.var.out=data.frame(tc=beta.out$cell.type,intercept=as.numeric(beta.out$intercept),
                        beta=as.numeric(beta.out$beta),
                        beta.std=as.numeric(beta.out$std),
                        PVE=as.numeric(beta.out$PVE),
                        var.ratio=as.numeric(beta.out$var.ratio),
                        pvalue=as.numeric(beta.out$p.Staa))
  dim(df.var.out) #38tc
  summary(df.var.out$pvalue); 
  sum(df.var.out$beta>0) #34 tc, increase with age
  
  head(df.var.out)
  if(used.stat=='cv2'){
    saveRDS(df.var.out,'LMM_cv2_out.rds')
  }
}

#######################################################################
## plot age var changing trend with age
if.plot=TRUE

res.LMM.gene=readRDS('LMM_cv2_out.rds')
res.LMM.gene$FDR=p.adjust(res.LMM.gene$pvalue,method='BH')

if(if.plot){
  plots1=lapply(keep.cell.types,function(tc){
    all.ages=DM.out[[tc]]
    all.ages=all.ages[!is.na(all.ages$mean_scran),]
    
    df=all.ages[all.ages$age %in% c('3m','24m'),]
    if(length(unique(df$age))==1){return(NULL)}
    
    df=df[df$gene %in% names(which(table(df$gene)==2)),] #have info in both ages
    
    ngene=length(unique(df$gene))
    df$age=factor(df$age,levels=c('3m','24m'))
    
    #ratio=res.log2ratio.gene[res.log2ratio.gene$cell.type==tc,]$ratio
    
    beta=res.LMM.gene[res.LMM.gene$tc==tc,]$beta
    FDR=res.LMM.gene[res.LMM.gene$tc==tc,]$FDR
    if(abs(FDR-0)<2.2e-16){FDR.sentence='FDR-adjust P value<2.2e-16'
    }else{FDR.sentence=paste0('FDR-adjust P value=',round(FDR,6))}  
    
    g1 <- ggplot(df,aes(x=age,y=log10(cv2_scran),group=gene,col=age) )+
      theme_bw()+geom_jitter(size=0.005)+scale_color_manual(values=c("#999999", "#E69F00"))+
      geom_line(lwd=0.03,col=grDevices::adjustcolor('grey10',alpha.f = 0.2))+
      ggtitle(paste0(tc,'\n#gene=',ngene,'\nbeta_age=',round(beta,5)))+
      ylab("cv2 (log10)")+
      theme(plot.title = element_text(size = 16),axis.title = element_text(size=14),
            legend.position = 'none')
    
    g1
  })
  plots1[[2]]
  plots1=Filter(Negate(is.null), plots1)
  length(plots1) #38

  jpeg(file="LMM_geneCV2.jpeg",width=2600, height=3600,quality = 20000,res= 120)
  grid.arrange(grobs=plots1,ncol=4)
  dev.off()
}



