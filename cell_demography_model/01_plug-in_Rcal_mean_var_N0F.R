
options(digits=16)

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

######################################################################
## define time unit 
cell.lifespan$human_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`
df.est=cell.lifespan[!duplicated(cell.lifespan$human_tc),]
dim(df.est) #21
min(df.est$lifespan) # there are some lifespan as 0.9 day
# per day survival p, then, mean life expectancy Mean=-1/log(p)
# thus, p = exp(-1/Mean)

plug.in.values=exp(-1/df.est$lifespan) #survival rate per day
summary(plug.in.values)
unique(sort(plug.in.values)) #21 values

df.est$plug.in.values=plug.in.values
df.est=df.est[order(df.est$lifespan),]

######################################################################
library(ggplot2)
library(gridExtra)

## C = F
f1_age_mean<-function(p,k){
  k*p^k - (k*p^(k+1) - p^(k+1) -k*p^k+p)/(p-1)
}

f1_age_var<-function(p,k){
  part1= -2*k^2*p^(k+1) + k^2*p^(k+2) + k^2*p^k + 2*k*p^(k+1) + p^(k+1)-2*k*p^(k+2)+p^(k+2)-p^2-p
  part2 = k*p^(k+1)-p^(k+1)-k*p^k+p
  k^2*p^k - part1/(p-1)^2 - (k*p^k - part2/(p-1))^2
}

f1_age_var<-function(p,k){
  k^2*p^k - 
    (-2*k^2*p^(k+1) + k^2*p^(k+2) + k^2*p^k + 2*k*p^(k+1) + p^(k+1)-2*k*p^(k+2)+p^(k+2)-p^2-p)/(p-1)^2 -
    (k*p^k - (k*p^(k+1)-p^(k+1)-k*p^k+p)/(p-1))^2 
}

######################################################################
## look at 3m and 24m alone
if(!file.exists('plug-in_Rcal_mean_var_N0F.rds')){
  
  k=c(30,660);  # day unit: (3-2)*30=30, (24-2)*30=660 
  #k=c(60,1320);  # half-day unit: (3-2)*30*2=60, (24-2)*30*2=1320 
  #k=c(4,88); #3m=12week, 24m=24*4=96week, maturation=2m 12-8=4week, 96-8=88week.
  tmp=df.est[,c('human_tc','plug.in.values')]
  para_df=rbind(tmp,tmp)
  para_df$k=c(rep(k[1],nrow(df.est)),rep(k[2],nrow(df.est)))
  para_df$month=c(rep(3,nrow(df.est)),rep(24,nrow(df.est)))
  #para_df$hours= para_df$k * 12; #time.unit=half-day, 12hr
  para_df$week = para_df$k
  
  summary(para_df$plug.in.values)
  para_df$mean=0;
  for(i in 1:nrow(para_df)){
    para_df[i,'mean']=as.numeric(f1_age_mean(para_df[i,2],para_df[i,3]))
    cat('row',i,'is done',as.numeric(para_df[i,'mean']),'\n')
  }  
  para_df$var=0
  for(i in 1:nrow(para_df)){
    para_df[i,'var']=as.numeric(f1_age_var(para_df[i,2],para_df[i,3]))
    cat('row',i,'is done',as.numeric(para_df[i,'var']),'\n')
  }  
  summary(para_df$mean)
  summary(para_df$var)
  
  saveRDS(para_df,'plug-in_Rcal_mean_var_N0F.rds')
}





