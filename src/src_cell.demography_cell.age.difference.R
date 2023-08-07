
options(digits=16)

######################################################################
## read in mouse turnover rate data
cell.lifespan=data.table::fread('Dataset_S1.txt')
dim(cell.lifespan) 


######################################################################
## define time unit 
cell.lifespan$Ron_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`
df.est=cell.lifespan[!duplicated(cell.lifespan$Ron_tc),]
dim(df.est) #22
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
## look at pairwise mean cell age difference: 3m, 18m, 24m
# month to days:  (3-2)*30=30, (18-2)*30=480, (24-2)*30=660

age.groups=c(3,18,24)

for(a1 in age.groups){
  for(a2 in age.groups){
   
    age1=(a1-2)*30;
    age2=(a2-2)*30
    output.file=paste0("change.in.cell.age_",a1,'_vs_',a2,'.rds');
    cat(output.file,'\n')
    if(a2<=a1) next
    #k=c(30,480);  # day unit: (3-2)*30=30, (18-2)*30=480 
    #k=c(30,660);  # day unit: (3-2)*30=30, (24-2)*30=660 
    k=c(age1,age2)
    tmp=df.est[,c('Ron_tc','plug.in.values')]
    para_df=rbind(tmp,tmp)
    para_df$k=c(rep(k[1],nrow(df.est)),rep(k[2],nrow(df.est)))
    
    para_df$month=c(rep(a1,nrow(df.est)),rep(a2,nrow(df.est)))
    #para_df$month=c(rep(3,nrow(df.est)),rep(24,nrow(df.est)))
    #para_df$hours= para_df$k * 12; #time.unit=half-day, 12hr
    para_df$day = para_df$k
    
    summary(para_df$plug.in.values)
    para_df$mean=0;
    for(i in 1:nrow(para_df)){
      para_df[i,'mean']=as.numeric(f1_age_mean(para_df[i,2],para_df[i,3]))
      #cat('row',i,'is done',as.numeric(para_df[i,'mean']),'\n')
    }  
    para_df$var=0
    for(i in 1:nrow(para_df)){
      para_df[i,'var']=as.numeric(f1_age_var(para_df[i,2],para_df[i,3]))
      #cat('row',i,'is done',as.numeric(para_df[i,'var']),'\n')
    }  
    summary(para_df$mean)
    summary(para_df$var)
    
    saveRDS(para_df,output.file)
  }
}





