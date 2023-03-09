
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(viridis)
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


p=c(0.1,0.5,0.7,0.9,0.95,0.99,0.995,0.999,0.9999,0.99999) #survival prob
k=10^seq(1,4,len=20) # a range of organismal ages, in days
k=sort(unique(c(k,10^seq(1,4,1))))
#k=c(k,5e4)

para_df=expand.grid(p,k)
colnames(para_df)=c('p','k')
para_df$mean<-apply(para_df,1,function(x) f1_age_mean(x[1],x[2]))
para_df$var<-apply(para_df,1,function(x) f1_age_var(x[1],x[2]))
summary(para_df$mean)
summary(para_df$var)

para_df$turnover=1-para_df$p
para_df$turnover=sprintf('%.5f',para_df$turnover)

p1<-ggplot(para_df,aes(x=k,y=mean,col=factor(turnover)))+geom_point()+
  geom_line()+ scale_y_log10()+
  scale_x_log10(breaks=c(10,100,1000,10000),
                labels = c(10,100,1000,10000))+
  theme_bw(base_size = 16)+ 
  scale_color_viridis(name="turnover rate",option='turbo',discrete = T)+
  xlab('Organismal age since maturation (day)')


pdf('cal_mean.cell.age.pdf',useDingbats = T,height = 6)
print(p1+scale_color_viridis(name="turnover rate\nper day",option='turbo',alpha=0.6,discrete = T)+
        ylab('Mean cell age (day)'))
dev.off()


#####################################################
## plot population density for different p
ps=c(0.1,0.5,0.7,0.9,0.95,0.99,0.995,0.999,0.9999,0.99999) #survival.prob
ks=10^seq(1,4,1) # a range of organismal ages, in days
pop_density<-function(t,p){
  c((1-p)*p^seq(0,t-1,1),p^t)
}
df=expand.grid(ps,ks)
colnames(df)=c('p','k') #survival.rate, organismal.age
out=lapply(1:nrow(df),function(i){
  x=as.numeric(df[i,])
  y=pop_density(t=x[2],p=x[1]) #possible age: 1,...,k
  #less than 1e-10 set as 0
  tmp=data.frame(p=x[1],k=x[2],age=seq(0,x[2],1),ncell=y);
  tmp
})

df.out=as.data.frame(Reduce(`rbind`,out))
summary(df.out$ncell)
min(df.out[df.out$ncell!=0,]$ncell)
df.out$organismal.age=factor(df.out$k)
df.out$log1age=log10(df.out$age+1)

df.out[df.out$ncell==0,]$ncell=NA
df.out$turnover=1-df.out$p
df.out$turnover=sprintf('%.5f',df.out$turnover)

df.out2=df.out;
df.out2$organismal.age=paste0(df.out2$organismal.age,' days')
p2<-ggplot(df.out2,aes(x=log1age,y=ncell,col=factor(turnover)))+
  facet_grid(organismal.age~turnover)+#scale_y_log10()+  
  scale_y_continuous(trans = scales::pseudo_log_trans(1e-5, 10),breaks=c(0, 10^seq(-5,-1,1),  1),
                     labels=c(0, 10^seq(-5,-1,1),  1))+
  scale_color_viridis(name="turnover rate\nper day",option='turbo',discrete = T)+
  geom_point(size=0.9)+ylab('% Cell')+
  xlab('log10(1+cell.age), cell.age range from 0 to t')+
  theme_classic(base_size = 12)+theme(axis.text.y = element_text(size=8))

jpeg(file="cal_cell.age.distribution_N0F.jpeg",
    width=1400, height=800,quality = 20000,res= 120)
print(p2)
dev.off()

