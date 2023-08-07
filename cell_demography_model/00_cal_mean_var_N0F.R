
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
  #scale_colour_discrete(name="survival rate")+
  xlab('Organismal age since maturation (day)')

para_df[para_df$var<=0,] #check for numerical value resolution in R, sometimes, when survival.prob too high (0.9999999), stackoverflow
para_df$sd=para_df$var^0.5
summary(para_df$sd)

p2<-ggplot(para_df,aes(x=k,y=sd,col=factor(turnover)))+geom_point()+
  geom_line()+scale_y_log10()+
  scale_x_log10(breaks=c(10,100,1000,10000),
                labels = c(10,100,1000,10000))+
  #scale_y_continuous(trans = scales::pseudo_log_trans(1e-06, 10),breaks=c(0, 1e-06,1e-05,1e-04,1e-03,0.01, 0.1, 1))+
  theme_bw(base_size = 16)+ 
  scale_color_viridis(name="turnover rate",option='turbo',discrete = T)+
  #scale_colour_discrete(name="survival rate")+
  xlab('Organismal age since maturation (day)')
  #theme(legend.position = 'none')
p2


unique(para_df$k)
para_df[abs(para_df$k-30)<0.5 & (para_df$p==0.1 | para_df$p>=0.99),]
para_df[abs(para_df$k-545)<1 & (para_df$p==0.1 | para_df$p>=0.99),]

#####################################################
## plot population density for different p
ps=c(0.1,0.5,0.7,0.9,0.95,0.99,0.995,0.999,0.9999,0.99999) #survival.prob
#log10p=log10(p)
#p=10^log10p
ks=10^seq(1,4,1) # a range of organismal ages, in days
#ks=c(ks,5e4)
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
p3<-ggplot(df.out2,aes(x=log1age,y=ncell,col=factor(turnover)))+
  facet_grid(organismal.age~turnover)+#scale_y_log10()+
  #scale_y_continuous(trans = scales::pseudo_log_trans(1e-10, 10),
  #                   breaks=c(0, 10^seq(-10,-1,1),  1))+
  #scale_y_continuous(trans = scales::pseudo_log_trans(1e-6, 10),breaks=c(0, 10^seq(-6,-1,1),  1),
  #                   labels=c(0, 10^seq(-6,-1,1),  1))+
  scale_y_continuous(trans = scales::pseudo_log_trans(1e-5, 10),breaks=c(0, 10^seq(-5,-1,1),  1),
                     labels=c(0, 10^seq(-5,-1,1),  1))+
  #scale_y_continuous(trans = scales::pseudo_log_trans(1e-6, 10),breaks=10^seq(-6,-1,1))+
  #geom_bar(stat='identity')+
  #scale_color_viridis(name="survival rate",option='turbo',discrete = T)+
  scale_color_viridis(name="turnover rate\nper day",option='turbo',discrete = T)+
  #scale_colour_discrete(name="survival rate")+
  geom_point(size=0.9)+ylab('% Cell')+
  xlab('log10(1+cell.age), cell.age range from 0 to t')+
  theme_classic(base_size = 12)+theme(axis.text.y = element_text(size=8))

jpeg(file="cal_cell.age.distribution_N0F.jpeg",
    width=1400, height=800,quality = 20000,res= 120)
print(p3)
dev.off()

tmp=df.out[df.out$p==0.99999,]
#####################################################
## mean.diff VS sd.diff
library(tidyverse);library(viridis)
df.2age=para_df[para_df$k %in% c(10,1000),]
dim(df.2age) 
colnames(df.2age)[1]='survival.rate'
colnames(df.2age)[2]='organismal.age'

#df2=df.2age %>% group_by(survival.rate,organismal.age) %>% 
#              summarise(mean=sum(ncell*age,na.rm=T),
#                    var=sum(ncell*age^2,na.rm=T)-(mean(ncell*age,na.rm=T))^2)
df2=df.2age
df2$sd=df2$var^0.5

## mean, sd, difference  
df3=df2 %>% group_by(survival.rate) %>% summarise(
                mean.diff=mean[organismal.age==1000]-mean[organismal.age==10],
                sd.diff=sd[organismal.age==1000]-sd[organismal.age==10])
df3
sprintf('%.2f',df3$sd.diff)
summary(df3$mean.diff)   
summary(df3$sd.diff)   

df3$turnover=1-df3$survival.rate                                                                           
df3$turnover=sprintf('%.5f',df3$turnover)

p4=ggplot(df3,aes(x=mean.diff,y=sd.diff,col=factor(turnover)))+
  geom_point(size=10)+
  scale_x_continuous(trans = scales::pseudo_log_trans(1e-11, 1000),
                     breaks=c(10^seq(-11,0,2),  10,100,1000))+
  scale_y_continuous(trans = scales::pseudo_log_trans(1, 500),
                     breaks=c(0,2,5,10,20,50,100,200,500))+
  #scale_size(name='',range = c(0.1,12))+
  #scale_color_viridis(option='viridis')+
  #scale_color_viridis(option='magma')+
  scale_color_viridis(name="turnover rate\nper day",option='turbo',alpha=0.6,discrete = T)+
  theme_classic(base_size = 16)+
  #scale_colour_discrete(name="survival rate per day")
  theme(legend.text=element_text(size=10))
p4+xlab('mean(1000)-mean(10)')+ylab('sd(1000)-sd(10)')


## mean, sd, ratio
summary(df2$sd) #no zero
summary(df2$mean) #no zero
df4=df2 %>% group_by(survival.rate) %>% summarise(
  mean.ratio=mean[organismal.age==1000]/mean[organismal.age==10],
  sd.ratio=sd[organismal.age==1000]-sd[organismal.age==10])
df4
summary(df4$mean.ratio)   
summary(df4$sd.ratio)   

df4$turnover=1-df4$survival.rate                                                                           
df4$turnover=sprintf('%.5f',df4$turnover)

p5=ggplot(df4,aes(x=mean.ratio,y=sd.ratio,col=factor(turnover)))+
  geom_point(size=10)+scale_x_log10()+
  scale_y_continuous(trans = scales::pseudo_log_trans(1, 400),
                     breaks=c(0,2,5,10,20,50,100,200,400))+
  #scale_size(name='',range = c(0.1,12))+
  #scale_color_viridis(option='viridis')+
  #scale_color_viridis(option='magma')+
  scale_color_viridis(name="turnover rate",option='turbo',alpha=0.6,discrete = T)+
  theme_classic(base_size = 20)+
  #scale_colour_discrete(name="survival rate per day")
  theme(legend.text=element_text(size=10))
p5+xlab('mean(1000)/mean(10)')+ylab('sd(1000)/sd(10)')


pdf('cal_mean.cell.age.pdf',useDingbats = T,height = 6)
print(p1+scale_color_viridis(name="turnover rate\nper day",option='turbo',alpha=0.6,discrete = T)+
        ylab('Mean cell age (day)'))
dev.off()

