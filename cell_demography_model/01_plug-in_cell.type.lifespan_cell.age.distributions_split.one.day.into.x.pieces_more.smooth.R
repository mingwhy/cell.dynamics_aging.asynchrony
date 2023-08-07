
library(gridExtra)
library(grid)
library(tidyverse);library(viridis)
options(digits=16)

######################################################################
## read in mouse turnover rate data
cell.lifespan=data.table::fread('../src/Dataset_S1.txt')
cell.lifespan$Ron_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`

tmp=cell.lifespan[!duplicated(cell.lifespan$Ron_tc),]
dim(tmp); #unique 22 Ron_tc
table(tmp$species) #1 human and 20 rodents

tc.orders=cell.lifespan[order(cell.lifespan$lifespan),]$`cell type annotation in TMS`
tc.orders #39 tc

cell.lifespan[cell.lifespan$Ron_tc %in% names(which(table(cell.lifespan$Ron_tc)>1)),]
cell.lifespan[grep('endothelial',cell.lifespan$Ron_tc),] #4 Ron_tc, unique estimates
cell.lifespan[grep('endothelial',cell.lifespan$`cell type annotation in TMS`),]

mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.5)),
  colhead = list(fg_params=list(cex = 0.5)),
  rowhead = list(fg_params=list(cex = 0.5)))

#####################################################################
## define time unit 
df.est=cell.lifespan[!duplicated(cell.lifespan$Ron_tc),]
dim(df.est) #21
min(df.est$lifespan) # there are some lifespan as 0.9 day
# per day survival p, then, mean life expectancy Mean=-1/log(p)
# thus, p = exp(-1/Mean)

## add data points to make the curve smooth
sort(df.est$lifespan)
pseudo.lifespan=c(seq(10,100,5),seq(101,20000,len=500),max(df.est$lifespan))
pseudo.lifespan=c(sort(df.est$lifespan[df.est$lifespan<10]),pseudo.lifespan)
pseudo.lifespan=pseudo.lifespan[!duplicated(pseudo.lifespan)]
df.est=data.frame('Ron_tc'=paste('pseudo.tc',1:length(pseudo.lifespan)), lifespan=pseudo.lifespan)

split.one.day.into.x.pieces=20;

plug.in.values_per.day=exp(-1/(df.est$lifespan)) #finite survival rate 
plug.in.values=exp(-1/(df.est$lifespan*split.one.day.into.x.pieces)) #use half-half-day as unit, finite survival rate per half-day
log(plug.in.values_per.day)/log(plug.in.values) #split.one.day.into.x.pieces
summary(plug.in.values)
sort(plug.in.values) #min 0.32, max, 0.9999965
length(unique(plug.in.values)) #21

df.est$plug.in.values=plug.in.values
df.est$plug.in.values_per.day=plug.in.values_per.day
df.est=df.est[order(df.est$lifespan),]
dim(df.est) #200 pseudo tc

# turnover rate b= mortality.rate = 1-survival.rate= 1- exp(-1/T)
summary(df.est$lifespan) #ranges from 0.9 to 292000
simu.lifespan=c(0.1,0.5,0.9,seq(1,1000,100))
plot(simu.lifespan,1-exp(-1/simu.lifespan))

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
#p=c(1e-05,1e-03,1e-02,1e-01,0.3,0.5,0.7,0.8,0.9,0.95,0.99,0.9999)
p=plug.in.values;

k=split.one.day.into.x.pieces*seq(30,690,30) #2m=mature,(3-2)m*30=30 day; (24-2)*30=660 days
#k=c(k,1200,1400,1600,1800,2000)
#k = k -2*30*2; #use 2m as the onset of maturation in mosue
k=ceiling(k[k>0])
#k=seq(1,96,2) #3m=12week, 24m=24*4=96week, maturation=2m 12-8=4week, 96-8=88week.

para_df=expand.grid(p,k)
colnames(para_df)=c('p','k')
#para_df$hours= para_df$k * 12; #time.unit=half-day, 12hr
para_df$day= para_df$k

para_df$mean<-apply(para_df,1,function(x) f1_age_mean(x[1],x[2]))
para_df$var<-apply(para_df,1,function(x) f1_age_var(x[1],x[2]))
summary(para_df$mean)
summary(para_df$var)
para_df[is.infinite(para_df$var),]

sort(unique(para_df$p),decreasing = T) #largest survival rate
para_df[para_df$var<=0,] #check for numerical value resolution in R, sometimes, when survival.prob too high (0.9999999), stackoverflow
#para_df[para_df$var<=0,]$var=0 

para_df$sd=sqrt(para_df$var)
summary(para_df$var)
head(para_df[order(para_df$var,decreasing = T),])
#para_df[para_df$var<=0,]$var=0; #for 0.999999
# when larger than 0.9, there is numerical value problem in R

# get tc orders in plotting
df.est$turnover=1- df.est$plug.in.values #plug in are survival rate values
df.est$turnover=sprintf('%.7f',df.est$turnover)

df.est$survival.rate_per.day=df.est$plug.in.values_per.day
df.est$survival.rate=df.est$plug.in.values
#df.est$survival.rate=sprintf('%.7f',df.est$plug.in.values)
df.est=df.est[order(df.est$plug.in.values),]

df.est$tc.survival.rate=paste(df.est$Ron_tc,df.est$survival.rate)
df.est$tc.turnover.rate=paste(df.est$Ron_tc,df.est$turnover)
#tc.orders=df.est$tc.survival.rate
tc.orders=df.est$tc.turnover.rate

colnames(para_df)
para_df2=merge(df.est[,c('Ron_tc',
                         "survival.rate_per.day" ,'survival.rate','tc.turnover.rate')],para_df,
               by.x='survival.rate',by.y='p')
para_df2$tc.turnover.rate=factor(para_df2$tc.turnover.rate,levels=tc.orders)


summary(para_df$day)
my_theme=list(
  #scale_x_log10(breaks=c(10,100,1000,2000,10000,20000,30000),labels = c(10,100,1000,2000,10000,20000,30000))+
  scale_x_log10(breaks=c(30,50,100,200,400,600,700),labels = c(30,50,100,200,400,600,700)),
  #scale_x_log10(breaks=c(10,100,200,500,1000,1500,2000,5000,10000,15000,20000,25000),
  #              labels = c(10,100,200,500,1000,1500,2000,5000,10000,15000,20000,25000)),
  theme_classic(base_size = 20),
  scale_color_viridis(name="turnover rate per day",option='turbo',alpha=0.6,discrete = T),
  #scale_colour_discrete(name="survival rate per day"),
  xlab('organismal age since maturation (days)'),
  theme(axis.text.x = element_text(angle = 45,hjust=1,size=15),
        legend.text=element_text(size=10)),
  #geom_vline(xintercept =c(720,11520,15840))
  #(3-2)*30*24=720, (18-2)*30*24=11520,(24-2)*30*24=15840
  geom_vline(xintercept =split.one.day.into.x.pieces*c(30,660),lty='dashed')
)

#p1<-ggplot(para_df2,aes(x=day,y=mean,group=factor(tc.turnover.rate),col=factor(tc.turnover.rate)))+
#  geom_point()+geom_line()+ scale_y_log10()+my_theme+theme(legend.position = 'none')
#p1

#p2<-ggplot(para_df2,aes(x=day,y=sd,group=factor(tc.turnover.rate),col=factor(tc.turnover.rate)))+
#  geom_point()+geom_line()+ scale_y_log10()+my_theme+theme(legend.position = 'none')
#p2  

######################################################################
## get mean diff and sd diff
split.one.day.into.x.pieces
df.2age=para_df2[para_df2$day %in% (split.one.day.into.x.pieces*c(30,660)),]
dim(df.2age) #21x2age=42 rows
colnames(df.2age)
#df.2age$day=df.2age$day/split.one.day.into.x.pieces;
#table(df.2age$day) #30day and 660 day

summary(df.2age$mean) # time unit = 1day / split.one.day.into.x.pieces
df.2age$mean.in.days=df.2age$mean/split.one.day.into.x.pieces
#infinite case, mean=1/u
df.2age$mortality=(1-df.2age$survival.rate)
df.2age$mortality_per.day=(1-df.2age$survival.rate_per.day)
df.2age$expected= -1/log(1-df.2age$mortality) #or turnover rate
df.2age$expected.in.days=df.2age$expected/split.one.day.into.x.pieces

1/(1-df.2age$survival.rate_per.day); #long lifespan are fine
df.2age$expected.in.days #some discrepancy for small lifespan 
log(df.2age$survival.rate_per.day)/log(df.2age$survival.rate) #split.one.day.into.x.pieces
log(1-df.2age$mortality_per.day)/log(1-df.2age$mortality)

ggplot(df.2age,aes(x=mortality,y=mean.in.days,group=day))+
  facet_wrap(.~day)+geom_point()+theme_classic()+
  scale_x_log10()+scale_y_log10()+
  geom_point(aes(x=survival.rate,y=expected),col='red')

summary(df.2age$mortality)
df.2age$log10mortality=log10(df.2age$mortality)
#df.2age$mortality_per.day=log10(1-df.2age$survival.rate_per.day)
x_breaks = 10^seq(-6,1,1)/ split.one.day.into.x.pieces
x_labels = x_breaks *split.one.day.into.x.pieces
p=ggplot(df.2age,aes(x=log10mortality,y=mean.in.days,group=day,col=factor(day)))+
  #p=ggplot(df.2age,aes(x=mortality_per.day,y=mean.in.days,group=day,col=factor(day)))+
  #geom_point()+
  geom_line(lwd=2)+theme_classic(base_size = 20)+
  scale_x_continuous(breaks=log10(x_breaks),labels=x_labels)+
  scale_y_log10()+
  #geom_point(aes(x=log10mortality,y=expected.in.days),col='black')+
  geom_line(aes(x=log10mortality,y=expected.in.days),col='black',lwd=2)+
  #geom_point(aes(x=mortality_per.day,y=expected.in.days),col='black')+
  #geom_line(aes(x=mortality_per.day,y=expected.in.days),col='black')+
  xlab('Turnover rate per day')+ylab('Mean cell age (day)')
p
grid.arrange(p,
             p+scale_x_continuous(trans=scales::pseudo_log_trans(10^-6,2),
                                  breaks=c(0,log10(x_breaks),1),labels=c(0,x_breaks,1))+
               theme(axis.text.x = element_text(angle = 45,hjust=1,vjust = 1,size=10)),
             ncol=1)

pdf('mean.cell.age_truncated_expo_smooth.pdf',useDingbats = T,height = 6,width = 8)
print(p+scale_color_manual(name='Organismal age',values=c("#56B4E9","#E69F00"),
                           labels=c('3m','24m'))+
        theme(axis.text= element_text(hjust=0.5,vjust = 0,size=15)))
dev.off()



