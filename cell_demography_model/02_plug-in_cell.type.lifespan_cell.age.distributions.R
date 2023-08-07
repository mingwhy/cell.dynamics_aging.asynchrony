
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
dim(df.est) #22
min(df.est$lifespan) # there are some lifespan as 0.9 day
# per day survival p, then, mean life expectancy Mean=-1/log(p)
# thus, p = exp(-1/Mean)

plug.in.values=exp(-1/df.est$lifespan) #survival rate per day
summary(plug.in.values)
sort(plug.in.values) #min 0.32, max, 0.9999965
length(unique(plug.in.values)) #21

df.est$plug.in.values=plug.in.values
df.est=df.est[order(df.est$lifespan),]
dim(df.est) #21 unique Ron_tc

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


k=seq(30,690,30) #2m=mature,(3-2)m*30=30 day; (24-2)*30=660 days
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
df.est$turnover=1- df.est$plug.in.values
df.est$turnover=sprintf('%.7f',df.est$turnover)

df.est$survival.rate=df.est$plug.in.values
#df.est$survival.rate=sprintf('%.7f',df.est$plug.in.values)
df.est=df.est[order(df.est$plug.in.values),]

df.est$tc.survival.rate=paste(df.est$Ron_tc,df.est$survival.rate)
df.est$tc.turnover.rate=paste(df.est$Ron_tc,df.est$turnover)
#tc.orders=df.est$tc.survival.rate
tc.orders=df.est$tc.turnover.rate

colnames(para_df)
para_df2=merge(df.est[,c('Ron_tc',"cell type annotation in TMS",'survival.rate','tc.turnover.rate')],para_df,
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
    xlab('Organismal age since maturation (day)'),
    theme(axis.text.x = element_text(angle = 45,hjust=1,size=15),
          legend.text=element_text(size=10)),
    #geom_vline(xintercept =c(720,11520,15840))
  #(3-2)*30*24=720, (18-2)*30*24=11520,(24-2)*30*24=15840
    geom_vline(xintercept =c(30,660),lty='dashed')
)

p1<-ggplot(para_df2,aes(x=day,y=mean,group=factor(tc.turnover.rate),col=factor(tc.turnover.rate)))+
  geom_point()+geom_line()+ scale_y_log10()+my_theme
p1

p2<-ggplot(para_df2,aes(x=day,y=sd,group=factor(tc.turnover.rate),col=factor(tc.turnover.rate)))+
  geom_point()+geom_line()+ scale_y_log10()+my_theme
p2  

######################################################################
## get mean diff and sd diff
df.2age=para_df2[para_df2$day %in% c(30,660),]
dim(df.2age) #21x2age=42 rows

df3=df.2age %>% group_by(Ron_tc,survival.rate,tc.turnover.rate) %>% 
    summarise(mean.diff=mean[day==660]-mean[day==30],
              sd.diff=sd[day==660]-sd[day==30],
              mean.ratio=(mean[day==660]/mean[day==30]),    
              sd.ratio=(sd[day==660]/sd[day==30]))

head(df3[order(df3$mean.diff,decreasing = T),])
head(df3[order(df3$mean.ratio,decreasing = T),])
head(df3[order(df3$sd.diff,decreasing = T),])
head(df3[order(df3$sd.ratio,decreasing = T),])
df.2age[df.2age$Ron_tc=='Nervous system:Neurons',]
df.2age[df.2age$Ron_tc=='Brain:endothelial cell',]

summary(df3$mean.diff)
summary(df3$sd.diff)

options(digits = 7)
p4=ggplot(df3,aes(x=mean.diff,y=sd.diff,col=tc.turnover.rate))+
  geom_point(size=10)+ 
  scale_x_continuous(trans = scales::pseudo_log_trans(1e-5, 800),
                                          breaks=c(0,10^seq(-5,1,2),10,20,50,100,200,500,800))+
  scale_y_continuous(trans = scales::pseudo_log_trans(1, 300),
                     breaks=c(0,2,5,10,20,50,100,200,300))+
  scale_color_viridis(name="turnover rate per day\ncell type annotation from Sender and Milo (2021)",option='turbo',alpha=0.6,discrete = T)+
  theme_classic(base_size = 16)+
#scale_colour_discrete(name="survival rate per day")
  theme(legend.text=element_text(size=10),
        axis.text.x = element_text(angle=45,hjust=1))
p4+xlab('mean(24m)-mean(3m)')+ylab('sd(24m)-sd(3m)')

pdf('plug_in_Rcal_mean_var_N0F.pdf',useDingbats = T,height=6,width = 16)
grid.arrange(p1+ylab('Mean cell age (day)')+
               guides(color = guide_legend(override.aes = list(size = 3) ) ),
             #p2+guides(color = guide_legend(override.aes = list(size = 3) ) ),
             ncol=1)
dev.off()


