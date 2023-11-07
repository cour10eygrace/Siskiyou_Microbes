library(ggplot2)
library(dplyr)

#read in data 
om<-read.csv("DATA/Soil Enviro final.csv")
om$Site<-as.factor(om$Site)
str(om)

#pull out sites from microbes MS 
om<-subset(om, Site=="1"|Site=="9"|Site=="10"|Site=="3"|Site=="6"|Site=="7"|Site=="8"|Site=="11"|Site=="5")
om<-mutate(om, homeaway=ifelse(Site=="1"|Site=="9"|Site=="10", "home", "away"))%>%
  mutate(success=case_when(Site=="5"&homeaway=="away"~"AwaySuccessful",
                           Site!="5"&homeaway=="away"~"AwayUnsuccessful", 
                           TRUE~"Home"))
names(om)



#supplementary figures & stats
#OM 
hist(om$OM)
hist(om$Log_OM) #not really normally distributed...

ggplot(om,aes(x=Site, y=Log_OM, fill=success))+ 
  geom_boxplot()+ theme_bw() + 
  scale_fill_manual(values=c( "springgreen3", "firebrick1", "steelblue")) 

omplot<-ggplot(om,aes(x=homeaway, y=Log_OM, fill=success))+ 
  geom_boxplot() + theme_bw()+ theme(legend.position = 'none')+xlab("")+
  scale_fill_manual(values=c( "#00BA38", "#F8766D", "#619CFF")) 

legend<-ggplot(om,aes(x=homeaway, y=Log_OM, fill=success))+ 
  geom_boxplot() + theme_bw()+xlab("")+
  scale_fill_manual(values=c( "#00BA38", "#F8766D", "#619CFF")) 

summary(aov(Log_OM~success, om)) #NS
ommod<-aov(Log_OM~success, om) #NS
TukeyHSD(ommod)

#N
hist(om$N)
hist(om$Log_N)#normal ish

ggplot(om,aes(x=Site, y=Log_N, fill=success))+ 
  geom_boxplot()+ theme_bw() + 
  scale_fill_manual(values=c( "springgreen3", "firebrick1", "steelblue")) 

Nplot<-ggplot(om,aes(x=homeaway, y=Log_N, fill=success))+ 
  geom_boxplot() + theme_bw() + theme(legend.position = 'none') +xlab("")+
  scale_fill_manual(values=c( "#00BA38", "#F8766D", "#619CFF")) 

summary(aov(Log_N~success, om)) #NS
Nmod<-aov(Log_N~success, om) #NS
TukeyHSD(Nmod)

#pH
hist(om$PH) #not really normally distributed... range 6-7.4
hist(om$Log_PH)#worse

ggplot(om,aes(x=Site, y=PH, fill=success))+ 
  geom_boxplot()+ theme_bw() + 
  scale_fill_manual(values=c( "springgreen3", "firebrick1", "steelblue")) 

pHplot<-ggplot(om,aes(x=homeaway, y=PH, fill=success))+ 
  geom_boxplot() + theme_bw() + theme(legend.position = 'none') +xlab("")+
  scale_fill_manual(values=c( "#00BA38", "#F8766D", "#619CFF")) 


summary(aov(PH~success, om)) #lower pH (more acidic) in successful soils than unsuccessful and home 
phmod<-aov(PH~success, om)
TukeyHSD(phmod)


#phosphorus 
#weak bray extraction- The Bray-P1 test works well for most soils that are slightly alkaline to highly acidic (pH of 7.4 or less). 
hist(om$P_WB)
hist(om$Log_P_WB)#not really normally distributed... 

ggplot(om,aes(x=Site, y=Log_P_WB, fill=success))+ 
  geom_boxplot()+ theme_bw() + 
  scale_fill_manual(values=c( "springgreen3", "firebrick1", "steelblue")) 

Pplot<-ggplot(om,aes(x=homeaway, y=Log_P_WB, fill=success))+ 
  geom_boxplot()+  theme_bw() + theme(legend.position = 'none')+xlab("")+
  scale_fill_manual(values=c( "#00BA38", "#F8766D", "#619CFF")) 


summary(aov(Log_P_WB~success, om)) #NS
Pmod<-aov(Log_P_WB~success, om)
TukeyHSD(Pmod)

#K
hist(om$K)
hist(om$Log_K)#normal ish

ggplot(om,aes(x=Site, y=Log_K, fill=success))+ 
  geom_boxplot()+ theme_bw() + 
  scale_fill_manual(values=c( "springgreen3", "firebrick1", "steelblue")) 

Kplot<-ggplot(om,aes(x=homeaway, y=Log_K, fill=success))+ 
  geom_boxplot() +  theme_bw() + theme(legend.position = 'none')+xlab("")+
  scale_fill_manual(values=c( "#00BA38", "#F8766D", "#619CFF")) 

summary(aov(Log_K~success, om)) #NS
Kmod<-aov(Log_K~success, om)
TukeyHSD(Kmod)

#CA:MG
hist(om$CA.MG)
hist(om$Log_CA.MG)#normal ish

ggplot(om,aes(x=Site, y=Log_CA.MG, fill=success))+ 
  geom_boxplot()+ theme_bw() + 
  scale_fill_manual(values=c( "springgreen3", "firebrick1", "steelblue")) 

camgplot<-ggplot(om,aes(x=homeaway, y=Log_CA.MG, fill=success))+ 
  geom_boxplot() +  theme_bw() + theme(legend.position = 'none')+xlab("")+
  scale_fill_manual(values=c( "#00BA38", "#F8766D", "#619CFF")) 

summary(aov(Log_CA.MG~success, om))  
camgmod<-aov(Log_CA.MG~success, om)
TukeyHSD(camgmod) #successful and unsuccessful same 

#CEC
hist(om$CEC)
hist(om$Log_CEC)#normal ish

ggplot(om,aes(x=Site, y=Log_CEC, fill=success))+ 
  geom_boxplot()+ theme_bw() + 
  scale_fill_manual(values=c( "springgreen3", "firebrick1", "steelblue")) 

cecplot<-ggplot(om,aes(x=homeaway, y=Log_CEC, fill=success))+ 
  geom_boxplot() +  theme_bw() + theme(legend.position = 'none')+xlab("")+
  scale_fill_manual(values=c( "#00BA38", "#F8766D", "#619CFF")) 


summary(aov(Log_CEC~success, om)) #lower pH (more acidic) in successful soils than unsuccessful and home 
cecmod<-aov(Log_CEC~success, om)
TukeyHSD(cecmod) #lower CEC in successful (because of low pH)


#joint plot 
gridExtra::grid.arrange(Nplot, Pplot, Kplot, camgplot, pHplot, cecplot, omplot, legend, nrow=3)

#principal components
ggplot(om,aes(x=Site, y=Prin1, fill=success))+ 
  geom_boxplot()+ theme_bw() + 
  scale_fill_manual(values=c( "springgreen3", "firebrick1", "steelblue")) 

ggplot(om,aes(x=homeaway, y=Prin1, fill=success))+ 
  geom_boxplot()

prin1mod<-aov(Prin1~success, om)
summary(prin1mod)
TukeyHSD(prin1mod) #successful & unsuccessful NS

ggplot(om,aes(x=Site, y=Prin2, fill=success))+ 
  geom_boxplot()+ theme_bw() + 
  scale_fill_manual(values=c( "springgreen3", "firebrick1", "steelblue")) 

ggplot(om,aes(x=homeaway, y=Prin2, fill=success))+ 
  geom_boxplot()

prin2mod<-aov(Prin2~success, om)
summary(prin2mod)
TukeyHSD(prin2mod) #successful loads negatively on pc2 (higher CA:MG & lower CEC)

