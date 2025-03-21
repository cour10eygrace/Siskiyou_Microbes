#load packages----
library(ggplot2)
library(phyloseq)
library(microbiome)
library(RColorBrewer)
library(dplyr)
library(knitr)
library(tidyselect)
library(gridExtra)
library(ggpubr) 
library(extrafont) 
library(lme4)
library(lmerTest)
library(DescTools)
library(knitr)
library(tidyverse)
library(rngtools)


#load stored phyloseq objects----
#Bacteria
load("DATA/physeq.bac.Rdata")

#Fungi
load("DATA/physeq.fungi.Rdata")
gc()

#First let's consider differences between sites where Horkelia 
#seeds were collected (Home) and the one successful transplant 
#site (Away)

#comparing site 5 (away) with sites 1,9,and 10 (seed home) for 
#both fungi and Bacteria. Basically want to ask if microbes where 
#seeds were collected (sites 1,9,10) differ from site 5 where the
#new population is. 

#subset physeq & metadata objects 
physeq.bac.sub<-subset_samples(physeq.bac, Site=="1"| Site=="9"|Site=="10"| Site=="5")
physeq.bac.away<-subset_samples(physeq.fungi,  Site=="5")
  
physeq.fungi.sub<-subset_samples(physeq.fungi, Site=="1"| Site=="9"|Site=="10"| Site=="5")
physeq.fungi.away<-subset_samples(physeq.fungi,  Site=="5")

metabsub<-subset(metab, Site=="1"| Site=="9"|Site=="10"| Site=="5")
metafsub<-subset(metaf, Site=="1"| Site=="9"|Site=="10"| Site=="5")
metabaway<-subset(metab,  Site=="5")
metafaway<-subset(metaf,  Site=="5")

metafsub<-subset(metafsub, Code!="W05X05XA" & Code!="W10X10XA")#remove 2 samples after rarefying fungi 
metafaway<-subset(metafaway, Code!="W05X05XA")#remove 1 sample after rarefying fungi  


##Alpha diversity---- 
#calculate
alphadiv_bac<-phyloseq::estimate_richness(physeq.bac)
alphadiv_bac<-merge(alphadiv_bac, metab, by=0)

alphadiv_bac_sub<-phyloseq::estimate_richness(physeq.bac.sub)
alphadiv_bac_sub<-merge(alphadiv_bac_sub, metabsub, by=0)

alphadiv_fun<-phyloseq::estimate_richness(physeq.fungi)
alphadiv_fun<-merge(alphadiv_fun, metaf, by=0)

alphadiv_fun_sub<-phyloseq::estimate_richness(physeq.fungi.sub)
alphadiv_fun_sub<-merge(alphadiv_fun_sub, metabsub, by=0)

#plot----
#away no success (site 3,6,7,8,11) vs away success (site 5)
c<-ggplot(subset(alphadiv_fun,Loc=="Away", Code!="W05X05XA"), aes(x = Hsuccess, y=Shannon, fill=Spp)) + geom_boxplot()+theme_classic()+ xlab(" ")+
  scale_x_discrete(labels=c("Away unsuccessful", "Away successful"))+
  scale_fill_manual(values=c('darkgrey', 'black'), name="Source", 
                     labels = c("Horkelia","Bare ground"))+ ylab("Shannon diversity ")

#home (site 1,9,10) vs away success (site 5)
#reorder so home gets plotted first 
alphadiv_fun_sub<-group_by(alphadiv_fun_sub, Loc)%>%mutate(order=min(Site))
alphadiv_fun_sub$Loc<-with(alphadiv_fun_sub, reorder(Loc, order))

a<-ggplot(alphadiv_fun_sub,  aes(x = Loc, y=Shannon, fill=Spp)) + geom_boxplot()+theme_classic()+ xlab("Site")+            
  scale_x_discrete(labels=c("Home", "Away successful"))+
  scale_fill_manual(values=c('darkgrey', 'black'), name="Source", 
                     labels = c("Horkelia","Bare ground"))+ ylab("Shannon diversity")
          
d<-ggplot(subset(alphadiv_bac,Loc=="Away"), aes(x = Hsuccess, y=Shannon, fill=Spp)) + geom_boxplot()+theme_classic()+ xlab(" ")+
  scale_x_discrete(labels=c("Away unsuccessful", "Away successful"))+
  scale_fill_manual(values=c('darkgrey', 'black'), name="Source", 
                    labels = c("Horkelia","Bare ground"))+ ylab("Shannon diversity ")

#reorder so home gets plotted first 
alphadiv_bac_sub<-group_by(alphadiv_bac_sub, Loc)%>%mutate(order=min(Site))
alphadiv_bac_sub$Loc<-with(alphadiv_bac_sub, reorder(Loc, order))

b<-ggplot(alphadiv_bac_sub,  aes(x = Loc, y=Shannon, fill=Spp)) + geom_boxplot()+theme_classic()+ xlab("Site")+            
  scale_x_discrete(labels=c("Home", "Away successful"))+
  scale_fill_manual(values=c('darkgrey', 'black'), name="Source", 
                    labels = c("Horkelia","Bare ground"))+ ylab("Shannon diversity")

#plot 2 panel Fig alpha div 
#pdf("Fig2.pdf", width = 12, height = 8)
#gridExtra::grid.arrange(a, b, c,d)  
#dev.off()

#Anovas----
aov_alphadiv_bac<-aov(log(Shannon) ~ Hsuccess*Spp, 
                      subset(alphadiv_bac, Loc=="Away"))
summary(aov_alphadiv_bac)
TukeyHSD(aov_alphadiv_bac)
   
aov_alphadiv_bac2<-aov(log(Shannon) ~ Hsuccess*Spp, alphadiv_bac_sub)
summary(aov_alphadiv_bac2)
TukeyHSD(aov_alphadiv_bac2)

aov_alphadiv_fun<-aov(log(Shannon) ~ Hsuccess*Spp, 
subset(alphadiv_fun, Loc=="Away", 
       Code!="W05X05XA"))#remove samples after rarefying fungi ))
summary(aov_alphadiv_fun)
TukeyHSD(aov_alphadiv_fun)

aov_alphadiv_fun2<-aov(log(Shannon) ~ Hsuccess*Spp, alphadiv_fun_sub)
summary(aov_alphadiv_fun2)
TukeyHSD(aov_alphadiv_fun2)

#Beta diversity- permanovas----
#Bacteria
bac_bray_sub<-distance(physeq.bac.sub, method = "bray")#slow
bac_bray_sub<-as.matrix(bac_bray_sub)
vegan::adonis2(bac_bray_sub~Loc*Spp, data=metabsub, permutations = 999)# Loc sig 

#Fungi 
fungi_bray_sub<-distance(physeq.fungi.sub, method = "bray")
fungi_bray_sub<-as.matrix(fungi_bray_sub)
vegan::adonis2(fungi_bray_sub~Loc*Spp, data=metafsub, permutations = 999)#all sig

#now test away soils only for H vs bare 
fungi_bray_away<-distance(physeq.fungi.away, method = "bray")
fungi_bray_away<-as.matrix(fungi_bray_away)

vegan::adonis2(fungi_bray_away~Spp, data=metafaway, permutations = 999)#not different 

##Plot ordinations----
# pull out legend function
get_only_legend <- function(plot) { 
    
# get tabular interpretation of plot 
plot_table <- ggplot_gtable(ggplot_build(plot))  
    
#  Mark only legend in plot 
legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")  
                              
# extract legend 
legend <- plot_table$grobs[[legend_plot]] 
                              
# return legend 
return(legend)  
}

#Bacteria
bac.ord.sub<-ordinate(physeq.bac.sub, "PCoA", "bray")
fig3b<-plot_ordination(physeq.bac.sub,bac.ord.sub, color = "Loc",
                shape="Spp")+theme_classic()+geom_point(size=2)+
  scale_color_manual(values=c('black', "darkgrey"), name="Site", 
                    labels = c("Away successful", "Home"))+
  scale_shape_manual(values=c('square', 'circle'), name="Source", 
                labels = c("Horkelia","Bare ground")) 

legend<-get_only_legend(fig3b) #pull out legend

fig3b<-fig3b+theme( legend.position = "none" ) #rewrite w/ no legend  
#Fungi 
fung.ord.sub<-ordinate(physeq.fungi.sub, "PCoA", "bray")

fig3a<-plot_ordination(physeq.fungi.sub,fung.ord.sub, color = "Loc", 
                shape="Spp")+theme_classic()+geom_point(size=2)+
  scale_color_manual(values=c('black', 'darkgrey'), name="Site")+  
                      #labels = c("Away successful", "Home"))+
  scale_shape_manual(values=c('square', 'circle'), name="Source", 
                     labels = c("Horkelia","Bare ground"))+
                   theme( legend.position = "none" )

#So microbial communities differ between home and away 
#for Horkelia soils for Bacteria and fungi and between Horkelia
#and bare for fungi but not between Horkleia and bare at away sites. 


#Now let's consider differences between away sites where Horkelia 
#was and was not successful

#For site 8 there was only one sample, so it was dropped because we can’t compare it to anything statistically.
physeq.bac.away<-subset_samples(physeq.bac, Loc=="Away")#& Site!="8") 
physeq.fungi.away<-subset_samples(physeq.fungi, Loc=="Away")#& Site!="8") 

metabaway<-subset(metab, Loc=="Away")#& Site!="8") 
metafaway<-subset(metaf, Loc=="Away") #& Site!="8") 
metafaway<-subset(metafaway, Code!="W06N05HA"&
  Code!="W07N04HC"& Code!="W07N06XA"& Code!="W05X05XA") #remove missing samples 

##Run permanovas
#Bacteria
bac_bray_away2<-distance(physeq.bac.away, method = "bray")#slow
bac_bray_away2<-as.matrix(bac_bray_away2)
vegan::adonis2(bac_bray_away2~Hsuccess*Spp, data=metabaway, permutations = 999)

#Fungi 
fungi_bray_away2<-distance(physeq.fungi.away, method = "bray")
fungi_bray_away2<-as.matrix(fungi_bray_away2)
vegan::adonis2(fungi_bray_away2~Hsuccess*Spp, data=metafaway, permutations = 999)


##Plot ordinations

#Bacteria
bac.ord.away<-ordinate(physeq.bac.away, "PCoA", "bray")#slow

fig3d<-plot_ordination(physeq.bac.away,bac.ord.away, color = "Hsuccess", 
                       shape="Spp")+theme_classic()+geom_point(size=2)+
  scale_color_manual(values=c('darkgrey', 'black'), name="Site", 
                     labels = c("Away unsuccessful","Away successful"))+
  scale_shape_manual(values=c('square', 'circle'), name="Source", 
                     labels = c("Horkelia","Bare ground"))+
  guides(color = guide_legend(order = 2),
         shape = guide_legend(order = 1))


legend2<-get_only_legend(fig3d) #pull out legend
fig3d<-fig3d+ theme( legend.position = "none" ) #rewrite w/ no legend

#Fungi 
fung.ord.away<-ordinate(physeq.fungi.away, "PCoA", "bray")
fig3c<-plot_ordination(physeq.fungi.away,fung.ord.away, color="Hsuccess",
                       shape="Spp")+theme_classic()+geom_point(size=2)+
  scale_color_manual(values=c('darkgrey', 'black'), name="Site")+
  scale_shape_manual(values=c('square', 'circle'), name="Source", 
                     labels = c("Horkelia","Bare ground"))+
  theme( legend.position = "none" ) #rewrite w/ no legend


#plot 4 panel fig 3
#pdf("Fig3.pdf", width = 12, height = 10)
#grid.arrange(fig3a, fig3b, legend, fig3c, fig3d, legend2, ncol=3)
#dev.off()


#So microbial communities differ between succuessful and 
#unsuccessful away sites as well 

#Funtional guilds----
#What about them differs? Look at Funguild data
guilds<- read.csv("DATA/SerpentineITS.guilds.csv")
guilds<-tidyr::separate(guilds, col=Taxonomy,into=c("x","x1", "x2", "x3","x4", "x5", "k", "k1", "p", "p1","c", "c1","o", "o1","f", "f1","g", "g1","s", "s1", "s2", "s3"))

#filter out anything not fungi #3562
ITS_Fungi<-filter(guilds,k1=="Fungi")%>%dplyr::select(-x,-x1,-x2, -x3,-x4,-x5, -k, -p, -c, -o,-f,-g,-s)
#calculate only funguild assigned OTUs #2321
ITS_Funguild<-subset(ITS_Fungi, Guild!="-")

#Calculate pathogens
path<-filter(ITS_Fungi, Guild=="Plant Pathogen")
path2<-path[,c(1:149)]
path2<-path2[, -1]
path_by_sample<-path2 %>%tidyr::gather(key = 'sample', value='value')%>%group_by(sample)%>%dplyr::mutate(rich=ifelse(value>0, 1, 0))%>%
  dplyr::summarise(pathabund=sum(value), pathrich=sum(rich))

#calculate mutualists 
mut<-filter(ITS_Fungi, Guild=="Arbuscular Mycorrhizal")
mut2<-mut[,c(1:149)]
mut2<-mut2[, -1]
mut_by_sample<-mut2 %>%tidyr::gather(key = 'sample', value='value')%>%group_by(sample)%>%mutate(rich=ifelse(value>0, 1, 0))%>%
  summarise(mutabund=sum(value), mutrich=sum(rich))

funguild_abund_rich<-left_join(mut_by_sample, path_by_sample)
funguild_abund_rich$mpratio<-funguild_abund_rich$mutabund/funguild_abund_rich$pathabund
funguild_abund_rich$mprich<-funguild_abund_rich$mutrich/funguild_abund_rich$pathrich
funguild_abund_rich<-rename(funguild_abund_rich,Code=sample)%>% left_join(metaf, .)%>%filter(!is.na(mpratio))

funguild_abund_richNS<-subset(funguild_abund_rich, Site!="5") #everything except successful site 

#calculate # of taxa classified to each----
muttax<-subset(guilds, Guild=="Arbuscular Mycorrhizal")
muttax<-muttax[, c(156:174)]
muttax<-distinct(muttax)

unique(muttax$c1)
unique(muttax$o1)
unique(muttax$f1)
unique(muttax$g1)
unique(muttax$s1)
unique(muttax$s2)

pathtax<-filter(guilds, Guild=="Plant Pathogen")
pathtax<-pathtax[, c(156:174)]
pathtax<-distinct(pathtax)

unique(pathtax$c1)
unique(pathtax$o1)
unique(pathtax$f1)
unique(pathtax$g1)
unique(pathtax$s1)
unique(pathtax$s2)

funguild_taxtable<-rbind(muttax, pathtax)
#write.csv(funguild_taxtable, "tableS2.csv")

#plot different comparisons----
#reorder factors so home gets plotted first 
funguild_abund_rich<-group_by(funguild_abund_rich, Loc)%>%mutate(order=min(Site))
funguild_abund_rich$Loc<-with(funguild_abund_rich, reorder(Loc, order))

fig4c<-ggplot(data=subset(funguild_abund_rich, Site=="1"| Site=="9"|Site=="10"| Site=="5"),
              aes(x = Loc, y=log(mpratio), fill=Spp)) + geom_boxplot()+theme_classic()+ xlab("Site")+            scale_x_discrete(labels=c("Home", "Away successful"))+  scale_fill_manual(values=c('darkgrey', 'black'), name="Source", 
                    labels = c("Horkelia","Bare ground"))+ ylab("mutualist:pathogen ratio")

fig4d<-ggplot(data=subset(funguild_abund_rich, Site=="1"| Site=="9"|Site=="10"| Site=="5"),
              aes(x = Loc, y=log(mprich), fill=Spp)) + geom_boxplot()+theme_classic()+ xlab("Site")+     scale_x_discrete(labels=c( "Home", "Away successful"))+
  scale_fill_manual(values=c('darkgrey', 'black'), name="Source", 
                    labels = c("Horkelia","Bare ground"))+ ylab("mutualist:pathogen richness")

#save dataset for anovas later
funguild_abund_rich0<- subset(funguild_abund_rich, Site=="1"| Site=="9"|Site=="10"| Site=="5")

#update to retain only away sites 
funguild_abund_rich<-subset(funguild_abund_rich,Loc=="Away")

fig4a<-ggplot(data=funguild_abund_rich, aes(x = Hsuccess, y=log(mpratio), fill=Spp)) + geom_boxplot()+theme_classic()+ xlab(" ")+
  scale_x_discrete(labels=c("Away unsuccessful", "Away successful"))+ scale_fill_manual(values=c('darkgrey', 'black'), name="Source", 
                                                                                        labels = c("Horkelia","Bare ground"))+ ylab("mutualist:pathogen ratio")
fig4b<-ggplot(data=funguild_abund_rich, aes(x = Hsuccess, y=log(mprich), fill=Spp)) + geom_boxplot()+theme_classic()+ xlab(" ")+
  scale_x_discrete(labels=c("Away unsuccessful", "Away successful"))+ scale_fill_manual(values=c('darkgrey', 'black'), name="Source", 
                                                                                        labels = c("Horkelia","Bare ground"))+ ylab("mutualist:pathogen richness")

#plot 2 panel fig4
#pdf("Fig4.pdf", width = 12, height = 8)
#ggpubr::ggarrange(fig4a, fig4b, fig4c, fig4d, common.legend = T, legend="right")
#dev.off()


#do home and away soils differ?
ggplot(data=funguild_abund_rich, aes(x = Hsuccess, y=log(mprich))) + facet_wrap(~Loc)+ geom_boxplot()+theme_classic()+xlab("Site")+ scale_fill_manual(values=c('darkgrey', 'black'), name="Source")+ ylab("mutualist:pathogen")
#make new factor 
funguild_abund_rich<-mutate(funguild_abund_rich, homeaway=ifelse(Site=="1"|Site=="9"|Site=="10", "home", "away"))%>%
  mutate(success=case_when(Site=="5"&homeaway=="away"~"AwaySuccessful",
                           Site!="5"&homeaway=="away"~"AwayUnsuccessful", 
                           TRUE~"Home"))

ggplot(data=funguild_abund_rich, aes(x = success, y=log(mprich))) + geom_boxplot()+theme_classic()+xlab("Site")+ scale_fill_manual(values=c('darkgrey', 'black'), name="Source")+ ylab("mutualist:pathogen") #Make Fig S2

aov_mpratio_home<-aov(log(mpratio+1)~ success, funguild_abund_rich)
summary(aov_mpratio_home)
TukeyHSD(aov_mpratio_home)

aov_mprich_home<-aov(log(mprich+1)~ success, funguild_abund_rich)
summary(aov_mprich_home)
TukeyHSD(aov_mprich_home)

#away unsuccessful vs home - not different 
ggplot(data=funguild_abund_richNS, aes(x = Loc, y=log(mpratio))) + geom_boxplot()+theme_classic()+ xlab("Site")+     scale_x_discrete(labels=c("Away unsuccessful", "Home"))+
  scale_fill_manual(values=c('darkgrey', 'black'), name="Source", 
                    labels = c("Horkelia","Bare ground"))+ ylab("mutualist:pathogen")
ggplot(data=funguild_abund_richNS, aes(x = Loc, y=log(mpratio), fill=Spp)) + geom_boxplot()+theme_classic()+ xlab("Site")+     scale_x_discrete(labels=c("Away unsuccessful", "Home"))+
  scale_fill_manual(values=c('darkgrey', 'black'), name="Source", 
                    labels = c("Horkelia","Bare ground"))+ ylab("mutualist:pathogen")

#So mutualist to pathogen ratios (abundance and richness) are 
#higher at the successful Horkelia sites, but they don't differ
#between Horkelia and bare ground at the successful site, 
#which suggests that the mutualists are not in *response* to 
#Horkelia presence

#let's test whether these are singificant 

#ANOVAS----
library(rcompanion)
#successful vs unsuccessful away (fig 3a, b)
funguild_abund_rich$mprichT<-transformTukey(funguild_abund_rich$mprich)#lambda=0.35
funguild_abund_rich$mpratioT<-transformTukey(funguild_abund_rich$mpratio)#lambda=0.2

hist(log(funguild_abund_rich$mprich+1))
aov_mprich<-aov(log(mprich+1)~ Hsuccess*Spp, funguild_abund_rich)
summary(aov_mprich)
TukeyHSD(aov_mprich)

#transformTukey(funguild_abund_rich$mpratio)#lambda 0.125
hist(log(funguild_abund_rich$mpratio+1))
aov_mpratio<-aov(log(mpratio+1)~ Hsuccess*Spp, funguild_abund_rich)
summary(aov_mpratio)
TukeyHSD(aov_mpratio)

aov_mpratio<-aov(mpratioT~ Hsuccess*Spp, funguild_abund_rich)
summary(aov_mpratio)
TukeyHSD(aov_mpratio)

aov_mprich<-(aov(mprichT~ Hsuccess*Spp, funguild_abund_rich))
summary(aov_mprich)
TukeyHSD(aov_mprich)

#home vs successful (fig3 c,d)
funguild_abund_rich0$mprichT<-transformTukey(funguild_abund_rich0$mprich)#lambda=0.25
funguild_abund_rich0$mpratioT<-transformTukey(funguild_abund_rich0$mpratio)#lambda=0.125

aov_mpratio0<-aov(mpratioT~ Loc*Spp, funguild_abund_rich0)
summary(aov_mpratio0)
TukeyHSD(aov_mpratio0)

aov_mprich0<-(aov(mprichT~ Loc*Spp, funguild_abund_rich0))
summary(aov_mprich0)
TukeyHSD(aov_mprich0)

#home vs unsuccessful (no fig)
funguild_abund_richNS$mprichT<-transformTukey(funguild_abund_richNS$mprich)#lambda=0.325
funguild_abund_richNS$mpratioT<-transformTukey(funguild_abund_richNS$mpratio)#lambda=0.15

aov_mpratioNS<-aov(mpratioT~ Loc*Spp, funguild_abund_richNS)
summary(aov_mpratioNS)
TukeyHSD(aov_mpratioNS)

aov_mprichNS<-(aov(mprichT~ Loc*Spp, funguild_abund_richNS))
summary(aov_mprichNS)
TukeyHSD(aov_mprichNS)

#So M:P ratios were higher at the successful site and also higher
#in bare soils overall. 
#M:P ratios are not different between Horkelia and bare ground at 
#the sucessful site but at the unsuccessful site, Horkelia has 
#lower M:P than bare soil. Also M:P in unsuccessful bare soil at 
#was not different than in successful Horkelia soil, 
#but is slightly lower than successful bare ground.  

#Do taxa differ in these soils across a similar pattern?
  
#Taxa plots---- 
  
#MUST RUN R VERSION 4.0.3 for this!!! 
#Bacteria----
getPalette = colorRampPalette(brewer.pal(5, "Set1"))
#group by phylum and remove NAs
physeq.bac.nona.phylum = subset_taxa(physeq.bac, Phylum != "NA") 
speciesPalette = getPalette(9)

#split into 4 groups away Y success (H and bare), away N success (H and bare)
physeq.awayY.phylum= subset_samples(physeq.bac.nona.phylum, Site=="5")
physeq.awayN.phylum= subset_samples(physeq.bac.nona.phylum, Loc== "Away"& Site!="5")

#Microbiome plots 
ps.b.phy.awayY <- microbiome::aggregate_top_taxa(physeq.awayY.phylum,"Phylum", top =5)
ps.b.phy.awayY.rel <- microbiome::transform(ps.b.phy.awayY, "compositional")

ps.b.phy.awayN <- microbiome::aggregate_top_taxa(physeq.awayN.phylum,"Phylum", top =5)
ps.b.phy.awayN.rel <- microbiome::transform(ps.b.phy.awayN, "compositional")

plot.awayY <- microbiome::plot_composition(ps.b.phy.awayY.rel, 
                                           group_by = "Spp")+ 
  ggtitle("Away soils successful") + xlab(" ")+ ylab(" ")+
  scale_fill_manual(values= speciesPalette)+ 
  theme(axis.ticks.y =element_blank(),axis.text.x =element_blank(),
        axis.text.y =element_blank(), plot.title = element_text(size = 10))

plot.awayN <- microbiome::plot_composition(ps.b.phy.awayN.rel,
                                           group_by = "Spp")+ 
  ggtitle("Away soils unsuccessful") + xlab(" ")+ ylab(" ")+
  scale_fill_manual(values= speciesPalette)+ 
  theme( axis.ticks.y =element_blank(), axis.text.x =element_blank(), 
         axis.text.y =element_blank(), plot.title = element_text(size = 10))

grid.arrange(plot.awayN,plot.awayY)

#Fungi ---- 
#group by phylum and remove NAs
physeq.fungi.nona.phylum = subset_taxa(physeq.fungi, Phylum != "NA") 
#create color pallete 
#speciesList = unique(tax_table(physeq.fungi.prune.nona.phylum)[,"Phylum"])

#split into 4 groups away Y success (H and bare), away N success (H and bare)
physeq.awayY.phylum= subset_samples(physeq.fungi.nona.phylum, Site=="5")
physeq.awayN.phylum= subset_samples(physeq.fungi.nona.phylum, Loc== "Away"& Site!="5")

#Microbiome plots 
ps.f.phy.awayY <- microbiome::aggregate_top_taxa(physeq.awayY.phylum,"Phylum", top =5)
ps.f.phy.awayY.rel <- microbiome::transform(ps.f.phy.awayY, "compositional")

ps.f.phy.awayN <- microbiome::aggregate_top_taxa(physeq.awayN.phylum,"Phylum", top =5)
ps.f.phy.awayN.rel <- microbiome::transform(ps.f.phy.awayN, "compositional")

plot.awayY <- microbiome::plot_composition(ps.f.phy.awayY.rel, 
                                           group_by = "Spp")+ 
  ggtitle("Away soils successful") + xlab(" ")+ ylab(" ")+
  scale_fill_manual(values= speciesPalette)+ 
  theme(axis.ticks.y =element_blank(),
        axis.text.y =element_blank(),
        axis.text.x =element_blank(), plot.title = element_text(size = 10))

plot.awayN <- microbiome::plot_composition(ps.f.phy.awayN.rel,
                                           group_by = "Spp")+ 
  ggtitle("Away soils unsuccessful") + xlab(" ")+ ylab(" ")+
  scale_fill_manual(values= speciesPalette)+ 
  theme( axis.ticks.y =element_blank(), 
         axis.text.x =element_blank(),
         axis.text.y =element_blank(), plot.title = element_text(size = 10))

grid.arrange(plot.awayN,plot.awayY)


#Somewhat hard to tell from the plots but let's run some tests 
#for differential abundance of phlya 

#relative abundance tests----
#Bacteria
physeq.bac.nona.phylum= subset_samples(physeq.bac.nona.phylum, Site=="5"|
                                           Site=="1"|Site=="10"|Site=="9")

ps.bac.phy <- microbiome::aggregate_top_taxa(physeq.bac.nona.phylum, "Phylum", top = 5)
bphyreads<-as.data.frame(t(abundances(ps.bac.phy)))
bphyreads$Code<-row.names(bphyreads)
bphyreads<-left_join(bphyreads, metab)

shapiro.test(log(bphyreads$p__Acidobacteria))
acido<-(aov((log(p__Acidobacteria)~ Hsuccess*Spp), bphyreads))
summary(acido)
TukeyHSD(acido)#yes H, bare < no bare 

shapiro.test(log(bphyreads$p__Actinobacteria))
actino<-aov(log(p__Actinobacteria)~ Hsuccess*Spp, bphyreads)
summary(actino)
TukeyHSD(actino)#yes bare > no bare 

shapiro.test(log(bphyreads$p__Verrucomicrobia))
verruco<-aov(log(p__Verrucomicrobia)~ Hsuccess*Spp, bphyreads)
summary(verruco)
TukeyHSD(verruco)#yes H, bare> no H, bare


shapiro.test(log(bphyreads$p__Planctomycetes))
plancto<-aov(log(p__Planctomycetes)~ Hsuccess*Spp, bphyreads)
summary(plancto)
TukeyHSD(plancto)#yes H, bare> no H, bare

shapiro.test(log(bphyreads$p__Proteobacteria))
proteo<-aov(log(p__Proteobacteria)~ Hsuccess*Spp, bphyreads)
summary(proteo)
TukeyHSD(proteo)#yes H < no bare, #yes bare< no H 

#Fungi
physeq.fungi.nona.phylum= subset_samples(physeq.fungi.nona.phylum, Site=="5"|
                                           Site=="1"|Site=="10"|Site=="9")

ps.fung.phy <- microbiome::aggregate_top_taxa(physeq.fungi.nona.phylum, "Phylum", top = 5)
fphyreads<-as.data.frame(t(abundances(ps.fung.phy)))
fphyreads$Code<-row.names(fphyreads)
fphyreads<-left_join(fphyreads, metaf)

#fungi data are not normally distributed run KW tests 
kruskal.test(Ascomycota~ Hsuccess , fphyreads)#sig
NemenyiTest(fphyreads$Ascomycota, as.factor(fphyreads$Hsuccess), dist="tukey") #Y>N
kruskal.test(Ascomycota~ successx , fphyreads)#sig
NemenyiTest(fphyreads$Ascomycota, as.factor(fphyreads$successx), dist="tukey") #yes H, bare > no bare 

kruskal.test(Basidiomycota~ Hsuccess , fphyreads)#NS
kruskal.test(Basidiomycota~ successx , fphyreads)#NS
NemenyiTest(fphyreads$Basidiomycota, as.factor(fphyreads$successx), dist="tukey") #no bare > no H 

kruskal.test(Glomeromycota~ Hsuccess, fphyreads)#sig
NemenyiTest(fphyreads$Glomeromycota, as.factor(fphyreads$Hsuccess), dist="tukey") #Y>N
kruskal.test(Glomeromycota~ successx, fphyreads)#sig
NemenyiTest(fphyreads$Glomeromycota, as.factor(fphyreads$successx), dist="tukey") #Yes H and bare > no H 
 #yes bare > no bare  

kruskal.test(Mucoromycota~ Hsuccess , fphyreads)#sig
NemenyiTest(fphyreads$Mucoromycota, as.factor(fphyreads$Hsuccess), dist="tukey") #Y > N
kruskal.test(Mucoromycota~ successx , fphyreads)#sig
NemenyiTest(fphyreads$Mucoromycota, as.factor(fphyreads$successx), dist="tukey") #

kruskal.test(Mortierellomycota~ Hsuccess , fphyreads)#NS
NemenyiTest(fphyreads$Mortierellomycota, as.factor(fphyreads$Hsuccess), dist="tukey") kruskal.test(Mortierellomycota~ successx , fphyreads)#NS
NemenyiTest(fphyreads$Mortierellomycota, as.factor(fphyreads$successx), dist="tukey") 


#For Bacteria,  successful Horkelia away sites have lower 
#Acidobacteria, but higher Verrucomicrobia and Planctomyctes 
#than unsuccussful away sites. All of these phyla are associated 
#with slower C mineralization i.e. oligotrophic strategy (Ho et al 2017 FEMS)

#For Fungi, successful Horkelia away sites have lower 
#Mortierellomycota, but higher Ascomycota, Glomeromycota and 
#Mucoromycota than unsuccussful away sites. Glomeromycota (AMF)
#may be important for the establishment and growth of Horkelia 
#seedlings. This confirms the pattern we see in M:P ratios. 

#Let's zoom in on the Glomeromycota results 
kruskal.test(Glomeromycota~ Hsuccess, fphyreads)#sig
NemenyiTest(fphyreads$Glomeromycota, as.factor(fphyreads$Hsuccess), dist="tukey") #Y>N
kruskal.test(Glomeromycota~ successx, fphyreads)#sig
NemenyiTest(fphyreads$Glomeromycota, as.factor(fphyreads$successx), dist="tukey") #Yes H and bare > no H 
#yes bare > no bare  

#So again this confirms that Glomeromycota (AMF) were higher in 
#successful Horkelia sites than non-successful Horkelia sites. 
#This is true for both Bare and Horeklia soils, suggesting the AMF
#were there before (and thus contributed to) the Horkelia seedling
#establishment. 


#Calculate composition of entire community for results 
#calculate relative abundances of major phlya across all samples 
ps.bac.phy <- microbiome::aggregate_top_taxa(physeq.bac.nona.phylum, "Phylum", top = 10)
ps.bac.phy.rel <- microbiome::transform(ps.bac.phy, "compositional")
bacterial_rel<-as.data.frame(ps.bac.phy.rel@otu_table)
bacterial_rel<-bacterial_rel%>%mutate(tot=rowMeans(.))%>%select(tot)

ps.f.phy <- microbiome::aggregate_top_taxa(physeq.fungi.nona.phylum, "Phylum", top = 7)
ps.f.phy.rel <- microbiome::transform(ps.f.phy, "compositional")
fungal_rel<-as.data.frame(ps.f.phy.rel@otu_table)
fungal_rel<-fungal_rel%>%mutate(tot=rowMeans(.))%>%select(tot)

kable(bacterial_rel)
kable(fungal_rel)

getPalette = colorRampPalette(brewer.pal(5, "Set1"))
speciesPalette = getPalette(9)

#bacteria
physeq.bac.nona.phylum = subset_taxa(physeq.bac, Phylum != "NA") 

taxa_names(physeq.bac.nona.phylum) <- gsub(taxa_names(physeq.bac.nona.phylum), pattern = "p_", replacement = "")

#remove phylum prefix 
ps<-physeq.bac.nona.phylum
tax_table(ps) <- gsub(tax_table(ps), pattern = "p_", replacement = "")

#Microbiome plots 
plotdatB <- microbiome::aggregate_top_taxa(ps,"Phylum", top =8)
plotdatB.rel <- microbiome::transform(plotdatB, "compositional")

plotdatB.relH<-subset_samples(plotdatB.rel, Site=="1"|Site=="9"|Site=="10") 
plotdatB.relAS<-subset_samples(plotdatB.rel, Site=="5")
plotdatB.relAU<-subset_samples(plotdatB.rel, Site=="3"|Site=="6"|Site=="7"|Site=="8"|Site=="11")

#change other value to grey 
speciesPalette[9]<-"#666666"

fig4b1<-microbiome::plot_composition(plotdatB.relH,
                                     average_by = "Site")+ 
  ggtitle("Home soils") + xlab(" ")+ ylab("Relative abundance")+
  scale_fill_manual(values= speciesPalette)+ 
  theme( #axis.ticks.y =element_blank(), 
    #axis.text.x =element_blank(),
    #axis.text.y =element_blank(), 
    plot.title = element_text(size = 10), 
    legend.position = "none")
fig4b2<-microbiome::plot_composition(plotdatB.relAU,
                                     average_by = "Site")+ 
  ggtitle("Away soils unsuccessful") + xlab(" ")+ ylab("")+
  scale_fill_manual(values= speciesPalette)+ 
  theme( axis.ticks.y =element_blank(), 
         #axis.text.x =element_blank(),
         axis.text.y =element_blank(), 
         plot.title = element_text(size = 10), 
         legend.position = "none")

fig4b3<-microbiome::plot_composition(plotdatB.relAS,
                                     average_by = "Site")+ 
  ggtitle("Away soils successful") + xlab(" ")+ ylab(" ")+
  scale_fill_manual(values= speciesPalette)+ 
  theme( axis.ticks.y =element_blank(), 
         #axis.text.x =element_blank(),
         axis.text.y =element_blank(), 
         plot.title = element_text(size = 10))

#fungi
physeq.fungi.nona.phylum = subset_taxa(physeq.fungi, Phylum != "NA") 

#Microbiome plots 
plotdatF <- microbiome::aggregate_top_taxa(physeq.fungi.nona.phylum,"Phylum", top =5)
plotdatF.rel <- microbiome::transform(plotdatF, "compositional")

plotdatF.relH<-subset_samples(plotdatF.rel, Site=="1"|Site=="9"|Site=="10") 
plotdatF.relAS<-subset_samples(plotdatF.rel, Site=="5")
plotdatF.relAU<-subset_samples(plotdatF.rel, Site=="3"|Site=="6"|Site=="7"|Site=="8"|Site=="11")

#change other value to grey 
speciesPalette[6]<-"#666666"

fig4a1<-microbiome::plot_composition(plotdatF.relH,
                                     average_by = "Site")+ 
  ggtitle("Home soils") + xlab(" ")+ ylab("Relative abundance")+
  scale_fill_manual(values= speciesPalette)+ 
  theme( #axis.ticks.y =element_blank(), 
    #axis.text.x =element_blank(),
    #axis.text.y =element_blank(), 
    plot.title = element_text(size = 10), 
    legend.position = "none")
fig4a2<-microbiome::plot_composition(plotdatF.relAU,
                                     average_by = "Site")+ 
  ggtitle("Away soils unsuccessful") + xlab(" ")+ ylab("")+
  scale_fill_manual(values= speciesPalette)+ 
  theme( axis.ticks.y =element_blank(), 
         #axis.text.x =element_blank(),
         axis.text.y =element_blank(), 
         plot.title = element_text(size = 10), 
         legend.position = "none")

fig4a3<-microbiome::plot_composition(plotdatF.relAS,
                                     average_by = "Site")+ 
  ggtitle("Away soils successful") + xlab(" ")+ ylab(" ")+
  scale_fill_manual(values= speciesPalette)+ 
  theme( axis.ticks.y =element_blank(), 
         #axis.text.x =element_blank(),
         axis.text.y =element_blank(), 
         plot.title = element_text(size = 10))


pdf("Fig4.pdf", width = 8, height = 12)
ggpubr::ggarrange(fig4a1, fig4a2, fig4a3, fig4b1, fig4b2, fig4b3, ncol=3, nrow = 2)
dev.off()






