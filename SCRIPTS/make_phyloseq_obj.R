#make phyloseq objects
library(phyloseq)
library(tidyr)
library(vegan)

#BACTERIA
#metadata
metab=read.csv("DATA/meta_phyloseq_Courtney.csv")
metab<-mutate(metab, Hsuccess=
                         ifelse(Site=="5"&Loc=="Away","Y","N"))%>%
  unite(col = "successx", Hsuccess, Spp, remove=F)%>%select(-Code2)
                    
rownames(metab) <- metab$Code# Assign rownames to be Sample IDs
metabx <- sample_data(metab)

#OTU data (rarefied with tax) 
otub<-unzip("DATA/Serp_16s_otu_table_wTaxa_rare253827.zip")
otub = read.csv(otub)
rownames(otub)<-otub$Row.names
otub$Row.names<-NULL
str(otub)
otumatb <- as(as.matrix(otub[, c(8:160)]), "matrix")
OTUB = otu_table(otumatb, taxa_are_rows = TRUE)

#taxa 
taxb<-otub[, c(1:8)]#pull out tax 
taxb$Row.names<-NULL
#determine % unassigned 
unassigned<-subset(taxb, Domain=="Unassigned")
(nrow(unassigned)/nrow(taxb))*100
#PULL OUT ONLY BACTERIA 
taxb<-filter(taxb, Domain=="k__Bacteria")
#make martrix
taxmatb <- as(as.matrix(taxb),"matrix")
TAXB = tax_table(taxmatb)

#make phyloseq object
physeq.bac = phyloseq(OTUB,TAXB, metabx)
#remove singletons 
physeq.bac = prune_taxa(taxa_sums(physeq.bac) > 1, physeq.bac)
save(physeq.bac, metab, file="DATA/physeq.bac.Rdata")

#FUNGI-----
#metadata
metaf = read.csv("DATA/meta_phyloseq_Courtney.csv")
metaf= metaf[, c(1:9)]
metaf$elev<-as.factor(metaf$elev)
metaf<-mutate(metaf, Hsuccess=
                ifelse(Site=="5"&Loc=="Away","Y","N"))%>%
  unite(col = "successx", Hsuccess, Spp, remove=F)
rownames(metaf) <- metaf$Code # Assign rownames to be Sample IDs 
metafx <- sample_data(metaf)

#OTU data (rarefied no tax) 
otu = read.csv("DATA/SerpentineITS.final.csv")
rownames(otu)<-otu$X.OTU.ID
otu$X.OTU.ID<-NULL
otumat <- as(as.matrix(otu), "matrix")
OTU = otu_table(otumat, taxa_are_rows = TRUE)

#taxa
tax<-read.csv("DATA/SerpentineITS.taxonomy.fix.csv")
row.names(tax)<-tax$X
tax$X<-NULL
#determine % unassigned 
unassigned<-subset(tax, is.na(Kingdom))
(nrow(unassigned)/nrow(tax))*100

#PULL OUT ONLY FUNGI
tax<-filter(tax, Kingdom=="Fungi")
#make matrix 
taxmat <- as(as.matrix(tax),"matrix")
TAX = tax_table(taxmat)


#make phyloseq object----
physeq.fungi = phyloseq(OTU,TAX, metafx)
#remove singletons 
physeq.fungi = prune_taxa(taxa_sums(physeq.fungi) > 1, physeq.fungi)

#rarefy fungi 
#rarefy to 100% of the miniumum sample depth in the dataset 
sort(sample_sums(physeq.fungi))#lowest is 22822 

physeq.fungi.rare=ps.rarefied = rarefy_even_depth(physeq.fungi, 
          sample.size=min(sample_sums(physeq.fungi)), replace=F)
sample_sums(physeq.fungi.rare)#2282 2reads per sample 

physeq.fungi<-physeq.fungi.rare
save(physeq.fungi, metaf, file="DATA/physeq.fungi.Rdata")
