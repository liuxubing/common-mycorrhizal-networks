#########################################################################
###
### R script to run analyses of the in-growth core GROWTH experiment 2017
###
###
### History:
### - 18-08-2018 file created
### - 05-11-2018 version for final test
### - 20-04-2019 minor corrections
### - 28-04-2019 final revision

#########################################################################
###
### Required R libraries and extra software
###

library(ggplot2)
library(gridExtra)
library(lmerTest)
library(lme4)

options(digits=4)        # set number of post-decimal digits in output to 4

rm(list=ls())            # clear workspace


######################################################################
###
### Prepare all files used in this script
### See Data Availability to download the data files
###

Biomass<-read.csv("Tree_Seedling_Growth.csv",header=T)
PotEnv<-read.csv("In-growth_Core_Enviroment.csv",header=T)

Biomass$Pot<-as.factor(Biomass$Pot)
Biomass$MeshSizeNum<-as.factor(Biomass$MeshSizeNum)
Biomass$Diversity<-as.factor(Biomass$Diversity)
Biomass$SpCode<-as.factor(Biomass$SpCode)
Biomass$ID<-as.factor(Biomass$ID)

PotEnv$Pot<-as.factor(PotEnv$Pot)
PotEnv$Diversity<-as.factor(PotEnv$Diversity)

summary(Biomass)
summary(PotEnv)

monobiomass<-subset(Biomass,Diversity==1)
mixbiomass<-subset(Biomass,Diversity==6)

# Note that in the paper "Soil fungal networks maintain local dominance of ectomycorrhizal trees"
# we focused on singleton pots, i.e. monobiomass

##############################################################################################
###
### Compare biomass between monoculture cores with 35um vs 0.5um for each species at each site 
### 

## Calculate mycorrhizal hyphae effect on real biomass

# For each species (monoculture treatment) at each site
# biomass.effect: ln(mean biomass 35um / mean biomass 0.5um) 

biomass.real<-function(Biomass,n=12){
  
  biomass.effect<-vector()
  biomassAB<-vector()
  
  N <- length(Biomass)/n
  
  for (i in 1:N){
    n1 <- 1+(i-1)*n
    n2 <- i*n
    biomassAB<-Biomass[n1:n2]
    biomass.effect[i]<-log(mean(biomassAB[(n/2+1):n],na.rm = T)/mean(biomassAB[1:(n/2)],na.rm = T))
  }
    
  return(biomass.effect)
}

realbiomass<- biomass.real(monobiomass$TotalBiomass)
realbiomass

## Estimate SE of the hyphae effect on real biomass using bootstrap

# n - num of experimental replications (seedling individuals) in 35um + 0.5um cores
# N - num of hyphae effect index (length(TotalBiomass)/n)
# M - num of bootsrap cycles

biomass.bootstrap<-function(Biomass,n=12,M=9999){
  
  N <- length(Biomass)/n
  
  biomassAB<-vector()
  biomassbootsrap<-vector()
  biomass.mean<-vector()
  biomass.CI95<-matrix(NA,N,2)
  
  for (i in 1:N){
    
    n1 <- 1+(i-1)*n
    n2 <- i*n
    
   for (j in 1:M){
    biomassAB <- Biomass[n1:n2]
    
    biomassbootsrap[j] <- log(mean(sample(biomassAB[(n/2+1):n],n/2,replace=T),na.rm = T)
                              /mean(sample(biomassAB[1:(n/2)],n/2,replace = T),na.rm = T))
    biomassbootsrap[M+1] <- log(mean(biomassAB[(n/2+1):n],na.rm = T)
                                /mean(biomassAB[1:(n/2)],na.rm = T))
    
    biomassbootsrap.mean <- mean(biomassbootsrap,na.rm = T)
    biomassbootsrap.CI95 <- quantile(biomassbootsrap, probs = c(0.025, 0.975),na.rm = T)
    
   }
    
    biomass.mean[i] <- biomassbootsrap.mean 
    biomass.CI95[i,] <- biomassbootsrap.CI95
  }
  
  bootbiomass<-data.frame(bootbiomass=biomass.mean,bootbiomass.CI95=biomass.CI95)
  return(bootbiomass)
}

bootbiomass<-biomass.bootstrap(monobiomass$TotalBiomass,12,9999)
bootbiomass

# Combine real biomasseffects and bootbiomass effects
biomass.real.boot<-data.frame(realbiomass,bootbiomass)
biomass.real.boot

# Read header of pot treatments and get whole data frame
monoculture <- read.csv("MonoculturePotTreatments.csv",header = T)
monoculture$SpCode2<-as.numeric(monoculture$SpCode2)
monoculture$Diversity<-as.factor(monoculture$Diversity)

biomass.hyphaeeffect<-data.frame(monoculture,biomass.real.boot)
biomass.hyphaeeffect

write.csv(biomass.hyphaeeffect,file = "biomass.hyphaeeffect.csv")

## plot the results of all sites

# plot results of ECM sites

site.ECM <- subset(biomass.hyphaeeffect, SiteMycorrhizalType == "ECM")
site.ECM$Species <- reorder(site.ECM$Species,site.ECM$SpCode2)
site.ECM$Site <- reorder(site.ECM$Site,site.ECM$SiteCode)

ECM.plot <- ggplot(site.ECM, aes(x=Species, y=realbiomass,fill=SpMycorrhizalType)) +
  geom_bar(colour="black",stat="identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin=bootbiomass.CI95.1,ymax=bootbiomass.CI95.2),width=0.1) +
  facet_grid(Site~.) +
  scale_fill_manual(values = c("#FFDDDD","#CCEEFF","red"),limits=c("ECM","AM","Conspecific")) +
  xlab("Seedling species") +
  ylab("Mycorrhizal hyphae effect on biomass \n ln(35 µm cores / 0.5 µm cores)") +
  ylim(-1.0,1.5) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  ggtitle("ECM dominated sites") +
  theme(plot.title = element_text(size = 16)) + 
  theme(plot.title = element_text(hjust = 0.5,face="bold")) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(axis.text = element_text(size=rel(1.6)), axis.title = element_text(size=rel(1.8)), 
        strip.text = element_text(size=rel(1.8))) +
  theme(legend.position=c(0.99,0.99), legend.justification=c(1,1),legend.direction="horizontal",legend.text = element_text(size = 12)) +
  guides(fill=guide_legend(title = NULL))
  
# plot results of AM sites

site.AM <- subset(biomass.hyphaeeffect, SiteMycorrhizalType == "AM")
site.AM$Species <- reorder(site.AM$Species,site.AM$SpCode2)
site.AM$Site <- reorder(site.AM$Site,site.AM$SiteCode)

AM.plot <- ggplot(site.AM, aes(x=Species, y=realbiomass,fill=SpMycorrhizalType)) +
  geom_bar(colour="black",stat="identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin=bootbiomass.CI95.1,ymax=bootbiomass.CI95.2),width=0.1) +
  facet_grid(Site~.) +
  scale_fill_manual(values = c("#FFDDDD","#CCEEFF","red"),limits=c("ECM","AM","Conspecific")) +
  xlab("Seedling species") +
  ylab("Mycorrhizal hyphae effect on biomass \n ln(35 µm cores / 0.5 µm cores)") +
  geom_hline(yintercept = 0) +
  theme_bw() +
  ggtitle("AM dominated sites") +
  ylab("") +
  theme(plot.title = element_text(size = 16)) +
  theme(plot.title = element_text(hjust = 0.5,face="bold")) +
  theme(axis.text = element_text(size=rel(1.6)), axis.title = element_text(size=rel(1.8)), 
        axis.text.y=element_blank(),strip.text = element_text(size=rel(1.8))) +
  theme(legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# plot ECM and AM sites to a single figure
# library("gridExtra")

grid.arrange(ECM.plot, AM.plot, nrow=1, ncol=2)
# PNG 1420 * 650, PDF 30.5 * 14.2 

############################################################################################
###
### Compare biomass of each focal species between conspecific site vs heterospecific sites 
###

# Subset monocluture pots and remove NAs

monobiomass1 <- subset(monobiomass, TotalBiomass!="NA")
monobiomass1$SiteType <- as.factor(monobiomass1$Site==monobiomass1$Species)

focalspecies <- data.frame(species=c("Cfab","Cfis","Efen","Ssup","Ccon","Asty"))

### Calculate mean biomass (se) in conspecific site vs heterospecific sites for all species

SiteMeanBiomass <- data.frame()

for (i in 1:length(focalspecies$species)) {
   
   sitebiomass <- subset(monobiomass1, Species==focalspecies$species[i])
   
   SiteMeanBiomassi <- aggregate(sitebiomass$TotalBiomass,list(SiteType=sitebiomass$SiteType,
                                                            MeshSize=sitebiomass$MeshSize),mean)
   SiteMeanBiomassSDi<- aggregate(sitebiomass$TotalBiomass,list(SiteType=sitebiomass$SiteType,
                                                             MeshSize=sitebiomass$MeshSize),sd)
   SiteMeanBiomassLengthi <- aggregate(sitebiomass$TotalBiomass,list(SiteType=sitebiomass$SiteType,
                                                                 MeshSize=sitebiomass$MeshSize),length)
   SiteMeanBiomassi$se <- SiteMeanBiomassSDi$x/SiteMeanBiomassLengthi$x
   SiteMeanBiomassi$n <- SiteMeanBiomassLengthi$x
   SiteMeanBiomassi$species <- focalspecies$species[i]
   
   SiteMeanBiomass<-rbind(SiteMeanBiomass, SiteMeanBiomassi)
   
   rm(SiteMeanBiomassi,SiteMeanBiomassSDi,SiteMeanBiomassLengthi)

}

print(SiteMeanBiomass)

SiteMeanBiomass$SpCode2 <- rep(1:6, each=4)
SiteMeanBiomass$species <- reorder(SiteMeanBiomass$species,SiteMeanBiomass$SpCode2)
SiteMeanBiomass$MycorrhizalType <- rep(c("ECM","AM"), each=12)


### plot the biomass conspecific site vs heterospecific sites for each species
## aes(x=MeshSize, y=x, fill=rev(SiteType)) ##

# bar plot results of ECM species

SiteMeanBiomass.ECM <- subset(SiteMeanBiomass, MycorrhizalType == "ECM")

SiteMeanBiomass.ECM.plot <- ggplot(SiteMeanBiomass.ECM,aes(x=MeshSize, y=x, fill=rev(SiteType))) +
  geom_bar(stat = "identity",colour="black", width=0.5, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin=x-se, ymax=x+se),position = position_dodge(0.7), width=0.2) +
  #geom_text(aes(y=-Inf,vjust=-1.5,label=paste("n =", format(SiteMeanBiomass.ECM$n))),position = position_dodge(0.7)) +
  facet_grid(species~.,scales = "free_y") +
  theme_bw() +
  xlab("Mesh Size") +
  ylab("Seedling biomass (g)") +
  ggtitle("ECM tree species") +
  theme(plot.title=element_text(hjust = 0.5, face="bold")) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(legend.position=c(0.99,0.99), legend.justification=c(1,1), legend.direction="horizontal",legend.text = element_text(size = 10)) +
  theme(axis.text = element_text(size=rel(1.2)), axis.title = element_text(size=rel(1.2)), 
        strip.text = element_text(size=rel(1.2))) +
  scale_x_discrete(limits=c("B","A"),labels=c("35 µm", "0.5 µm")) +
  scale_fill_manual(values=c("white","grey"), labs(""),labels=c("Conspecific site", "Heterospecific sites"))

# bar plot results of AM species

SiteMeanBiomass.AM <- subset(SiteMeanBiomass, MycorrhizalType == "AM")

SiteMeanBiomass.AM.plot <- ggplot(SiteMeanBiomass.AM,aes(x=MeshSize, y=x, fill=rev(SiteType))) +
  geom_bar(stat="identity",colour="black", width=0.5, position=position_dodge(0.7)) +
  geom_errorbar(aes(ymin=x-se, ymax=x+se),position = position_dodge(0.7), width=0.2) +
  # geom_text(aes(y=-Inf,vjust=-1.5,label=paste("n =", format(SiteMeanBiomass.AM$n))),position = position_dodge(0.7)) +
  facet_grid(species~.,scales="free_y") +
  theme_bw() +
  xlab("Mesh Size") +
  ylab("") +
  ggtitle("AM tree species") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(legend.position="none") +
  theme(axis.text = element_text(size=rel(1.2)), axis.title = element_text(size=rel(1.2)), 
        strip.text = element_text(size=rel(1.2))) +
  scale_x_discrete(limits=c("B","A"),labels=c("35 µm", "0.5 µm")) +
  scale_fill_manual(values=c("white","grey"), labs("Site Type"),labels=c("Conspecific", "Heterospecific"))

# plot ECM and AM sites to a single figure

# library("gridExtra")

grid.arrange(SiteMeanBiomass.ECM.plot, SiteMeanBiomass.AM.plot, nrow=1, ncol=2)
# png 1440 * 838

### plot the biomass conspecific site vs heterospecific sites for each species
## aes(x=SiteType, y=x, fill=rev(MeshSize)) ##

# plot results of ECM species

SiteMeanBiomass.ECM <- subset(SiteMeanBiomass, MycorrhizalType == "ECM")

SiteMeanBiomass.ECM.plot <- ggplot(SiteMeanBiomass.ECM,aes(x=SiteType, y=x, fill=rev(MeshSize))) +
  geom_bar(stat = "identity",colour="black", width=0.5, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin=x-se, ymax=x+se),position = position_dodge(0.7), width=0.2) +
  geom_text(aes(y=-Inf,vjust=-1.5,label=paste("n =", format(SiteMeanBiomass.ECM$n))),position = position_dodge(0.7)) +
  facet_grid(species~.,scales = "free_y") +
  theme_bw() +
  xlab("Site Type") +
  ylab("Seedling biomass (g)") +
  ggtitle("ECM tree species") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(legend.position = c(0.99,0.99), legend.justification = c(1,1), legend.direction = "horizontal") +
  scale_x_discrete(limits=c("TRUE","FALSE"),labels=c("Conspecific site", "Heterospecific sites")) +
  scale_fill_manual(values = c("white","grey"), labs("Mesh Size"),labels=c("35 µm", "0.5 µm"))

# plot results of AM species

SiteMeanBiomass.AM <- subset(SiteMeanBiomass, MycorrhizalType == "AM")

SiteMeanBiomass.AM.plot <- ggplot(SiteMeanBiomass.AM,aes(x=SiteType, y=x, fill=rev(MeshSize))) +
  geom_bar(stat = "identity",colour="black", width=0.5, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin=x-se, ymax=x+se),position = position_dodge(0.7), width=0.2) +
  geom_text(aes(y=-Inf,vjust=-1.5,label=paste("n =", format(SiteMeanBiomass.AM$n))),position = position_dodge(0.7)) +
  facet_grid(species~.,scales = "free_y") +
  theme_bw() +
  xlab("Site Type") +
  ylab("") +
  ggtitle("AM tree species") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  scale_x_discrete(limits=c("TRUE","FALSE"),labels=c("Conspecific site", "Heterospecific sites")) +
  scale_fill_manual(values = c("white","grey"), labs("Mesh Size"),labels=c("35 µm", "0.5 µm"))

# plot ECM and AM sites to a single figure

grid.arrange(SiteMeanBiomass.ECM.plot, SiteMeanBiomass.AM.plot, nrow=1, ncol=2)


########################################################################################################
###
### Compare biomass between singleton cores with conspecific vs heterospecific for ECM and AM respectively
###

## Calculate site effect on real biomass

# subset the monoculture cores
monobiomass2<-subset(Biomass,Diversity==1)
monobiomass2$SiteType <- as.factor(monobiomass2$Site==monobiomass2$Species)

# get the results of own site vs other sites: SiteEffect = ln(conspecific site / heterospecific sites) 
# for each species with each mesh size

Site.hyphaeeffect <- data.frame(species=rep(c("Cfab","Cfis","Efen","Ssup","Ccon","Asty"),each=2),
                                MeshSize=c("A","B"),MycorrhizalType=rep(c("ECM","AM"),each=6))

biomass.real.2<-function(Biomass){
  
  biomass<-vector()
  biomassAB<-vector()
  biomassSpecies<-vector()
  
  for (i in 1:length(Site.hyphaeeffect$species)) {
    
    biomassSpecies<-subset(Biomass, Species==Site.hyphaeeffect$species[i])
    biomassAB<-subset(biomassSpecies, MeshSize==Site.hyphaeeffect$MeshSize[i])
    
    biomass[i]<-log(mean(biomassAB$TotalBiomass[biomassAB$SiteType=="TRUE"],na.rm = T)
                    /mean(biomassAB$TotalBiomass[biomassAB$SiteType=="FALSE"],na.rm = T))
  }
  
  print(biomass)
}

Site.hyphaeeffect$SiteEffect <- biomass.real.2(monobiomass2)
Site.hyphaeeffect

# combine the results into two mycorrhizal types: ECM vs AM

Site.hyphaeeffect.ECMvsAM <- aggregate(SiteEffect~MycorrhizalType+MeshSize,data=Site.hyphaeeffect,mean)
Site.hyphaeeffect.ECMvsAM$se <- aggregate(SiteEffect~MycorrhizalType+MeshSize,
                                          data=Site.hyphaeeffect,sd)$SiteEffect/sqrt(3)
Site.hyphaeeffect.ECMvsAM  # SiteEffect = mean of ln(conspecific site / heterospecific sites) for each mycorrhizal type with se

# plot the own site vs other sites effects ECM vs AM, 35 vs 0.5

ggplot(Site.hyphaeeffect.ECMvsAM, aes(x=MeshSize, y=SiteEffect,fill=rev(MycorrhizalType))) +
  geom_bar(colour="black",stat="identity", width=0.5,position=position_dodge(0.7)) +
  geom_errorbar(aes(ymin=SiteEffect-se,ymax=SiteEffect+se),width=0.2, position=position_dodge(0.7)) +
  xlab("Mesh Size") +
  ylab("Plant-soil feedback on biomass \n ln(Conspecific site / Heterospecific sites)") +
  ylim(-1.2,0.5) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(legend.position = c(0.99,0.99), legend.justification = c(1,1), legend.direction = "horizontal") +
  theme(axis.text = element_text(size=rel(1.2)), axis.title = element_text(size=rel(1.2)), 
        strip.text = element_text(size=rel(1.2))) +
  guides(fill=guide_legend(title = NULL)) +
  scale_x_discrete(limits=c("B","A"),labels=c("35 µm", "0.5 µm")) +
  scale_fill_manual(values = c("white","grey"), labs(""),labels=c("ECM", "AM"))

ggsave("Fig3.pdf")

######################################################################
###
### Regression Models for seedling growth 
###

### linear models

# For the original biomass data

lm1 <- aov(TotalBiomass ~ SiteMycorrhizalType * SpMycorrhizalType + MeshSize 
           + Species + Diversity + Site + Pot, data = Biomass)
summary(lm1)

# For the hyphae effect data (ln(35µm cores / 0.5µm cores))

lm2 <- aov(realbiomass ~ SiteMycorrhizalType * SpMycorrhizalType + Species + Site,
           data=biomass.hyphaeeffect)
summary(lm2)

### Mixed-effects Models

# Mixed-model analysis with 'lmer' (library: lme4)

mm1 <- lmer(TotalBiomass ~ SiteMycorrhizalType * SpMycorrhizalType + MeshSize 
            + Diversity + (1|Species) + (1|Site) + (1|Pot), data = Biomass)
anova(mm1, ddf="lme4")
summary(mm1)

mm2 <- lmer(realbiomass~SiteMycorrhizalType*SpMycorrhizalType + (1|Site) + (1|Species),
            data=biomass.hyphaeeffect)
anova(mm2, ddf="lme4")
summary(mm2)

# select the best fitted model by AIC

lm3<-step(lm1)
summary(lm3)

lm4<-step(lm2)
summary(lm4)

step(mm1)
step(mm2)

### get the results of best model

# For the original biomass data
mm.best.1 <- lmer(TotalBiomass ~ SiteMycorrhizalType * SpMycorrhizalType + MeshSize + 
                    (1 | Species) + (1 | Site) + (1 | Pot), data=Biomass)
summary(mm.best.1)

# For the hyphae effect data (ln(35µm cores / 0.5µm cores))

mm.best.2 <- lmer(realbiomass~SiteMycorrhizalType*SpMycorrhizalType + (1|Site) + (1|Species),
                  data=biomass.hyphaeeffect)
summary(mm.best.2)

lm.best.2 <- lm(realbiomass ~ SiteMycorrhizalType * SpMycorrhizalType + Species + Site,
                data=biomass.hyphaeeffect)
summary.aov(lm2)


####################################################################################
###
### Compare enviromental variables among pots with different mesh sizes (Figure S1)
###

## Soil temperature

pairwise.t.test(PotEnv$Temp20170526,PotEnv$MeshSize)
Temp0526.aov<-aov(Temp20170526~MeshSize+Site,data=PotEnv)
summary(Temp0526.aov)

pairwise.t.test(PotEnv$Temp20170926,PotEnv$MeshSize)
Temp0926.aov<-aov(Temp20170926~MeshSize+Site,data=PotEnv)
summary(Temp0926.aov)

pairwise.t.test(PotEnv$Temp20180117,PotEnv$MeshSize)
Temp0117.aov<-aov(Temp20180117~MeshSize+Site,data=PotEnv)
summary(Temp0117.aov)

## Soil Humidity

pairwise.t.test(PotEnv$Humidity20170526,PotEnv$MeshSize)
Humidity0526.aov<-aov(Humidity20170526~MeshSize+Site,data=PotEnv)
summary(Humidity0526.aov)

pairwise.t.test(PotEnv$Humidity20170926,PotEnv$MeshSize)
Humidity0926.aov<-aov(Humidity20170926~MeshSize+Site,data=PotEnv)
summary(Humidity0926.aov)

pairwise.t.test(PotEnv$Humidity20180117,PotEnv$MeshSize)
Humidity0117.aov<-aov(Humidity20180117~MeshSize+Site,data=PotEnv)
summary(Humidity0117.aov)

Humidity0117.aov<-aov(Humidity20180117~MeshSize+Site+MycorrhizalType+Diversity,data=PotEnv)
summary(Humidity0117.aov)

## Plot the results of enviromental variables

## bar plot ## 

# Calculate the mean (se) Temp and Humidity for different MeshSize (35 vs 0.5) at each site

PotMeanTemp <- aggregate(Temp20180117~Site+MeshSize,data=PotEnv,mean)
PotMeanTempSD <- aggregate(Temp20180117~Site+MeshSize,data=PotEnv,sd)

PotMeanHumidity <- aggregate(PotEnv$Humidity20180117*100,list(Site=PotEnv$Site,MeshSize = PotEnv$MeshSize),mean)
PotMeanHumiditySD <- aggregate(PotEnv$Humidity20180117*100,list(Site=PotEnv$Site,MeshSize = PotEnv$MeshSize),sd)

PotMeanEnv <- data.frame(Site=PotMeanTemp$Site,MeshSize=PotMeanTemp$MeshSize,PotMeanTemp=PotMeanTemp$Temp20180117,
                        PotMeanTempSE=PotMeanTempSD$Temp20180117/sqrt(9),PotMeanHumidity=PotMeanHumidity$x,
                        PotMeanHumiditySE=PotMeanHumiditySD$x/sqrt(9))

Temp.plot <- ggplot(PotMeanEnv,aes(x=Site,y=PotMeanTemp,fill=rev(MeshSize))) +
  geom_bar(stat = "identity",colour="black", width=0.7, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin=PotMeanTemp-PotMeanTempSE, ymax=PotMeanTemp+PotMeanTempSE),
                position = position_dodge(0.9), width=0.2) +
  xlab("Site") +
  ylab("Pot Mean Temperature (°C)") +
  ylim(0,21) +
  ggtitle("Effect of Mesh Size on Pot Temperature at Each Site") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,face="bold")) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=rel(1.2)), axis.title = element_text(size=rel(1.2)), 
        strip.text = element_text(size=rel(1.2))) +
  scale_fill_manual(values = c("white","grey"))

Humidity.plot <- ggplot(PotMeanEnv,aes(x=Site,y=PotMeanHumidity,fill=rev(MeshSize))) +
  geom_bar(stat = "identity",colour="black", width=0.7, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin=PotMeanHumidity-PotMeanHumiditySE, ymax=PotMeanHumidity+PotMeanHumiditySE),
                position = position_dodge(0.9), width=0.2) +
  xlab("") +
  ylab("Pot Mean Moisture Content (%)") +
  ylim(0,25) +
  ggtitle("Effect of Mesh Size on Pot Moisture Content at Each Site") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,face="bold")) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(legend.position = c(0.01,0.88), legend.justification = c(0,0), legend.direction = "horizontal") +
  theme(axis.text = element_text(size=rel(1.2)), axis.title = element_text(size=rel(1.2)), 
        strip.text = element_text(size=rel(1.2))) +
  scale_fill_manual(values = c("white","grey"), labs("Mesh Size"),labels=c("35 µm","0.5 µm"))

# plot Temp and Humidity sites to a single figure

grid.arrange(Humidity.plot, Temp.plot,nrow=2, ncol=1)
# png 1440*838  pdf 8.3*11.8

## box-plot ##

Temp.plot <- ggplot(PotEnv,aes(x=Site,y=Temp20180117,fill=rev(MeshSize))) +
  geom_boxplot()+
  xlab("Site") +
  ylab("Pot Mean Temperature (°C)") +
  ggtitle("Effect of Mesh Size on Pot Temperature at Each Site") +
  ylim(15,21) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,face="bold")) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=rel(1.2)), axis.title = element_text(size=rel(1.2)), 
        strip.text = element_text(size=rel(1.2))) +  
  scale_fill_manual(values = c("white","grey")); Temp.plot

Humidity.plot <- ggplot(PotEnv,aes(x=Site,y=Humidity20180117*100,fill=rev(MeshSize))) +
  geom_boxplot()+
  xlab("") +
  ylab("Pot Mean Moisture Content (%)") +
  ggtitle("Effect of Mesh Size on Pot Moisture Content at Each Site") +
  ylim(0,31) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,face="bold")) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(legend.position = c(0.01,0.88), legend.justification = c(0,0), legend.direction = "horizontal") +
  theme(axis.text = element_text(size=rel(1.2)), axis.title = element_text(size=rel(1.2)), 
        strip.text = element_text(size=rel(1.2))) +
  scale_fill_manual(values = c("white","grey"), labs("Mesh Size"),labels=c("35 µm","0.5 µm")); Humidity.plot

# plot Temp and Humidity sites to a single figure

grid.arrange(Humidity.plot, Temp.plot,nrow=2, ncol=1)
# png 875*838

#####################################################################################################
###
### Compare root colonization of each focal species between conspecific site vs heterospecific sites
###

# import data

colonization <- read.csv("Tree_Seedling_Root_Colonization.csv",header = T)

focalspecies <- data.frame(species=c("Cfab","Cfis","Efen","Ssup","Ccon","Asty"))

### Calculate mean colonization rates (se) in conspecific site vs heterospecific sites for all species

SiteMeanColonization <- data.frame()

for (i in 1:length(focalspecies$species)) {
  
  sitecolonization <- subset(colonization, Species==focalspecies$species[i])
  
  SiteMeanColonizationi <- aggregate(sitecolonization$Rate,list(SiteType=sitecolonization$SiteType,
                                                              MeshSize=sitecolonization$MeshSize),mean)
  SiteMeanColonizationSDi<- aggregate(sitecolonization$Rate,list(SiteType=sitecolonization$SiteType,
                                                                 MeshSize=sitecolonization$MeshSize),sd)
  SiteMeanColonizationLengthi <- aggregate(sitecolonization$Rate,list(SiteType=sitecolonization$SiteType,
                                                                      MeshSize=sitecolonization$MeshSize),length)
  SiteMeanColonizationi$se <- SiteMeanColonizationSDi$x/SiteMeanColonizationLengthi$x
  SiteMeanColonizationi$n <- SiteMeanColonizationLengthi$x
  SiteMeanColonizationi$species <- focalspecies$species[i]
  
  SiteMeanColonization<-rbind(SiteMeanColonization, SiteMeanColonizationi)
  
  rm(SiteMeanColonizationi,SiteMeanColonizationSDi,SiteMeanColonizationLengthi)
  
}

print(SiteMeanColonization)

SiteMeanColonization$SpCode2 <- rep(1:6, each=4)
SiteMeanColonization$species <- reorder(SiteMeanColonization$species,SiteMeanColonization$SpCode2)
SiteMeanColonization$MycorrhizalType <- rep(c("ECM","AM"), each=12)


### plot the root colonization conspecific site vs heterospecific sites for each species
## aes(x=MeshSize, y=x, fill=rev(SiteType)) ##

# plot results of ECM species

SiteMeanColonization.ECM <- subset(SiteMeanColonization, MycorrhizalType == "ECM")

SiteMeanColonization.ECM.plot <- ggplot(SiteMeanColonization.ECM,aes(x=MeshSize, y=x, fill=SiteType)) +
  geom_bar(stat = "identity",colour="black", width=0.5, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin=x-se, ymax=x+se),position = position_dodge(0.7), width=0.2) +
  #geom_text(aes(y=-Inf,vjust=-1.5,label=paste("n =", format(SiteMeanColonization.ECM$n))),position = position_dodge(0.7)) +
  facet_grid(species~.,scales = "free_y") +
  theme_bw() +
  xlab("Mesh Size") +
  ylab("Seedling root colonization") +
  ggtitle("ECM tree species") +
  theme(plot.title=element_text(hjust = 0.5, face="bold")) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(legend.position="none") +
  theme(axis.text = element_text(size=rel(1.2)), axis.title = element_text(size=rel(1.2)), 
        strip.text = element_text(size=rel(1.2))) +
  scale_x_discrete(limits=c("B","A"),labels=c("35 µm", "0.5 µm")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_manual(values=c("white","grey"), labs(""),labels=c("Conspecific site", "Heterospecific sites"))

# plot results of AM species

SiteMeanColonization.AM <- subset(SiteMeanColonization, MycorrhizalType == "AM")

SiteMeanColonization.AM.plot <- ggplot(SiteMeanColonization.AM,aes(x=MeshSize, y=x, fill=SiteType)) +
  geom_bar(stat="identity",colour="black", width=0.5, position=position_dodge(0.7)) +
  geom_errorbar(aes(ymin=x-se, ymax=x+se),position = position_dodge(0.7), width=0.2) +
  #geom_text(aes(y=-Inf,vjust=-1.5,label=paste("n =", format(SiteMeanColonization.AM$n))),position = position_dodge(0.7)) +
  facet_grid(species~.,scales="free_y") +
  theme_bw() +
  xlab("Mesh Size") +
  ylab("") +
  ggtitle("AM tree species") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(legend.position=c(0.99,0.99), legend.justification=c(1,1), legend.direction="horizontal",legend.text = element_text(size = 14)) +
  theme(axis.text = element_text(size=rel(1.2)), axis.title = element_text(size=rel(1.2)), 
        strip.text = element_text(size=rel(1.2))) +
  scale_x_discrete(limits=c("B","A"),labels=c("35 µm", "0.5 µm")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_manual(values=c("white","grey"), labs(""),labels=c("Conspecific", "Heterospecific"))

# plot ECM and AM sites to a single figure

# library("gridExtra")

grid.arrange(SiteMeanColonization.ECM.plot, SiteMeanColonization.AM.plot, nrow=1, ncol=2)
# png 1440 * 838


### plot the root colonization conspecific site vs heterospecific sites for each species
## aes(x=SiteType, y=x, fill=rev(MeshSize)) ##

# plot results of ECM species

SiteMeanColonization.ECM <- subset(SiteMeanColonization, MycorrhizalType == "ECM")

SiteMeanColonization.ECM.plot <- ggplot(SiteMeanColonization.ECM,aes(x=SiteType, y=x, fill=rev(MeshSize))) +
  geom_bar(stat = "identity",colour="black", width=0.5, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin=x-se, ymax=x+se),position = position_dodge(0.7), width=0.2) +
  #geom_text(aes(y=-Inf,vjust=-1.5,label=paste("n =", format(SiteMeanColonization.ECM$n))),position = position_dodge(0.7)) +
  facet_grid(species~.,scales = "free_y") +
  theme_bw() +
  xlab("Site Type") +
  ylab("Seedling biomass (g)") +
  ggtitle("ECM tree species") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(legend.position = c(0.99,0.99), legend.justification = c(1,1), legend.direction = "horizontal") +
  scale_x_discrete(limits=c("Con","Hetero"),labels=c("Conspecific site", "Heterospecific sites")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_manual(values = c("white","grey"), labs(""),labels=c("35 µm", "0.5 µm"))

# plot results of AM species

SiteMeanColonization.AM <- subset(SiteMeanColonization, MycorrhizalType == "AM")

SiteMeanColonization.AM.plot <- ggplot(SiteMeanColonization.AM,aes(x=SiteType, y=x, fill=rev(MeshSize))) +
  geom_bar(stat = "identity",colour="black", width=0.5, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin=x-se, ymax=x+se),position = position_dodge(0.7), width=0.2) +
  #geom_text(aes(y=-Inf,vjust=-1.5,label=paste("n =", format(SiteMeanColonization.AM$n))),position = position_dodge(0.7)) +
  facet_grid(species~.,scales = "free_y") +
  theme_bw() +
  xlab("Site Type") +
  ylab("") +
  ggtitle("AM tree species") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  scale_x_discrete(limits=c("Con","Hetero"),labels=c("Conspecific site", "Heterospecific sites")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_manual(values = c("white","grey"), labs("Mesh Size"),labels=c("35 µm", "0.5 µm"))

# plot ECM and AM sites to a single figure

grid.arrange(SiteMeanColonization.ECM.plot, SiteMeanColonization.AM.plot, nrow=1, ncol=2)



###
### End of script
###
#######################################################