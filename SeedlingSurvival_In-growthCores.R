###########################################################################
###
### R script to run analyses of the in-growth core SURVIVAL experiment 2018
###
###
### History:
### - 18-08-2018 file created
### - 05-11-2018 version for final test
### - 28-06-2019 final revision

###########################################################################

library(ggplot2)
library(gridExtra)
library(lmerTest)
library(lme4)

######################################################################
###
### Prepare all files used in this script
### See Data Availability to download the data files
###

SeedlingSurvival<-read.csv("Tree_Seedling_Survival.csv")
#View(SeedlingSurvival)


### GLMM ####

SeedlingSurvival.glmm <- glmer((survival ~ MeshSize * MycorrhizalType + (1 | SPcode) + (1 | Core)), family = binomial, data=SeedlingSurvival)
summary(SeedlingSurvival.glmm)


### two-proportions z-test ###

## Calculate the mean survival rate ± se for each treatment

SeedlingSurvivalmean<-aggregate(survival~SPcode+MeshSize,data=SeedlingSurvival,mean)

# se for binary data sqrt(pq/n)
SeedlingSurvivalmean$se<-sqrt((SeedlingSurvivalmean$survival*(1-SeedlingSurvivalmean$survival))/48)

SeedlingSurvivalmean$MycorrhizalType<-c("AM","AM","ECM","ECM","ECM","ECM","AM","AM")
SeedlingSurvivalmean$SP <-rep(c(7,6,1,2,3,4,8,5),2)
SeedlingSurvivalmean$SPCode <- reorder(SeedlingSurvivalmean$SPcode,SeedlingSurvivalmean$SP)

anova(lm(survival~MycorrhizalType+SPcode*MeshSize,data=SeedlingSurvival))

ECM.SeedlingSurvival<-subset(SeedlingSurvivalmean,MycorrhizalType=="ECM")
AM.SeedlingSurvival<-subset(SeedlingSurvivalmean,MycorrhizalType=="AM")

CFA.L.SeedlingSurvival<-subset(SeedlingSurvival, SPcode=="Cfab"&MeshSize=="L")
CFA.S.SeedlingSurvival<-subset(SeedlingSurvival, SPcode=="Cfab"&MeshSize=="S")
prop.test(x=c(sum(CFA.L.SeedlingSurvival$survival1),sum(CFA.S.SeedlingSurvival$survival1)),n=c(48,48),alternative = "greater")

CFI.L.SeedlingSurvival<-subset(SeedlingSurvival, SPcode=="Cfis"&MeshSize=="L")
CFI.S.SeedlingSurvival<-subset(SeedlingSurvival, SPcode=="Cfis"&MeshSize=="S")
prop.test(x=c(sum(CFI.L.SeedlingSurvival$survival),sum(CFI.S.SeedlingSurvival$survival)),n=c(48,48),alternative = "greater")

CHU.L.SeedlingSurvival<-subset(SeedlingSurvival, SPcode=="Chui"&MeshSize=="L")
CHU.S.SeedlingSurvival<-subset(SeedlingSurvival, SPcode=="Chui"&MeshSize=="S")
prop.test(x=c(sum(CHU.L.SeedlingSurvival$survival),sum(CHU.S.SeedlingSurvival$survival)),n=c(48,48),alternative = "greater")

LHA.L.SeedlingSurvival<-subset(SeedlingSurvival, SPcode=="Lhai"&MeshSize=="L")
LHA.S.SeedlingSurvival<-subset(SeedlingSurvival, SPcode=="Lhai"&MeshSize=="S")
prop.test(x=c(sum(LHA.L.SeedlingSurvival$survival),sum(LHA.S.SeedlingSurvival$survival)),n=c(48,48),alternative = "greater")

CAL.L.SeedlingSurvival<-subset(SeedlingSurvival, SPcode=="Calb"&MeshSize=="L")
CAL.S.SeedlingSurvival<-subset(SeedlingSurvival, SPcode=="Calb"&MeshSize=="S")
prop.test(x=c(sum(CAL.L.SeedlingSurvival$survival),sum(CAL.S.SeedlingSurvival$survival)),n=c(48,48),alternative = "less")

CCO.L.SeedlingSurvival<-subset(SeedlingSurvival, SPcode=="Ccon"&MeshSize=="L")
CCO.S.SeedlingSurvival<-subset(SeedlingSurvival, SPcode=="Ccon"&MeshSize=="S")
prop.test(x=c(sum(CCO.L.SeedlingSurvival$survival),sum(CCO.S.SeedlingSurvival$survival)),n=c(48,48),alternative = "less")

OGL.L.SeedlingSurvival<-subset(SeedlingSurvival, SPcode=="Ogla"&MeshSize=="L")
OGL.S.SeedlingSurvival<-subset(SeedlingSurvival, SPcode=="Ogla"&MeshSize=="S")
prop.test(x=c(sum(OGL.L.SeedlingSurvival$survival),sum(OGL.S.SeedlingSurvival$survival)),n=c(48,48),alternative = "less")

SSU.L.SeedlingSurvival<-subset(SeedlingSurvival, SPcode=="Ssup"&MeshSize=="L")
SSU.S.SeedlingSurvival<-subset(SeedlingSurvival, SPcode=="Ssup"&MeshSize=="S")
prop.test(x=c(sum(SSU.L.SeedlingSurvival$survival),sum(SSU.S.SeedlingSurvival$survival)),n=c(48,48),alternative = "less")

ECM.L.SeedlingSurvival<-subset(SeedlingSurvival, MycorrhizalType=="ECM"&mesh=="L")
ECM.S.SeedlingSurvival<-subset(SeedlingSurvival, MycorrhizalType=="ECM"&mesh=="S")
prop.test(x=c(sum(ECM.L.SeedlingSurvival$survival),sum(ECM.S.SeedlingSurvival$survival)),n=c(192,192),alternative = "greater")

AM.L.SeedlingSurvival<-subset(SeedlingSurvival, MycorrhizalType=="AM"&mesh=="L")
AM.S.SeedlingSurvival<-subset(SeedlingSurvival, MycorrhizalType=="AM"&mesh=="S")
prop.test(x=c(sum(AM.L.SeedlingSurvival$survival),sum(AM.S.SeedlingSurvival$survival)),n=c(192,192),alternative = "less")

ECM.SeedlingSurvival$treatment<-rep(c("35 µm","0.5 µm"),each=4)
AM.SeedlingSurvival$treatment<-rep(c("35 µm","0.5 µm"),each=4)

library("ggplot2")

p.SeedlingSurvival.ecm<-ggplot(ECM.SeedlingSurvival,aes(x=SP,y=survival,fill=treatment))+
  geom_bar(stat="identity",width=0.7,position=position_dodge(0.8),colour="black")+
  geom_errorbar(aes(ymax=ECM.SeedlingSurvival$survival1+ECM.SeedlingSurvival$se,ymin=ECM.SeedlingSurvival$survival1-ECM.SeedlingSurvival$se),
                position=position_dodge(0.8),width=0.2)+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank())+
  xlab("ECM species") + 
  ylab("Seedling survival rate") + 
  theme(axis.text = element_text(size=rel(2)), axis.title = element_text(size=rel(2)), 
        strip.text = element_text(size=rel(2))) +
  guides(fill=guide_legend(title=NULL))+
  theme(legend.position=c(0.99,0.99), legend.justification=c(1,1),legend.direction="horizontal",legend.text = element_text(size = 24))+
  scale_fill_manual(values=c("white","grey"),guide=FALSE)+  
  annotate("text",x=0.6,y=1.05,label="(a)",size=10)+
  annotate("text",x=1,y=0.95,label="*",size=16)+
  annotate("text",x=2,y=1.00,label="*",size=16)+
  annotate("text",x=3,y=0.95,label="*",size=16);p.SeedlingSurvival.ecm

p.SeedlingSurvival.am<-ggplot(AM.SeedlingSurvival,aes(x=SP,y=survival,fill=treatment))+
  geom_bar(stat="identity",width=0.7,position=position_dodge(0.8),colour="black")+
  geom_errorbar(aes(ymax=AM.SeedlingSurvival$survival1+AM.SeedlingSurvival$se,ymin=AM.SeedlingSurvival$survival1-AM.SeedlingSurvival$se),
                position=position_dodge(0.8),width=0.2)+
  scale_fill_manual(values=c("white","grey"),guide=FALSE)+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank())+
  xlab("AM species") + 
  ylab("Seeding survival rate") +
  theme(axis.text = element_text(size=rel(2)), axis.title = element_text(size=rel(2)), 
        strip.text = element_text(size=rel(2))) +
  guides(fill=guide_legend(title=NULL))+
  theme(legend.position="none")+
  annotate("text",x=0.6,y=1.05,label="(b)",size=10)+
  annotate("text",x=1,y=0.93,label="*",size=16);p.SeedlingSurvival.am

#require(gridExtra)
#ggsave(grid.arrange(p.SeedlingSurvival.ecm,p.SeedlingSurvival.am,heights=c(1.7/4, 1.5/4),nrow=2), file="SeedlingSurvivaldensity.png", width=6, height=8)


####Combine ECM species or AM species together###########

SeedlingSurvivalTotal<-aggregate(survival~SPcode+Core+MeshSize,data=SeedlingSurvival,mean)

ECM.SeedlingSurvival.total<-rbind(subset(SeedlingSurvivalTotal,SP=="Cfab"),subset(SeedlingSurvivalTotal,SP=="Cfis"),subset(SeedlingSurvivalTotal,SP=="Chui"),subset(SeedlingSurvivalTotal,SP=="Lhai"))
ECM.SeedlingSurvival.total.L<-subset(ECM.SeedlingSurvival.total,mesh=="L")
ECM.SeedlingSurvival.total.S<-subset(ECM.SeedlingSurvival.total,mesh=="S")
ECM.L.mean = mean(ECM.SeedlingSurvival.total.L$survival);ECM.L.se = sd(ECM.SeedlingSurvival.total.L$survival)/sqrt(24)
ECM.S.mean = mean(ECM.SeedlingSurvival.total.S$survival);ECM.S.se = sd(ECM.SeedlingSurvival.total.S$survival)/sqrt(24)

total.ecm<-data.frame(treatment=c("35 µm", "0.5 µm"), mean=c(ECM.L.mean,ECM.S.mean), se = c(ECM.L.se, ECM.S.se))

p.total.ecm<-ggplot(total.ecm,aes(x=treatment,y=mean,fill=c("white","grey50"))) +
  geom_bar(position="dodge",stat="identity",width=0.8,colour="black") +
  geom_errorbar(aes(ymax=total.ecm$mean + total.ecm$se,ymin=total.ecm$mean - total.ecm$se),
                position=position_dodge(0.8),width=0.2)+
  scale_fill_manual(values=c("white","grey"),guide=FALSE)+
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank())+
  xlab("All ECM species") + 
  theme(axis.text = element_text(size=rel(2)), axis.title = element_text(size=rel(2)), 
        strip.text = element_text(size=rel(2))) +
  guides(fill=guide_legend(title=NULL)) +
  annotate("text",x=0.6,y=1,label="(c)",size=10) +
  annotate("text",x=1.5,y=0.93,label="*",size=16) +
  theme(legend.position="none") +
  ylab(" ") + 
  scale_y_continuous(breaks=NULL);p.total.ecm

AM.SeedlingSurvival.total<-rbind(subset(SeedlingSurvivalTotal,SP=="Calb"),subset(SeedlingSurvivalTotal,SP=="Ccon"),subset(SeedlingSurvivalTotal,SP=="Ogla"),subset(SeedlingSurvivalTotal,SP=="Ssup"))
AM.SeedlingSurvival.total.L<-subset(AM.SeedlingSurvival.total,mesh=="L")
AM.SeedlingSurvival.total.S<-subset(AM.SeedlingSurvival.total,mesh=="S")
AM.L.mean=mean(AM.SeedlingSurvival.total.L$survival1);AM.L.se=sd(AM.SeedlingSurvival.total.L$survival1)/sqrt(24)
AM.S.mean=mean(AM.SeedlingSurvival.total.S$survival1);AM.S.se=sd(AM.SeedlingSurvival.total.S$survival1)/sqrt(24)

total.am<-data.frame(treatment=c("35 µm", "0.5 µm"), mean=c(AM.L.mean,AM.S.mean), se = c(AM.L.se, AM.S.se))

p.total.am<-ggplot(total.am,aes(x=treatment,y=mean,fill=c("white","grey50"))) +
  geom_bar(position="dodge",stat="identity",width=0.8,colour="black") +
  geom_errorbar(aes(ymax=total.am$mean + total.am$se,ymin=total.am$mean - total.am$se),
                position=position_dodge(0.8),width=0.2) +
  scale_fill_manual(values=c("white","grey"),guide=FALSE) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank())+
  xlab("All AM species") + 
  ylab(" ") + 
  theme(axis.text = element_text(size=rel(2)), axis.title = element_text(size=rel(2)), 
        strip.text = element_text(size=rel(2))) +
  guides(fill=guide_legend(title=NULL))+
  annotate("text",x=0.6,y=1,label="(d)",size=10)+
  theme(legend.position="none")+
  scale_y_continuous(breaks=NULL);p.total.am

require(gridExtra)
grid.arrange(p.SeedlingSurvival.ecm,p.SeedlingSurvival.am,p.total.ecm,p.total.am,heights=c(1.7/4,1.5/4),
             layout_matrix = rbind(c(1,1,1,1,3),c(2,2,2,2,4)))

# PNG 1420 * 800, PDF 8.6*15.2

#############################################
### DNA sequence data ###

### Plot the pathogen / mycorrhizal ratio

PMratio<-read.csv("PathogenMycorrhizalRatio.csv", header = T, row.names = 1)

ecm<-subset(PMratio,mtype=="ecm")
am<-subset(PMratio,mtype=="am")

## PMratio

t.test(pm~mtype,data = PMratio)

se <- function(x) {sd(x)/sqrt(length(x))}
PMratiomean<-data.frame(mtype=c("ECM species","AM species"),
                        mean=c(mean(ecm$pm),mean(am$pm)),
                        se=c(se(ecm$pm),se(am$pm)))

p.PMratio<-ggplot(data=PMratiomean,aes(x=mtype,y=mean))+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),data=PMratiomean,
                position = "dodge", width = 0.1)+
  geom_point(size=4,shape=21,fill="white")+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank())+
  ylim(0,3)+
  xlab(" ") + theme(axis.title.x = element_text(size = 15, face = "bold"))+
  ylab("Pathogen / Mycorrhizas ratio") + 
  scale_x_discrete(limits=c("ECM species","AM species")) +
  annotate("text",x=1.5,y=1.2,label="*",size=10) +
  annotate("text",x=0.6,y=3,label="(b)",size=6)+
  theme(axis.text = element_text(size=rel(1.2)), axis.title = element_text(size=rel(1.2)), 
        strip.text = element_text(size=rel(1.2)));p.PMratio

# PNG 460 * 440, pdf 4.4*4.6

## Abundances of pathogenic fungi

t.test(Pathogen.Absolute~mtype,data = PMratio)

t.test(Pathogen.Relative~mtype,data = PMratio)

se <- function(x) {sd(x)/sqrt(length(x))}
Pathogen.Absolute.mean<-data.frame(mtype=c("ECM species","AM species"),
                        mean=c(mean(ecm$Pathogen.Absolute),mean(am$Pathogen.Absolute)),
                        se=c(se(ecm$Pathogen.Absolute),se(am$Pathogen.Absolute)))

p.Pathogen.Absolute<-ggplot(data=Pathogen.Absolute.mean,aes(x=mtype,y=mean))+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),data=Pathogen.Absolute.mean,
                position = "dodge", width = 0.1)+
  geom_point(size=4,shape=21,fill="white")+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ylim(0,3000)+
  xlab(" ") + theme(axis.title.x = element_text(size = 15, face = "bold"))+
  ylab("Pathogen read abundance") + 
  scale_x_discrete(limits=c("ECM species","AM species")) +
  annotate("text",x=0.6,y=3000,label="P < 0.05",size=6) +
  #annotate("text",x=0.6,y=3000,label="(a)",size=6)+
  theme(axis.text = element_text(size=rel(1.2)), axis.title = element_text(size=rel(1.2)), 
        strip.text = element_text(size=rel(1.2)));p.Pathogen.Absolute
