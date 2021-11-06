##################################################################################################################
##### R code belonging to Bongers et al 2021 Nature Ecology and Evolotion
##### Functional diversity effects on productivity increase with age in a forest biodiveristy experiment
#################################################################################################################

rm(list=ls())
{
  library(ggplot2)
  library(colorspace)## color panels for ggplot
  library(plyr)#count()
  library(dplyr)
  library(reshape2)#dcast() / melt() data
  library(gridExtra)
  library(lme4) #lmer
  library(lmerTest)#lmer with pvalue
  library(MuMIn)#r.squaredGLMM
  library(readxl)## import sheets of excel files
}


#########################################################################################
######## ------------------------LOAD THE DATA ---------------------------------#########
setwd()
setwd("c:/Users/Thinkpad/Documents/IBCAS/P2 Community and traits/20200515 Community trait_FINAL Final/manuscript/Submission NEE/Revision/Manuscript R1/Revision R2/")

## trait descriptions
df.traits<-as.data.frame(read_excel("Bongers et al NEE_species trait values.xlsx",sheet="Table S1"))
head(df.traits)
df.traits[,c("TRAITS","TRAITS.new","TRAITgroup")]

###meta data
meta.data<-read.csv("Meta data_original.csv")
str(meta.data)
table(meta.data[meta.data$MODEL=="All",c("AGE","RESPONSE","CAT")])
#Estimate --> slope of the relationships
# Pr...t.. --> P value of the relationships
# R2 --> R2 of the relationship
# RESPONSE -> is the response variable the cumulative stand-volume (VOL) or the yearly stand-volume increment (INCR)
# MODEL -> grouping of the collective dataset and the 10 subsets
# CAT -> the four different community indices = CWM, FD
# TRAIT -> contain the original 38 trait names
# AGE -> numeric value of the year of the analyses, ranging from 1-10 for VOL and 1-9 for RESPONSE = INCR

########################################################################################
### list with the traits that need to be corrected with -1 to get a general positive effect
t.cor<-c("LCC","Ca","CN","LDMC","LEAFTHICK","LEAFTOUGH","LOG10RATIO",
         "MEANROUND","P","PALIS","PSI50","SPONGY","STOMDENS","UPPEREPI",
         "VEINDENS","VPDMAX","VPDPOI")
#### correct the estimate values with -1
meta.data$Estimate.corr<-ifelse(meta.data$TRAIT%in%t.cor&meta.data$CAT=="CWM", meta.data$Estimate*-1,meta.data$Estimate)

# TRAIT.new -> the trait names used in publication, 
TRAITS<-as.vector(df.traits$TRAITS)
meta.data$TRAIT.order <- factor(meta.data$TRAIT, levels = TRAITS)
TRAITS.new<-as.vector(df.traits$TRAITS.new)
meta.data$TRAIT.new<-df.traits$TRAITS.new[match(meta.data$TRAIT,df.traits$TRAITS)]
meta.data$TRAIT.new <- factor(meta.data$TRAIT.new, levels = TRAITS.new)
names(meta.data)

## create trait group in the right order
TRAITgr.group<-c("categorical (2)","LES (6)","leaf structure (6)","leaf chemicals (6)","conductance (9)","hydraulics (9)")
meta.data$TRAITgr<-df.traits$TRAITgroup[match(meta.data$TRAIT,df.traits$TRAITS)]
meta.data$TRAITgr <- factor(meta.data$TRAITgr, levels = TRAITgr.group )
table(meta.data[,c("TRAIT.new","TRAITgr")])

#####################################################################################################################
#### ---------------- FIGURE 2 --> effect sizes (SLOPES)
##################################################################################################################
## test if the effect is different among CWM vs FD and change with AGE
## USE THE CORRECTED ESTIMATES TO HAVE AN OVERALL POSITIVE EFFECT

### subset data for analyses
dat<-droplevels(subset(meta.data, MODEL=="All"&RESPONSE=="VOL.a")) ## accumulated stand volume
dat<-droplevels(subset(meta.data, MODEL=="All"&RESPONSE=="VOLgr")) ## annual stand volume increment

## center AGE 
dat$AGE2<-dat$AGE-mean(dat$AGE)
dat$AGE2.sq<-dat$AGE2^2

### mixed-effects model
m2<-lmer(Estimate.corr~CAT*(AGE2+AGE2.sq)
         +(1+(AGE2+AGE2.sq)|CAT:TRAIT),data=dat)
anova(m2,type=1)
summary(m2)
r.squaredGLMM(m2)
dat$predict.Est.corr.m2<-predict(m2)

# calculate the mean value to include in the figure as point and get values for in the text
m<-aggregate(dat[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat$AGE,CAT=dat$CAT),FUN=mean)
sd<-aggregate(dat[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat$AGE,CAT=dat$CAT),FUN=sd)
n<-aggregate(dat[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat$AGE,CAT=dat$CAT),FUN=length)
sd$Estimate.corr.se<-sd$Estimate.corr/sqrt(n$Estimate.corr)
msd<-merge(m,sd[,c("AGE","CAT","Estimate.corr.se")],by=c("AGE","CAT"))
msd$CAT2<-""### create an empty grey bar above the plot ot make the plot same sizes

###### paired-t-test for CAT difference per year
test<-dcast(dat[dat$AGE==1,c("CAT","TRAIT","Estimate.corr")], TRAIT~CAT,fun=mean)
t.test(test$CWM,test$FD,paired=T)

## standing volume
pEff<-ggplot(
)+geom_point(data=msd, aes(x=AGE, y=Estimate.corr,shape=CAT,col=CAT),size=2,position=position_dodge(width=0.5)
)+geom_errorbar(data=msd, aes(x=AGE, ymin=Estimate.corr-Estimate.corr.se,ymax=Estimate.corr+Estimate.corr.se,col=CAT),position=position_dodge(width=0.5)
)+scale_color_manual(values=diverge_hcl(2,h=c(250,40),c=96)
)+theme_bw(
)+scale_y_continuous(limits = c(-11,60), breaks=c(0,20,40,60),labels=c(0,20,40,60)
)+scale_x_continuous(limits = c(0.5,10.5), breaks=c(1:10),labels=c(1:10)
)+theme(panel.grid = element_blank(),legend.position = c(0.15,0.9),legend.title = element_blank()
)+labs(x="Stand age (year)",y="Trait effect on stand volume (slope)",title="a"
)+theme (plot.title = element_text (hjust = -0.1, vjust=-5)
)+geom_line(data = msd, aes(x=AGE,y=predict.Est.corr.m2,col=CAT)
)+annotate(geom="text",x=1,y=-11,label="*")+annotate(geom="text",x=2,y=-11,label="*")+annotate(geom="text",x=3,y=-11,label="*"
)+annotate(geom="text",x=4,y=-11,label="*")+annotate(geom="text",x=5,y=-11,label="*")+annotate(geom="text",x=6,y=-11,label="*"
)+facet_grid(.~CAT2)
pEff

## increment
pEff<-ggplot(
)+geom_point(data=msd, aes(x=AGE, y=Estimate.corr,shape=CAT,col=CAT),size=2,position=position_dodge(width=0.5)
)+geom_errorbar(data=msd, aes(x=AGE, ymin=Estimate.corr-Estimate.corr.se,ymax=Estimate.corr+Estimate.corr.se,col=CAT),position=position_dodge(width=0.5)
)+scale_color_manual(values=diverge_hcl(2,h=c(250,40),c=96)
)+theme_bw(
)+scale_y_continuous(limits = c(-15,30), breaks=c(-10,0,10,20,30),labels=c(-10,0,10,20,30)
)+scale_x_continuous(limits = c(0.5,9.5), breaks=c(1:9),labels=c(1:9)
)+theme(panel.grid = element_blank(),legend.position = c(0.15,0.9),legend.title = element_blank()
)+labs(x="Stand age (year)",y="Trait effect on stand volume increment (slope)",title="a"
)+theme (plot.title = element_text (hjust = -0.1, vjust=-5)
)+geom_line(data = msd, aes(x=AGE,y=predict.Est.corr.m2,col=CAT)
)+annotate(geom="text",x=1,y=-15,label="*")+annotate(geom="text",x=2,y=-15,label="*")+annotate(geom="text",x=3,y=-15,label="*"
)+annotate(geom="text",x=4,y=-15,label="*")+annotate(geom="text",x=6,y=-15,label="*")+annotate(geom="text",x=7,y=-15,label="*"
)+annotate(geom="text",x=8,y=-15,label="*"
)+facet_grid(.~CAT2)
pEff

### trait group effects per CAT
cwm2<-lmer(Estimate.corr~TRAITgr*(AGE2+AGE2.sq)
           +(1+(AGE2+AGE2.sq)|TRAIT),data=subset(dat,CAT=="CWM"))
anova(cwm2,type=1)
r.squaredGLMM(cwm2)
AIC(cwm2)
fd2<-lmer(Estimate.corr~TRAITgr*(AGE2+AGE2.sq)
          +(1+(AGE2+AGE2.sq)|TRAIT),data=subset(dat,CAT=="FD"))
anova(fd2,type=1)
r.squaredGLMM(fd2)
AIC(fd2)
##
dat$predict.Est.corr.m2<-ifelse(dat$CAT=="CWM",predict(cwm2),predict(fd2))

## calculate the average for the estimates from the model
m<-aggregate(dat[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat$AGE,CAT=dat$CAT,TRAITgr=dat$TRAITgr),FUN=mean)
sd<-aggregate(dat[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat$AGE,CAT=dat$CAT,TRAITgr=dat$TRAITgr),FUN=sd)
n<-aggregate(dat[,c("Estimate.corr","predict.Est.corr.m2")],by=list(AGE=dat$AGE,CAT=dat$CAT,TRAITgr=dat$TRAITgr),FUN=length)
sd$Estimate.corr.se<-sd$Estimate.corr/sqrt(n$Estimate.corr)
msd2<-merge(m,sd[,c("AGE","CAT","TRAITgr","Estimate.corr.se")],by=c("AGE","CAT","TRAITgr"))

#### colour coding for the trait-groups
TRAITgr.col<-c("forestgreen","darksalmon","gold3","yellowgreen","steelblue4","firebrick")

### standing volume
pTraitGr<-ggplot(data=msd2, aes(x=AGE, y=Estimate.corr, col=TRAITgr,shape=TRAITgr)
)+geom_point(position=position_dodge(width=0.5)
)+geom_errorbar(aes(x=AGE,ymin=Estimate.corr-Estimate.corr.se,ymax=Estimate.corr+Estimate.corr.se),position=position_dodge(width=0.5)
)+geom_line(aes(x=(AGE),y=predict.Est.corr.m2,col=TRAITgr)
)+facet_grid(.~CAT
)+scale_color_manual(values=TRAITgr.col
)+theme_bw(
)+theme(panel.grid = element_blank(),legend.position = c(0.2,0.85),legend.title = element_blank(), legend.background = element_blank()
)+scale_y_continuous(limits = c(-11,100), breaks=c(0,25,50,75,100),labels=c(0,25,50,75,100)
)+scale_x_continuous(limits = c(0.5,10.5), breaks=c(1:10),labels=c(1:10)
)+labs(x="Stand age (year)",y="Trait effect on stand volume (slope)",title="b"
)+theme (plot.title = element_text (hjust = -0.12, vjust=-5))

#tiff("FIG 2 CWM and FD Effect across time and groups.tiff",unit="px",res=600,width=5000,height=3500)
pdf("FIG 2 CWM and FD Effect across time and groups.pdf",width=9,height=6)
grid.arrange(pEff,pTraitGr,nrow=1,ncol=2)
dev.off()

### INCREMENT
pTraitGr<-ggplot(data=msd2, aes(x=AGE, y=Estimate.corr, col=TRAITgr,shape=TRAITgr)
)+geom_point(position=position_dodge(width=0.5)
)+geom_errorbar(aes(x=AGE,ymin=Estimate.corr-Estimate.corr.se,ymax=Estimate.corr+Estimate.corr.se),position=position_dodge(width=0.5)
)+geom_line(aes(x=(AGE),y=predict.Est.corr.m2,col=TRAITgr)
)+facet_grid(.~CAT
)+scale_color_manual(values=TRAITgr.col
)+theme_bw(
)+theme(panel.grid = element_blank(),legend.position = c(0.2,0.85),legend.title = element_blank(), legend.background = element_blank()
)+scale_y_continuous(limits = c(-15,30), breaks=c(-10,0,10,20,30),labels=c(-10,0,10,20,30)
)+scale_x_continuous(limits = c(0.5,10.5), breaks=c(1:10),labels=c(1:10)
)+labs(x="Stand age (year)",y="Trait effect on stand volume increment (slope)",title="b"
)+theme (plot.title = element_text (hjust = -0.12, vjust=-5))
pTraitGr

tiff("SUPP CWM and FD effect across time and groups_INCREMENT.tiff",unit="px",res=600,width=5000,height=3500)
#pdf("SUPP CWM and FD Effect across time and groups_INCREMENT.pdf",width=9,height=6)
grid.arrange(pEff,pTraitGr,nrow=1,ncol=2)
dev.off()


#################### ------------------ SUPPLEMENTAL DATA ------------------------#########################
### subset data for analyses
dat<-droplevels(subset(meta.data, MODEL=="All"&RESPONSE=="VOL.a")) ## accumulated stand volume
## center AGE
dat$AGE2<-dat$AGE-mean(dat$AGE)
dat$AGE2.sq<-dat$AGE2^2

# mixed-effects model
m2<-lmer(Estimate.corr~CAT*(AGE2+AGE2.sq)
         +(1+(AGE2+AGE2.sq)|CAT:TRAIT),data=dat)
dat$predict.Est.corr.m2<-predict(m2)

## Figure with all trait specific data
pSI<-ggplot(data = dat
)+geom_point(aes(x=AGE,y=Estimate.corr,group=CAT,col=CAT,shape=CAT)
)+geom_line(aes(x=AGE,y=predict.Est.corr.m2,group=CAT,col=CAT)
)+facet_wrap(facets = ~ TRAIT.order ,nrow = 8, ncol=5
)+scale_x_continuous(breaks=c(1:10),labels=1:10
)+scale_color_manual(values=diverge_hcl(2,h=c(250,40),c=96)
)+theme_bw()+theme(legend.position = "none"
)+labs(y=expression("Fitted and observed trait effect on stand volume (slope)"),x="Stand age (year)"
)+theme(panel.grid = element_blank(),legend.position = c(0.70,0.05),legend.title = element_blank()
)
tiff("SUPP Fitted and observed effects across age per trait.tiff",res=600, unit="px",width=6000,height=8000)
#pdf("SUPP Fitted and observed Effect across age per trait.pdf",width=7,height=10)
print(pSI)
dev.off()


#############################################################################################
#### -------- FIGURE 3 ---> R2 across time, for CAT and TRAIT
###############################################################################################
dat<-droplevels(subset(meta.data, MODEL=="All"&RESPONSE=="VOL.a")) ## accumulated stand volume
dat<-droplevels(subset(meta.data, MODEL=="All"&RESPONSE=="VOLgr")) ## annual stand volume increment

### center AGE
dat$AGE2<-dat$AGE-mean(dat$AGE)
dat$AGE2.sq<-dat$AGE2^2

### mixed effects model
m2<-lmer(asin(sqrt(R2))~CAT*(AGE2+AGE2.sq)
         +(1+(AGE2+AGE2.sq)|CAT:TRAIT),data=dat)
anova(m2,type=1)
summary(m2)
r.squaredGLMM(m2)
AIC(m2)
dat$predict.R2.m2<-sin(predict(m2))^2
##### SEE --> backtransform asin(sqrt())
#   tr.asin<-asin(sqrt(0.1))
#   sin(tr.asin)^2

## calculate the mean value to include in the figure as point and get values for in the text
m<-aggregate(dat[,c("R2","predict.R2.m2")],by=list(AGE=dat$AGE,CAT=dat$CAT),FUN=mean)
sd<-aggregate(dat[,c("R2","predict.R2.m2")],by=list(AGE=dat$AGE,CAT=dat$CAT),FUN=sd)
n<-aggregate(dat[,c("R2","predict.R2.m2")],by=list(AGE=dat$AGE,CAT=dat$CAT),FUN=length)
sd$R2.se<-sd$R2/sqrt(n$R2)
msd<-merge(m,sd[,c("AGE","CAT","R2.se")],by=c("AGE","CAT"))
msd$CAT2<-""

###### paired-t-test for CAT difference per year
dat$R2.asinsq<-asin(sqrt(dat$R2))
### using the asin-squared transofrmation also for the paired t-test
test<-dcast(dat[dat$AGE==1,c("CAT","TRAIT","R2.asinsq")], TRAIT~CAT,fun=mean)
t.test(test$CWM,test$FD,paired=T)

### stand volume
pR2<-ggplot(
)+geom_point(data=msd, aes(x=(AGE), y=R2,col=CAT,shape=CAT),size=2,position=position_dodge(width=0.5)
)+geom_errorbar(data=msd, aes(x=(AGE), ymin=R2-R2.se,ymax=R2+R2.se,col=CAT),position=position_dodge(width=0.5)
)+scale_color_manual(values=diverge_hcl(2,h=c(250,40),c=96)
)+facet_grid(.~CAT2
)+geom_line(data = msd, aes(x=AGE,y=predict.R2.m2,col=CAT)
)+theme_bw(
)+scale_y_continuous(limits = c(-0.005,0.12)#, breaks=c(0,0.05,0.1,0.15),labels=c(0,0.05,0.1,0.15)
)+scale_x_continuous(limits = c(0.5,10.5), breaks=c(1:10),labels=c(1:10)
)+theme(panel.grid = element_blank(),legend.position = c(0.15,0.9),legend.title = element_blank()
)+labs(y=expression("Reliability to predict stand volume ("~R^2~")"),x="Stand age (year)",title="a"
)+theme (plot.title = element_text (hjust = -0.15, vjust=-5)
)+annotate(geom="text",x=1,y=-0.005,label="*")+annotate(geom="text",x=2,y=-0.005,label="*")+annotate(geom="text",x=3,y=-0.005,label="*"
)+annotate(geom="text",x=4,y=-0.005,label="*")+annotate(geom="text",x=5,y=-0.005,label="*")+annotate(geom="text",x=8,y=-0.005,label="*"
)+annotate(geom="text",x=9,y=-0.005,label="*")+annotate(geom="text",x=10,y=-0.005,label="*"
)
pR2

#### INCREMENT
pR2<-ggplot(
)+geom_point(data=msd, aes(x=(AGE), y=R2,col=CAT,shape=CAT),size=2,position=position_dodge(width=0.5)
)+geom_errorbar(data=msd, aes(x=(AGE), ymin=R2-R2.se,ymax=R2+R2.se,col=CAT),position=position_dodge(width=0.5)
)+scale_color_manual(values=diverge_hcl(2,h=c(250,40),c=96)
)+facet_grid(.~CAT2
)+geom_line(data = msd, aes(x=AGE,y=predict.R2.m2,col=CAT)
)+theme_bw(
)+scale_y_continuous(limits = c(-0.005,0.12)#, breaks=c(0,0.05,0.1,0.15),labels=c(0,0.05,0.1,0.15)
)+scale_x_continuous(limits = c(0.5,9.5), breaks=c(1:9),labels=c(1:9)
)+theme(panel.grid = element_blank(),legend.position = c(0.15,0.9),legend.title = element_blank()
)+labs(y=expression("Reliability to predict stand volume increment ("~R^2~")"),x="Stand age (year)",title="c"
)+theme (plot.title = element_text (hjust = -0.15, vjust=-5)
)+annotate(geom="text",x=1,y=-0.005,label="*")+annotate(geom="text",x=2,y=-0.005,label="*")+annotate(geom="text",x=3,y=-0.005,label="*"
)+annotate(geom="text",x=4,y=-0.005,label="")+annotate(geom="text",x=5,y=-0.005,label="*")+annotate(geom="text",x=6,y=-0.005,label="*"
)+annotate(geom="text",x=7,y=-0.005,label="*")+annotate(geom="text",x=8,y=-0.005,label="*"
)
pR2

#### test effect of TRAITgroup with AGE per category
### mixed-effects model 
cwm2<-lmer(asin(sqrt(R2))~TRAITgr*(AGE2+AGE2.sq)
           +(1+(AGE2+AGE2.sq)|TRAIT),data=subset(dat,CAT=="CWM"))
anova(cwm2,type=1)
r.squaredGLMM(cwm2)
AIC(cwm2)
fd2<-lmer(asin(sqrt(R2))~TRAITgr*(AGE2+AGE2.sq)
          +(1+(AGE2+AGE2.sq)|TRAIT),data=subset(dat,CAT=="FD"))
anova(fd2,type=1)
r.squaredGLMM(fd2)
AIC(fd2)
##
dat$predict.R2.group.m2<-ifelse(dat$CAT=="CWM",sin(predict(cwm2))^2,sin(predict(fd2))^2)

## calculate the average for the R2 from the model
m<-aggregate(dat[,c("R2","predict.R2.group.m2")],by=list(AGE=dat$AGE,CAT=dat$CAT,TRAITgr=dat$TRAITgr),FUN=mean)
sd<-aggregate(dat[,c("R2","predict.R2.group.m2")],by=list(AGE=dat$AGE,CAT=dat$CAT,TRAITgr=dat$TRAITgr),FUN=sd)
n<-aggregate(dat[,c("R2","predict.R2.group.m2")],by=list(AGE=dat$AGE,CAT=dat$CAT,TRAITgr=dat$TRAITgr),FUN=length)
sd$R2.se<-sd$R2/sqrt(n$R2)
msd2<-merge(m,sd[,c("AGE","CAT","TRAITgr","R2.se")],by=c("AGE","CAT","TRAITgr"))
msd2$meanse<-msd2$R2-msd2$R2.se

#### colour coding for the trait-groups
TRAITgr.col<-c("forestgreen","darksalmon","gold3","yellowgreen","steelblue4","firebrick")

## stand volume
pR2.group<-ggplot(data=msd2, aes(x=AGE, y=R2, col=TRAITgr,shape=TRAITgr)
)+geom_point(position=position_dodge(width=0.5)
)+geom_errorbar(aes(x=AGE,ymin=R2-R2.se,ymax=R2+R2.se),position=position_dodge(width=0.5)
)+facet_grid(.~CAT
)+geom_line(aes(x=(AGE),y=predict.R2.group.m2,col=TRAITgr)
)+scale_color_manual(values=TRAITgr.col
)+theme_bw(
)+theme(panel.grid = element_blank(),legend.position = c(0.7,0.85),legend.title = element_blank(), legend.background = element_blank()
)+scale_y_continuous(limits = c(-0.005,0.12)#, breaks=c(0,0.05,0.1,0.15),labels=c(0,0.05,0.1,0.15)
)+scale_x_continuous(limits = c(0.5,10.5), breaks=c(1:10),labels=c(1:10)
)+labs(y=expression("Reliability to predict stand volume ("~R^2~")"),x="Stand age (year)",title="b"
)+theme (plot.title = element_text (hjust = -0.15, vjust=-5)
)
pR2.group

tiff("FIG 3 R2 across time for CWM and FD and Groups.tiff",unit="px",res=600,width=5000,height=3500)
pdf("FIG 3 across time for CWM and FD and Groups.pdf",width=9,height=6)
grid.arrange(pR2,pR2.group,ncol=2,nrow=1)
dev.off()

### INCREMENT
pR2.group<-ggplot(data=msd2, aes(x=AGE, y=R2, col=TRAITgr,shape=TRAITgr)
)+geom_point(position=position_dodge(width=0.5)
)+geom_errorbar(aes(x=AGE,ymin=R2-R2.se,ymax=R2+R2.se),position=position_dodge(width=0.5)
)+facet_grid(.~CAT
)+geom_line(aes(x=(AGE),y=predict.R2.group.m2,col=TRAITgr)
)+scale_color_manual(values=TRAITgr.col
)+theme_bw(
)+theme(panel.grid = element_blank(),legend.position = c(0.7,0.85),legend.title = element_blank(), legend.background = element_blank()
)+scale_y_continuous(limits = c(-0.005,0.12)#, breaks=c(0,0.05,0.1,0.15),labels=c(0,0.05,0.1,0.15)
)+scale_x_continuous(limits = c(0.5,10.5), breaks=c(1:10),labels=c(1:10)
)+labs(y=expression("Reliability to predict stand volume increment ("~R^2~")"),x="Stand age (year)",title="d"
)+theme (plot.title = element_text (hjust = -0.15, vjust=-5)
)
pR2.group

tiff("SUPP Effect and R2 across time for CWM and FD and Groups_INCREMENT.tiff",unit="px",res=600,width=5500,height=7000)
pdf("SUPP Effect and R2 across time for CWM and FD and Groups_INCREMENT.pdf",width=9,height=12)
grid.arrange(pEff,pTraitGr,pR2,pR2.group,ncol=2,nrow=2)
dev.off()

###########################-------- SUPPLEMeNtAL FIGURE ----------##############################
dat<-droplevels(subset(meta.data, MODEL=="All"&RESPONSE=="VOL.a")) ## accumulated stand volume

### center AGE
dat$AGE2<-dat$AGE-mean(dat$AGE)
dat$AGE2.sq<-dat$AGE2^2

### mixed effects model
m2<-lmer(asin(sqrt(R2))~CAT*(AGE2+AGE2.sq)
         +(1+(AGE2+AGE2.sq)|CAT:TRAIT),data=dat)
dat$predict.R2.m2<-sin(predict(m2))^2

# figure that includes the observed and predicted values from the model, ordered per trait in separate panel
pSI<-ggplot(data = dat)+geom_point(aes(x=AGE,y=R2,group=CAT,col=CAT,shape=CAT)
)+geom_line(aes(x=AGE,y=predict.R2.m2,group=CAT,col=CAT)
)+facet_wrap(facets = ~ TRAIT.order ,nrow = 8, ncol=5
)+scale_x_continuous(breaks=c(1:10),labels=1:10
)+scale_color_manual(values=diverge_hcl(2,h=c(250,40),c=96)
)+theme_bw()+theme(legend.position = "none"
)+labs(y=expression("Fitted and observed reliability to predict stand volume ("~R^2~")"),x="Stand age (year)"
)+theme(panel.grid = element_blank(),legend.position = c(0.70,0.05),legend.title = element_blank()
)
tiff("SUPP Fitted and observed R2 per trait.tiff",res=600, unit="px",width=6000,height=8000)
#pdf("SUPP Fitted and observed R2 per trait.pdf",width=7,height=10)
print(pSI)
dev.off()


####################################################################################################
#######----- FIGURE 4 --> Effect sizes for 10 pools -----------------------------------------########
#####################################################################################################
names(meta.data)
plyr::count(meta.data$MODEL)

######################################################################################
###----- calculate how many species pools have negative or positive effect at AGE==10 or 9 for increment
POOLS<-c("A.A","A.B","A.C","B.A","B.B","B.C","A.SLA","A.Rar","B.SLA","B.Rar")

## set the correct dataset
dat<-droplevels(subset(meta.data, RESPONSE=="VOL.a"&MODEL%in%POOLS))

dat$Est.f<-ifelse(dat$Estimate.corr<0,"neg","pos")
dat$sig<-dat$Pr...t..
dat$Sig.f<-ifelse(dat$sig<0.05,"P<0.05","ns")

# TRAIT.new -> the trait names used in publication, 
names(df.traits)
head(df.traits)

dat$R2.f<-ifelse(dat$R2<0.15,"alow","bhigh")

TRAITgr.col<-c("forestgreen","darksalmon","gold3","yellowgreen","steelblue4","firebrick")
TRAITgr.38col<-c('forestgreen','forestgreen',
                      'darksalmon','darksalmon','darksalmon','darksalmon','darksalmon','darksalmon',
                      'gold3','gold3','gold3','gold3','gold3','gold3',
                      'yellowgreen','yellowgreen','yellowgreen','yellowgreen','yellowgreen','yellowgreen',
                      'steelblue4','steelblue4','steelblue4','steelblue4','steelblue4','steelblue4','steelblue4','steelblue4','steelblue4',
                      'firebrick','firebrick','firebrick','firebrick','firebrick','firebrick','firebrick','firebrick','firebrick')

dat2<-droplevels(subset(dat,AGE==10))
table(dat2[,c("TRAIT.new","Est.f","CAT")])

p1<-ggplot(data =dat2, aes(y=reorder(TRAIT.new, desc(TRAIT.new)), x=Estimate.corr, col=Est.f,shape=R2.f) 
)+geom_point()+facet_grid(.~CAT,scales = "free"
)+scale_color_manual(values= c("chocolate3","aquamarine4")
)+scale_shape_manual(values= c(1,16)
)+theme_bw()+theme(legend.position = "none",panel.grid = element_blank()
)+labs(x="Trait effect on stand volume (slope)",y="")+geom_vline(xintercept = 0,lty=2
)+theme(axis.text.y = element_text(color = rev(TRAITgr.38col)),strip.placement = "outside"
)
tiff("FIG 4 Species pool effects.tiff", unit="px",res=600,width=4000,height=3200)
pdf("FIG 4 Species pool effects.pdf",width=7.5,height=6)
print(p1)
dev.off()
