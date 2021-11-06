##################################################################################################################
##### R code belonging to Bongers et al 2021 Nature Ecology and Evolotion
##### Functional diversity effects on productivity increase with age in a forest biodiveristy experiment
#################################################################################################################

rm(list=ls())
{
library(ggplot2)
library(colorspace)# color panels 
library(plyr)#count # load this package first to keep functions well
library(dplyr)
library(reshape2)#dcast() / melt() data
library(psych)#principal() = PCA / cor.plot()
library(gridExtra)# make figures with panels
library(readxl)## import sheets of excel files
library(FD) ### calculate CWM and FDis values
library(MuMIn)# model selection but also r.squaredGLM()
require(vegan) # pca analyases - rda()
}

########################################################################################
#########--------------------------- LOAD DATA---------------------------------#########
########################################################################################
setwd()
#### original stand volume data
df.478<-as.data.frame(read_excel("Data_Bongers et al NEE_stand volume and increment.xlsx",sheet="Raw data", na="NA"))

##### designed relative species proportions per PTAG
df.design<-read.csv("Bongers et al NEE_species proportions.csv")
str(df.design)

#### trait data 
df.traits<-as.data.frame(read_excel("Bongers et al NEE_species trait values.xlsx",sheet="Species mean trait values",na=""))
head(df.traits)[1:10]
str(df.traits)

### SELECTED TRAITS
TRAITS<-c("D","AM", ## categorical traits
          "SLA","LOG10LA","LDMC","LNC","LCC","CN", # LES
          "LEAFTHICK","LEAFTOUGH","UPPEREPI", "PALIS","SPONGY", "LOG10RATIO", #structure
          "P","Phen","Tan", "Ca",  "K", "Mg",#chemicals
          "STOMDENS", "STOMSIZE", "STOMIND", "CONMEAN","CONMAX",
          "CONMAXFIT","VPDMAX",  "VPDMAXFIT",  "VPDPOI", # conductance traits
          "WD","HYDCOND","PSI50","WPOT","VEINDENS",
          "DIAMVEIN","MEANAREA","MEANROUND",  "DHYD") # hydraulic traits

TRAITS%in%names(df.traits)

####################################################################################################
####--------------------- PCA analysis to make PC axis per trait group------------------######
######################################################################################################
TRAITS
######## all 38 traits
##### run pca analysis
set<-df.traits[,TRAITS]
pca <- rda(set, choices=1:2,scale = TRUE) 
summary(pca)$cont

smry <- summary(pca)
df1  <- data.frame(smry$sites[,1:2])   # PC1 and PC2
df1$PC1<-df1$PC1*-1 ### flip the direction of the pc values beause original values have a negative effect on VOLUME
colnames(df1)<-c("PC1.38","PC2.38")
df.traits<-cbind(df.traits,df1)

df2  <- data.frame(smry$species[,1:2])     # loadings for PC1 and PC2
df2$PC1<-df2$PC1*-1
df.loading<-cbind(df2,"TRAIT"=rownames(df2),"TRAITgr"="All (38)")

rda.plot <- ggplot(df1, aes(x=PC1.les, y=PC2.les)) + 
  geom_point(size=3,col="black") +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme_bw()

rda.biplot <- rda.plot +
  geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2), 
               color="grey", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2, 
            aes(x=PC1*scalef,y=PC2*scalef,label=rownames(df2),
                hjust=0.5*(1-sign(PC1)),vjust=0.5*(1-sign(PC2))), 
            color="black", size=2)

#################### create PCA per trait.group
TRAITS.les<-TRAITS[3:8]#LES traits
##### run pca analysis
set<-df.traits[,TRAITS.les]
pca <- rda(set, choices=1:2,scale = TRUE) 
summary(pca)$cont

smry <- summary(pca)
df1  <- data.frame(smry$sites[,1:2])   # PC1 and PC2
df1$PC1<-df1$PC1*-1 ### flip the direction of the pc values beause original values have a negative effect on VOLUME
colnames(df1)<-c("PC1.les","PC2.les")
df.traits<-cbind(df.traits,df1)

df2  <- data.frame(smry$species[,1:2])     # loadings for PC1 and PC2
df2$PC1<-df2$PC1*-1
df.loading<-rbind(df.loading,cbind(df2,"TRAIT"=rownames(df2),"TRAITgr"="LES (6)"))

rda.plot <- ggplot(df1, aes(x=PC1.les, y=PC2.les)) + 
  geom_point(size=3,col="black") +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme_bw()

rda.biplot <- rda.plot +
  geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2), 
               color="grey", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2, 
            aes(x=PC1*scalef,y=PC2*scalef,label=rownames(df2),
                hjust=0.5*(1-sign(PC1)),vjust=0.5*(1-sign(PC2))), 
            color="black", size=2)

##### 
TRAITS.struc<-TRAITS[9:14]
##### run pca analysis
set<-df.traits[,TRAITS.struc]
pca <- rda(set, choices=1:2,scale = TRUE) 
summary(pca)$cont

smry <- summary(pca)
df1  <- data.frame(smry$sites[,1:2])   # PC1 and PC2
df1$PC1<-df1$PC1*-1 ### flip the direction of the pc values beause original values have a negative effect on VOLUME
colnames(df1)<-c("PC1.struc","PC2.struc")
df.traits<-cbind(df.traits,df1)

df2  <- data.frame(smry$species[,1:2])     # loadings for PC1 and PC2
df2$PC1<-df2$PC1*-1
df.loading<-rbind(df.loading,cbind(df2,"TRAIT"=rownames(df2),"TRAITgr"="leaf structure (6)"))

rda.plot <- ggplot(df1, aes(x=PC1.struc, y=PC2.struc)) + 
  geom_point(size=3,col="black") +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme_bw()

rda.biplot <- rda.plot +
  geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2), 
               color="grey", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2, 
            aes(x=PC1*scalef,y=PC2*scalef,label=rownames(df2),
                hjust=0.5*(1-sign(PC1)),vjust=0.5*(1-sign(PC2))), 
            color="black", size=2)
##### 
TRAITS.chem<-TRAITS[15:20]
##### run pca analysis
set<-df.traits[,TRAITS.chem]
pca <- rda(set, choices=1:2,scale = TRUE) 
summary(pca)$cont

smry <- summary(pca)
df1  <- data.frame(smry$sites[,1:2])   # PC1 and PC2
colnames(df1)<-c("PC1.chem","PC2.chem")
df.traits<-cbind(df.traits,df1)

df2  <- data.frame(smry$species[,1:2])     # loadings for PC1 and PC2
df.loading<-rbind(df.loading,cbind(df2,"TRAIT"=rownames(df2),"TRAITgr"="leaf chemicals (6)"))

rda.plot <- ggplot(df1, aes(x=PC1.chem, y=PC2.chem)) + 
  geom_point(size=3,col="black") +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme_bw()

rda.biplot <- rda.plot +
  geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2), 
               color="grey", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2, 
            aes(x=PC1*scalef,y=PC2*scalef,label=rownames(df2),
                hjust=0.5*(1-sign(PC1)),vjust=0.5*(1-sign(PC2))), 
            color="black", size=2)
##### 
TRAITS.conduc<-TRAITS[21:29]
##### run pca analysis
set<-df.traits[,TRAITS.conduc]
pca <- rda(set, choices=1:2,scale = TRUE) 
summary(pca)$cont

smry <- summary(pca)
df1  <- data.frame(smry$sites[,1:2])   # PC1 and PC2
colnames(df1)<-c("PC1.conduc","PC2.conduc")
df.traits<-cbind(df.traits,df1)

df2  <- data.frame(smry$species[,1:2])     # loadings for PC1 and PC2
df.loading<-rbind(df.loading,cbind(df2,"TRAIT"=rownames(df2),"TRAITgr"="conductance (9)"))

rda.plot <- ggplot(df1, aes(x=PC1.conduc, y=PC2.conduc)) + 
  geom_point(size=3,col="black") +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme_bw()

rda.biplot <- rda.plot +
  geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2), 
               color="grey", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2, 
            aes(x=PC1*scalef,y=PC2*scalef,label=rownames(df2),
                hjust=0.5*(1-sign(PC1)),vjust=0.5*(1-sign(PC2))), 
            color="black", size=2)

##### 
TRAITS.hydr<-TRAITS[30:38]
##### run pca analysis
set<-df.traits[,TRAITS.hydr]
pca <- rda(set, choices=1:2,scale = TRUE) 
summary(pca)$cont

smry <- summary(pca)
df1  <- data.frame(smry$sites[,1:2])   # PC1 and PC2
df1$PC1<-df1$PC1*-1 ### flip the direction of the pc values because original values have a negative effect on VOLUME
colnames(df1)<-c("PC1.hydr","PC2.hydr")
df.traits<-cbind(df.traits,df1)

df2  <- data.frame(smry$species[,1:2])     # loadings for PC1 and PC2
df2$PC1<-df2$PC1*-1
df.loading<-rbind(df.loading,cbind(df2,"TRAIT"=rownames(df2),"TRAITgr"="hydraulics (9)"))

rda.plot <- ggplot(df1, aes(x=PC1.hydr, y=PC2.hydr)) + 
  geom_point(size=3,col="black") +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme_bw()

rda.biplot <- rda.plot +
  geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2), 
               color="grey", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2, 
            aes(x=PC1*scalef,y=PC2*scalef,label=rownames(df2),
                hjust=0.5*(1-sign(PC1)),vjust=0.5*(1-sign(PC2))), 
            color="black", size=2)

##############################################################################################
############### ------- create the supplemental figure with the trait loadings --------############
df.loading
TRAITS%in%df.loading$TRAIT
#### make the order of the traits by setting them as factors with levels
df.loading$TRAIT.order <- factor(df.loading$TRAIT, levels = TRAITS)

TRAIT.group<-c("All (38)","LES (6)","leaf structure (6)","leaf chemicals (6)","conductance (9)","hydraulics (9)")
TRAIT.group%in%df.loading$TRAITgr
df.loading$TRAITgr.order <- factor(df.loading$TRAITgr, levels = TRAIT.group )

pp<-ggplot(data=df.loading,aes(x=PC1,y=reorder(TRAIT, desc(TRAIT.order)))
)+facet_grid(.~TRAITgr.order
)+geom_point(
)+geom_line(arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"), size = 3
)+theme_bw()+geom_vline(xintercept = 0,lty=2)+labs(y="",x="")

#tiff("20210625 loadings figure.tiff",width=3000,height=1500,unit="px",res=300)
print(pp)
dev.off()

#######################################################################################
## ------------- CALCULATE multivariate FD with multiple traits -----##############
### FIRST
str(df.traits)
rownames(df.traits)<-df.traits$species_name
SPECIES.traits<-df.traits$species_name
names(df.traits)
sp.tr.m<-as.matrix(df.traits[,3:length(df.traits)])
head(sp.tr.m)

######## SECOND
names(df.design)
abun.sp<-dcast(df.design[,c("PTAG","species_name","rel.abun.des")],PTAG~species_name)
head(abun.sp)
## save PTAg in rownames
rownames(abun.sp)<-abun.sp$PTAG
PTAG.list<-as.vector(abun.sp$PTAG)
PTAG.sp.m<-as.matrix(abun.sp[,SPECIES.traits])

###########################################################################################
#### calculate the multivariate trait value
### create dataset to save this data
FD.data<-as.data.frame(cbind(PTAG=PTAG.list))
head(FD.data)
{
## first use all 38 traits together for a FD using all 38 traits
sp.tr.m2<-as.matrix(sp.tr.m[,TRAITS])
fd<-dbFD(x=sp.tr.m2, a = PTAG.sp.m) ## note the dataset to include
a<-as.data.frame(fd$FDis)
colnames(a)<-"All38.FDis"
FD.data<-cbind(FD.data, a)
##### specific categorical
sp.tr.m2<-as.matrix(sp.tr.m[,TRAITS[1:2]])
fd<-dbFD(x=sp.tr.m2, a = PTAG.sp.m) ## note the dataset to include
a<-as.data.frame(fd$FDis)
colnames(a)<-"categorical.FDis"
FD.data<-cbind(FD.data, a)
##### specific LES selection
sp.tr.m2<-as.matrix(sp.tr.m[,TRAITS.les])
fd<-dbFD(x=sp.tr.m2, a = PTAG.sp.m) ## note the dataset to include
a<-as.data.frame(fd$FDis)
colnames(a)<-"LES.FDis"
FD.data<-cbind(FD.data, a)
####### structure
sp.tr.m2<-as.matrix(sp.tr.m[,TRAITS.struc])
fd<-dbFD(x=sp.tr.m2, a = PTAG.sp.m) ## note the dataset to include
a<-as.data.frame(fd$FDis)
colnames(a)<-"structure.FDis"
FD.data<-cbind(FD.data, a)
####### chemicals
sp.tr.m2<-as.matrix(sp.tr.m[,TRAITS.chem])
fd<-dbFD(x=sp.tr.m2, a = PTAG.sp.m) ## note the dataset to include
a<-as.data.frame(fd$FDis)
colnames(a)<-"chemicals.FDis"
FD.data<-cbind(FD.data, a)
####### conductance
sp.tr.m2<-as.matrix(sp.tr.m[,TRAITS.conduc])
fd<-dbFD(x=sp.tr.m2, a = PTAG.sp.m) ## note the dataset to include
a<-as.data.frame(fd$FDis)
colnames(a)<-"conductance.FDis"
FD.data<-cbind(FD.data, a)
####### hydraulics
sp.tr.m2<-as.matrix(sp.tr.m[,TRAITS.hydr])
fd<-dbFD(x=sp.tr.m2, a = PTAG.sp.m) ## note the dataset to include
a<-as.data.frame(fd$FDis)
colnames(a)<-"hydraulics.FDis"
FD.data<-cbind(FD.data, a)
################# calculate the CWM and FD value for the PC1 and PC2 
sp.tr.m2<-as.matrix(sp.tr.m[,"PC1.38"])
fd<-dbFD(x=sp.tr.m2, a = PTAG.sp.m) ## note the dataset to include
ab<-cbind(fd$FDis,fd$CWM)
colnames(ab)<-c("PC1.38.FDis","PC1.38.CWM")
FD.data<-cbind(FD.data, ab)
################# calculate the FD value for the PC1 and PC2 
sp.tr.m2<-as.matrix(sp.tr.m[,"PC1.les"])
fd<-dbFD(x=sp.tr.m2, a = PTAG.sp.m) ## note the dataset to include
ab<-cbind(fd$FDis,fd$CWM)
colnames(ab)<-c("PC1.les.FDis","PC1.les.CWM")
FD.data<-cbind(FD.data, ab)
################# calculate the FD value for the PC1 and PC2 
sp.tr.m2<-as.matrix(sp.tr.m[,"PC1.struc"])
fd<-dbFD(x=sp.tr.m2, a = PTAG.sp.m) ## note the dataset to include
ab<-cbind(fd$FDis,fd$CWM)
colnames(ab)<-c("PC1.struc.FDis","PC1.struc.CWM")
FD.data<-cbind(FD.data, ab)
################# calculate the FD value for the PC1 and PC2 
sp.tr.m2<-as.matrix(sp.tr.m[,"PC1.chem"])
fd<-dbFD(x=sp.tr.m2, a = PTAG.sp.m) ## note the dataset to include
ab<-cbind(fd$FDis,fd$CWM)
colnames(ab)<-c("PC1.chem.FDis","PC1.chem.CWM")
FD.data<-cbind(FD.data, ab)
################# calculate the FD value for the PC1 and PC2 
sp.tr.m2<-as.matrix(sp.tr.m[,"PC1.conduc"])
fd<-dbFD(x=sp.tr.m2, a = PTAG.sp.m) ## note the dataset to include
ab<-cbind(fd$FDis,fd$CWM)
colnames(ab)<-c("PC1.conduc.FDis","PC1.conduc.CWM")
FD.data<-cbind(FD.data, ab)
################# calculate the FD value for the PC1 and PC2 
sp.tr.m2<-as.matrix(sp.tr.m[,"PC1.hydr"])
fd<-dbFD(x=sp.tr.m2, a = PTAG.sp.m) ## note the dataset to include
ab<-cbind(fd$FDis,fd$CWM)
colnames(ab)<-c("PC1.hydr.FDis","PC1.hydr.CWM")
FD.data<-cbind(FD.data, ab)
}
names(FD.data)

######### create a replication of the FD data across all years to merge with the VOLUME dataset
df.FD<-FD.data
FD.all<-rbind(cbind("YEAR"=2010,df.FD),cbind("YEAR"=2011,df.FD),cbind("YEAR"=2012,df.FD),
              cbind("YEAR"=2013,df.FD),cbind("YEAR"=2014,df.FD),cbind("YEAR"=2015,df.FD),
              cbind("YEAR"=2016,df.FD),cbind("YEAR"=2017,df.FD),cbind("YEAR"=2018,df.FD),
              cbind("YEAR"=2019,df.FD),cbind("YEAR"=2020,df.FD))
table(FD.all[,c("YEAR")])

### merge the FD data to the df.478 dataset
names(df.478)
str(df.478)
df.new<-merge(df.478, FD.all, by=c("YEAR","PTAG"))
df.new$AGE<-ifelse(df.new$SITE=="A", df.new$YEAR-2009,df.new$YEAR-2010)
# Only select the plots with age 1--10 --> remove 2010 for SITE B and 2020 for SITE A - these have no data anyway
df.478.10<-droplevels(subset(df.new, AGE%in%c(1:10)))
table(df.478.10[,c("YEAR","SITE")])

###### ---- get scaled data for all CWM and FD values ranging from 0-1
## function to get normalized values between 0 and 1
norm.scale<-function(df){
  return( (df-min(df))/(max(df)-min(df)) )
}

#### calculate the normalized values of these new FD trait values to keep comparison
df.478.norm<-df.478.10
names(df.478.norm)
for(c in 19:length(df.478.10)){
  df.478.norm[,c]<-norm.scale(df.478.10[,c])
}

### get list with the multivariate FD traits
names(df.478.norm)
FD.traits<-as.vector(names(df.478.norm[19:length(df.478.norm)]))
FD.traits

###### perform all relationships to create the meta-data
## total 478 plots
df.analyses<-df.478.norm
plyr::count(FD.traits%in%names(df.analyses))
length(unique(df.analyses$PTAG))

# save data 
meta.table<-as.data.frame(matrix(nrow=0,ncol=15))
colnames(meta.table)<-c("Estimate","Std. Error","t value","Pr(>|t|)","R2","adj.R2","factor",
                        "Nptag","RESPONSE","MODEL","CAT","trait","trait.nr","TRAIT","AGE")
t=1;i=1;a=1;r=1;m=1
# run analyses for Standing volume (VOL.a) and for volume increment (VOLgr)
for(r in 1:1){ ### note to check all 3 response variables 
  if(r==1){RESPONSE<-"VOL.a";aa<-10}
  for(m in 1:1){
    if(m==1){MODEL<-"All"; set<-df.analyses }
    # run through the 2 different community indices
    for(i in 1:3){
      if(i==1){Tr<-FD.traits[1:7];CAT="multi.FD";Tr.short<-c("All (38)","categorical (2)","LES (6)","leaf structure (6)","leaf chemicals (6)","conductance (9)","hydraulics (9)")}
      if(i==2){Tr<-FD.traits[c(8,10,12,14,16,18)];CAT="PC1.FD";Tr.short<-c("All (38)","LES (6)","leaf structure (6)","leaf chemicals (6)","conductance (9)","hydraulics (9)")}
      if(i==3){Tr<-FD.traits[c(9,11,13,15,17,19)];CAT="PC1.CWM";Tr.short<-c("All (38)","LES (6)","leaf structure (6)","leaf chemicals (6)","conductance (9)","hydraulics (9)")}
      # run through the 38 traits
      for(t in 1:length(Tr)){
        sub<-set[,c("YEAR","AGE","SITE","TREE_R","PTAG",RESPONSE,Tr[t])]
        colnames(sub)[6:7]<-c("response","trait")# rename the columns
        # run through the 9 or 10 years
        for(a in 1:aa){  ## NOTE -> aa is set 10 or 9 depending on the response variable (VOL.a = standing volume, VOLgr = volume increment)
          sub1<-subset(sub,AGE==a)
          sub2<-sub1[complete.cases(sub1[,c("response","trait")]),]## cpmplete cases removes missing values for the response variable, due to missing data
          a1<-lm(response~trait, data = sub2)
          ### save the results of linear mode in datafile
          { 
            at<-as.data.frame(summary(a1)$coefficients)
            at$R2<-summary(a1)$r.squared
            at$adj.R2<-summary(a1)$adj.r.squared
            at$factor<-rownames(at)
            at$Nptag<-length(sub2$PTAG)
            at$RESPONSE<-RESPONSE
            at$MODEL<-MODEL
            at$CAT<-CAT
            at$trait<-Tr[t]
            at$trait.nr<-t
            at$TRAIT<-Tr.short[t]
            at$AGE<-a
            meta.table<-rbind(meta.table,at) 
          }
}}}}}
table(meta.table[meta.table$MODEL=="All",c("TRAIT","CAT")])
# remove the intercept data. no need for the meta-analyses
meta.data<-droplevels(subset(meta.table, factor=="trait"))
#write.csv(meta.data,"meta data Multivariate_original.csv")

####### --- use meta data of multivariate traits -----#########
names(meta.data)
table(meta.data[,c("TRAIT","CAT")])

TRAIT.group<-c("All (38)","categorical (2)","LES (6)","leaf structure (6)","leaf chemicals (6)","conductance (9)","hydraulics (9)")
meta.data$TRAIT.group <- factor(meta.data$TRAIT, levels = TRAIT.group )
table(meta.data[,c("TRAIT","TRAIT.group")])

### First only look at the FD value from mulitvariate trait using 38
dat<-subset(meta.data, MODEL=="All"&RESPONSE=="VOL.a")
### center AGE
dat$AGE2<-dat$AGE-mean(dat$AGE)
dat$AGE2.sq<-dat$AGE2^2
### effect size
m0<-lm(Estimate~CAT*TRAIT.group*(AGE2+AGE2.sq),data=dat)
anova(m0)
dat$Est.predict<-predict(m0)
### Reliability - R2
m1<-lm(asin(sqrt(R2))~CAT*TRAIT.group*(AGE2+AGE2.sq),data=dat)
anova(m1)
dat$predict.R2<-sin(predict(m1))^2

#### different kind of figure
TRAITS.group.col<-c("forestgreen","darksalmon","gold3","yellowgreen","steelblue4","firebrick")

sub<-droplevels(subset(dat, TRAIT!="All (38)"))
pEff.gr<-ggplot(
)+geom_point(data=sub, aes(x=(AGE), y=Estimate,col=TRAIT.group,shape=TRAIT.group),size=2,position=position_dodge(width=0.5)
)+geom_line(data = sub, aes(x=AGE,y=Est.predict,col=TRAIT.group)
)+facet_grid(.~CAT
)+theme_bw(
)+scale_color_manual(values=TRAITS.group.col
#)+scale_y_continuous(limits = c(-11,80), breaks=c(0,20,40,60),labels=c(0,20,40,60)
)+scale_x_continuous(limits = c(0.5,10.5), breaks=c(1:10),labels=c(1:10)
)+theme(panel.grid = element_blank(),legend.position = c(0.12,0.72),legend.title = element_blank()
)+labs(x="Stand age (year)",y="Trait effect on stand volume (slope)",title="b"
)+theme (plot.title = element_text (hjust = -0.08, vjust=-5)
)
pEff.gr

pR2.gr<-ggplot(
)+geom_point(data=sub, aes(x=(AGE), y=R2,col=TRAIT.group,shape=TRAIT.group),size=2,position=position_dodge(width=0.5)
)+geom_line(data = sub, aes(x=AGE,y=predict.R2,col=TRAIT.group)
)+facet_grid(.~CAT
)+theme_bw(
)+scale_color_manual(values=TRAITS.group.col
)+scale_y_continuous(limits = c(-0.005,0.15)#, breaks=c(0,0.05,0.1,0.15),labels=c(0,0.05,0.1,0.15)
)+scale_x_continuous(limits = c(0.5,10.5), breaks=c(1:10),labels=c(1:10)
)+theme(panel.grid = element_blank(),legend.position = "none",legend.title = element_blank()
)+labs(y=expression("Reliability to predict stand volume ("~R^2~")"),x="Stand age (year)",title="d"
)+theme (plot.title = element_text (hjust = -0.08, vjust=-5)
)
pR2.gr

######################################################################################
##### --------- create meta-data with TREE_R as factor -------########################
df.analyses<-df.478
df.analyses$TREE_log2<-log2(df.analyses$TREE_R)
### normalize the richness effect because to be able to compare
norm.scale<-function(df){
  return( (df-min(df))/(max(df)-min(df)) )
}
df.analyses$TREE_norm<-norm.scale(df.analyses$TREE_log2)

# save data in this plot
meta.table<-as.data.frame(matrix(nrow=0,ncol=15))
colnames(meta.table)<-c("Estimate","Std. Error","t value","Pr(>|t|)","R2","adj.R2","factor",
                        "Nptag","RESPONSE","MODEL","CAT","trait","trait.nr","TRAIT","AGE")
RESPONSE<-"VOL.a";aa<-10
MODEL<-"All"; set<-df.analyses
for(a in 1:aa){  ## NOTE -> aa is set 10 or 9 depending on the response variable (VOL.a = standing volume, VOLgr = volume increment)
      sub<-set[set$AGE==a,c("YEAR","AGE","SITE","TREE_norm","PTAG",RESPONSE)]
      colnames(sub)[6]<-c("response")# rename the columns
      sub2<-sub[complete.cases(sub[,c("response")]),]## complete cases removes missing values for the response variable, due to missing data
      a1<-lm(response~TREE_norm, data = sub2)
      ### save the model in datafile
      { 
        at<-as.data.frame(summary(a1)$coefficients)
        at$R2<-summary(a1)$r.squared
        at$adj.R2<-summary(a1)$adj.r.squared
        at$factor<-rownames(at)
        at$Nptag<-length(sub2$PTAG)
        at$RESPONSE<-RESPONSE
        at$MODEL<-MODEL
        at$CAT<-"Richness"
        at$trait<-"log2(Richness)"
        at$trait.nr<-t
        at$TRAIT<-"Richness"
        at$AGE<-a
        meta.table<-rbind(meta.table,at) 
      }
}
table(meta.table[meta.table$MODEL=="All",c("RESPONSE","CAT")])
# remove the intercept data. no need for the meta-analyses
meta.data_richness<-droplevels(subset(meta.table, factor!="(Intercept)"))
##############################################################################################

#############################################################################################
#### ------ make a figure which combines the ALL 38 trait data and the richness ------######
### all 38 traits data
head(meta.data)
dat.all38<-droplevels(subset(meta.data,RESPONSE=="VOL.a"&MODEL=="All"&TRAIT=="All (38)"))
### center AGE
dat.all38$AGE2<-dat.all38$AGE-mean(dat.all38$AGE)
dat.all38$AGE2.sq<-dat.all38$AGE2^2
### effect size
m0<-lm(Estimate~CAT*(AGE2+AGE2.sq),data=dat.all38)
anova(m0)
dat.all38$Est.predict<-predict(m0)
### Reliability - R2
m1<-lm(asin(sqrt(R2))~CAT*(AGE2+AGE2.sq),data=dat.all38)
anova(m1)
dat.all38$predict.R2<-sin(predict(m1))^2
head(dat.all38)

### RIchness effect data
head(meta.data_richness)
### center AGE
meta.data_richness$AGE2<-meta.data_richness$AGE-mean(meta.data_richness$AGE)
meta.data_richness$AGE2.sq<-meta.data_richness$AGE2^2
### effect size
m0<-lm(Estimate~AGE2+AGE2.sq,data=meta.data_richness)
anova(m0)
meta.data_richness$Est.predict<-predict(m0)
### Reliability - R2
m1<-lm(asin(sqrt(R2))~AGE2+AGE2.sq,data=meta.data_richness)
anova(m1)
meta.data_richness$predict.R2<-sin(predict(m1))^2

####Merge the data of richness with the All 38 data
dat.merge<-rbind(dat.all38[,c("Estimate","R2","Est.predict","predict.R2","AGE","AGE2","AGE2.sq","CAT","TRAIT")],
                 meta.data_richness[,c("Estimate","R2","Est.predict","predict.R2","AGE","AGE2","AGE2.sq","CAT","TRAIT")])
table(dat.merge[,c("TRAIT","CAT")])

pEff<-ggplot(data=dat.merge
)+geom_point(aes(x=AGE, y=Estimate,col=CAT,shape=CAT),size=2,position=position_dodge(width=0.5)
)+geom_line(aes(x=AGE,y=Est.predict,col=CAT,lty=CAT)
)+facet_grid(.~TRAIT
)+theme_bw(
)+scale_color_manual(values=c("brown",diverge_hcl(2,h=c(250,40),c=96),"black")
)+scale_shape_manual(values=c(2,16,17,15)
)+scale_linetype_manual(values=c(2,1,1,1)
)+scale_y_continuous(limits = c(-11,60), breaks=c(0,20,40,60),labels=c(0,20,40,60)
)+scale_x_continuous(limits = c(0.5,10.5), breaks=c(1:10),labels=c(1:10)
)+theme(panel.grid = element_blank(),legend.position = c(0.65,0.8),legend.title = element_blank()
)+labs(x="Stand age (year)",y="Trait effect on stand volume (slope)",title="a"
)+theme (plot.title = element_text (hjust = -0.1, vjust=-5)
)
pEff

pR2<-ggplot(data=dat.merge
)+geom_point(aes(x=AGE, y=R2,col=CAT,shape=CAT),size=2,position=position_dodge(width=0.5)
)+geom_line(aes(x=AGE,y=predict.R2,col=CAT,lty=CAT)
)+facet_grid(.~TRAIT
)+theme_bw(
)+scale_color_manual(values=c("brown",diverge_hcl(2,h=c(250,40),c=96),"black")
)+scale_shape_manual(values=c(2,16,17,15)
)+scale_linetype_manual(values=c(2,1,1,1)
)+scale_y_continuous(limits = c(-0.005,0.15)#, breaks=c(0,0.05,0.1,0.15),labels=c(0,0.05,0.1,0.15)
)+scale_x_continuous(limits = c(0.5,10.5), breaks=c(1:10),labels=c(1:10)
)+theme(panel.grid = element_blank(),legend.position = "none",legend.title = element_blank()
)+labs(y=expression("Reliability to predict stand volume ("~R^2~")"),x="Stand age (year)",title="c"
)+theme (plot.title = element_text (hjust = -0.15, vjust=-5)
)
pR2

tiff("SUPP Effect and R2 for multi FD PCA Richness.tiff",width=3500,height=2500,unit="px",res=300)
grid.arrange(pEff,pEff.gr,pR2,pR2.gr,layout_matrix = cbind(c(1,3), c(1,3),c(2,4),c(2,4),c(2,4)) )
dev.off()
