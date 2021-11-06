##################################################################################################################
##### R code belonging to Bongers et al 2021 Nature Ecology and Evolotion
##### Functional diversity effects on productivity increase with age in a forest biodiveristy experiment
####  Franca J Bongers et al
#################################################################################################################

rm(list=ls())
{
library(ggplot2)
library(plyr)#count()
library(dplyr)
library(reshape2)#dcast() / melt() data
library(gridExtra)
library(lme4)
library(lmerTest)
library(car)
library(readxl)## import sheets of excel files
}

### trait labels used in the dataset
TRAITS<-c("D","AM", ## catergorical traits
          "SLA","LOG10LA","LDMC","LNC","LCC","CN", # LES
          "LEAFTHICK","LEAFTOUGH","UPPEREPI", "PALIS","SPONGY", "LOG10RATIO", #structure
          "P","Phen","Tan", "Ca",  "K", "Mg",#chemicals
          "STOMDENS", "STOMSIZE", "STOMIND", "CONMEAN","CONMAX",
          "CONMAXFIT","VPDMAX",  "VPDMAXFIT",  "VPDPOI", # conductance traits
          "WD","HYDCOND","PSI50","WPOT","VEINDENS",
          "DIAMVEIN","MEANAREA","MEANROUND",  "DHYD") # hydraulic traits

TRAITS.list<-c(paste(TRAITS,"CWM",sep="."),paste(TRAITS,"FDis",sep="."))

########################################################################################
######## ------------------LOAD THE DATA ---------------------------------##############
setwd()
### VOLUME and INCREMENT data
df.478<-as.data.frame(read_excel("Data_Bongers et al NEE_stand volume and increment.xlsx",sheet="Raw data", na="NA"))
names(df.478)
str(df.478)

#### CWM and FD data
df.FD<-read.csv("CWM FD data 38 single traits_original.csv") #import the CWM and FD data you got from the FD package using the 38 traits
head(df.FD)[1:10]

#### merge the data of FD.data with the df.478 -> NOTE that FD.data needs to be replicated for 10 years
FD.all<-rbind(cbind("YEAR"=2010,df.FD),cbind("YEAR"=2011,df.FD),cbind("YEAR"=2012,df.FD),
              cbind("YEAR"=2013,df.FD),cbind("YEAR"=2014,df.FD),cbind("YEAR"=2015,df.FD),
              cbind("YEAR"=2016,df.FD),cbind("YEAR"=2017,df.FD),cbind("YEAR"=2018,df.FD),
              cbind("YEAR"=2019,df.FD),cbind("YEAR"=2020,df.FD))
table(FD.all[,c("YEAR")])

### merge the FD data to the df.478 dataset
df.new<-merge(df.478, FD.all, by=c("YEAR","PTAG"))
df.new$AGE<-ifelse(df.new$SITE=="A", df.new$YEAR-2009,df.new$YEAR-2010)
# Only select the plots with age 1--10 --> remove 2010 for SITE B and 2020 for SITE A
df.new.10<-droplevels(subset(df.new, AGE%in%c(1:10)))
table(df.new.10[,c("AGE","SITE")])

############# --------- create figure with the relationships for age 5 and 10 ------###########
tiff("SUPP VOLUME FD relationships for AGE 10.tiff",width=6000,height=9000,unit="px",res=300)
set<-droplevels(subset(df.new.10,AGE==10))
par(mfrow=c(8,5))
for(t in 1:length(TRAITS)){
  sub<-set[,c("PTAG","AGE","VOL.a",paste(TRAITS[t],".FDis",sep=""))]
  colnames(sub)<-c("PTAG","AGE","resp","trait")
  lm1<-lm(resp~trait,data=sub)
  plot(resp~trait,data=sub,ylab="Stand volume",xlab=paste("FD -",TRAITS[t]),main="Stand age 10 years")
  curve(expr=summary(lm1)$coeff[1,1]+summary(lm1)$coeff[2,1]*x,add=T,lty=ifelse(summary(lm1)$coeff[2,4]<0.05,1,2))
}
dev.off()
#######################################################################################

############### ---- get scaled data for all CWM and FD values ranging from 0-1-----------#############
## function to get normalized values between 0 and 1
norm.scale<-function(df){
  return( (df-min(df))/(max(df)-min(df)) )
}

## make a new dataframe in which all values of CWM, FDis, FRic and RaoQ are replaced with the normalized values
df.new.norm<-df.new.10
names(df.new.norm)
for(c in 20:length(df.new.10)){
  df.new.norm[,c]<-norm.scale(df.new.10[,c])
}
names(df.new.norm)

##########################################################################################################################
#### CREATING META-DATA -> Productivity - Trait relationships per year
########################################################################################################################### 
# select the dataset for analyses
## total 478 plots
df.analyses<-df.new.norm
plyr::count(TRAITS.list%in%names(df.analyses))
length(unique(df.analyses$PTAG))

# save data in this plot
meta.table<-as.data.frame(matrix(nrow=0,ncol=15))
colnames(meta.table)<-c("Estimate","Std. Error","t value","Pr(>|t|)","R2","adj.R2","factor",
                      "Nptag","RESPONSE","MODEL","CAT","trait","trait.nr","TRAIT","AGE")
t=1;i=1;a=1;r=1
# run analyses for Standing volume (VOL.a) and for volume increment (VOLgr)
for(r in 1:2){ 
  if(r==1){RESPONSE<-"VOL.a";aa<-10}
  if(r==2){RESPONSE<-"VOLgr";aa<-9}
for(m in 1:11){
  if(m==1){MODEL<-"All"; set<-df.analyses }
  if(m==2){MODEL<-"A.A";set<-subset(df.analyses,A.A==1)}
  if(m==3){MODEL<-"A.B";set<-subset(df.analyses,A.B==1)}
  if(m==4){MODEL<-"A.C";set<-subset(df.analyses,A.C==1)}
  if(m==5){MODEL<-"B.A";set<-subset(df.analyses,B.A==1)}
  if(m==6){MODEL<-"B.B";set<-subset(df.analyses,B.B==1)}
  if(m==7){MODEL<-"B.C";set<-subset(df.analyses,B.C==1)}
  if(m==8){MODEL<-"A.SLA";set<-subset(df.analyses,A.SLA==1)}
  if(m==9){MODEL<-"A.Rar";set<-subset(df.analyses,A.Rar==1)}
  if(m==10){MODEL<-"B.SLA";set<-subset(df.analyses,B.SLA==1)}
  if(m==11){MODEL<-"B.Rar";set<-subset(df.analyses,B.Rar==1)}  
  
# run through the 2 different community indices
for(i in 1:2){
if(i==1){Tr<-c(paste(TRAITS,"CWM",sep="."));CAT="CWM"} 
if(i==2){Tr<-c(paste(TRAITS,"FDis",sep="."));CAT="FD"}
# run through the 38 traits
for(t in 1:length(Tr)){
  sub<-set[,c("YEAR","AGE","SITE","TREE_R","PTAG",RESPONSE,Tr[t])]
  colnames(sub)[6:7]<-c("response","trait")# rename the columns
# run through the 9 or 10 years
for(a in 1:aa){  ## NOTE -> aa is set 10 or 9 depending on the response variable (VOL.a = standing volume, VOLgr = volume increment)
  sub1<-subset(sub,AGE==a)
  sub2<-sub1[complete.cases(sub1[,c("response","trait")]),]## cpmplete cases removes missing values for the response variable, due to missing data
  a1<-lm(response~trait, data = sub2)
  ### save the model in datafile
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
      at$TRAIT<-TRAITS[t]
      at$AGE<-a
      meta.table<-rbind(meta.table,at) 
  }
}}}}}
head(meta.table)
# remove the intercept data. no need for the meta-analyses
meta.data<-droplevels(subset(meta.table, factor=="trait"))
table(meta.data[meta.data$MODEL=="All",c("RESPONSE","CAT")])
write.csv(meta.data,"Meta data_original.csv")


