##### Project number 2
##### Community trait effects on productivity across time

######
rm(list=ls())
library(ggplot2)
library(plyr)#count()
library(dplyr)
library(reshape2)#dcast() / melt() data

#####################################################################################
## --------------------------------- LOAD DATA --------------------------------######
#####################################################################################
#### corrected individual data for site A and B
setwd("c:\\Users/Thinkpad/Documents/IBCAS/P2 Community and traits/")
SP2<-read.csv("Individual data 2010-2020_fjb.csv")

######################################################################################
### ----------  calculate the SUM of crown area  --------------------------------#####
###------------> measure for canopy closure -------------------------------------#####
#######################################################################################
names(SP2)
plyr::count(SP2$CELL_X)
#### use only the middel 36 plants for plot-level calculations
## middel 36 -> x and y coord from 7-12
SP2.sub<-droplevels(subset(SP2,middle36=="middle"))
plyr::count(SP2.sub$CELL_Y)

### estimate the canopy area per individual tree
{SP2.sub$CA10<-pi*(SP2.sub$CR10^2)
SP2.sub$CA11<-pi*(SP2.sub$CR11^2)
SP2.sub$CA12<-pi*(SP2.sub$CR12^2)
SP2.sub$CA13<-pi*(SP2.sub$CR13^2)
SP2.sub$CA14<-pi*(SP2.sub$CR14^2)
SP2.sub$CA15<-pi*(SP2.sub$CR15^2)
SP2.sub$CA16<-pi*(SP2.sub$CR16^2)
SP2.sub$CA17<-pi*(SP2.sub$CR17^2)
SP2.sub$CA18<-pi*(SP2.sub$CR18^2)
SP2.sub$CA19<-pi*(SP2.sub$CR19^2)
SP2.sub$CA20<-pi*(SP2.sub$CR20^2)
}

SP2.CA<-aggregate(SP2.sub[,c("CA10","CA11","CA12","CA13","CA14","CA15","CA16","CA17","CA18","CA19","CA20")],
                  by=list(SITE=SP2.sub$SITE,PTAG=SP2.sub$PTAG,TREE_R=SP2.sub$TREE_R),FUN=sum, na.rm=T)
head(SP2.CA)

#### put the different year columns in one single column
c1<-c(paste("CA",c(10:20),sep=""))
df.lai<-melt(SP2.CA[,c("SITE","PTAG","TREE_R",c1)],id=c("SITE","PTAG","TREE_R"),
          variable.name = "YEAR.f",value.name = "CA",all=T)
df.lai$LAI<-(df.lai$CA/(666.7/400*36))

### create a YEAR column
year<-data.frame("YEAR"=c(2010:2020),"c1"=c1)
df.lai$YEAR<-year$YEAR[match(df.lai$YEAR.f,year$c1)]
df.lai$AGE<-ifelse(df.lai$SITE=="A",df.lai$YEAR-2009,df.lai$YEAR-2010)

df.sub<-droplevels(subset(df.lai, AGE%in%c(1:10)))
df.sub$LAI.f<-ifelse(df.sub$LAI<1,0,1)
### create NA values for plots with no data -> use CA==NA
df.sub$LAI<-ifelse(is.na(df.sub$CA),NA,df.sub$LAI)
df.sub$LAI.f<-ifelse(df.sub$LAI<1,0,1)

nm<-glm(LAI.f~AGE,data = df.sub,family=binomial(link=logit))
summary(nm)
coef(nm)
AIC(nm)
df.sub$fitted.lai<-fitted(nm)

p<-ggplot(data = df.sub)+geom_line(aes(x=AGE,y=fitted.lai)
)+geom_hline(yintercept=0.5,lty=2)+geom_hline(yintercept=0.75,lty=2
)+theme_bw()+labs(x="Stand age (year)",y="Fitted probability of plots having a closed canopy"
)+scale_x_continuous(breaks=c(1:10))+ylim(0,1)

#tiff("FIG S6 Estimation of fraction of plots with closed canopy.tiff",res=600,unit="px",width=2500,height=2500)
print(p)
dev.off()
