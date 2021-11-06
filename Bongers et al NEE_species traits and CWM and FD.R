##################################################################################################################
##### R code belonging to Bongers et al 2021 Nature Ecology and Evolotion
##### Functional diversity effects on productivity increase with age in a forest biodiveristy experiment
#################################################################################################################

rm(list=ls())
#required libraries
{
library(ggplot2)
library(plyr)#count # load this package first to keep functions well
library(dplyr)
library(reshape2)#dcast() / melt() data
library(psych)#principal() = PCA / cor.plot()
library(gridExtra)
library(readxl)## import sheets of excel files
library(FD)# functional diversity function
}
####################################################################################
######## ------------------LOAD THE DATA ---------------------------------##########
setwd()
df.design<-read.csv("Bongers et al NEE_species proportions.csv")
str(df.design)

df.traits<-as.data.frame(read_excel("Bongers et al NEE_Species trait values.xlsx",sheet="Species mean trait values",na=""))
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
TRAITS%in%names(df.traits)[3:40]

####################################################################################################
######## TRAIT CORRELATIONS AMONG SPECIES
#### panel.cor and hist functions
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r1=cor(x,y,use="pairwise.complete.obs")
  r =0.8 ## all values should have the same size
  txt <- format(c(r1, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.9/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * r)
}
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="white", ...)
}
tiff("SPECIES trait correlations matrix.tiff",unit="px",res=300,width=12000,height=12000)
print(pairs(df.traits[,TRAITS], gap=2/10, upper.panel=panel.cor, lower.panel=points,diag.panel = panel.hist))
dev.off()

###################################################################################################
####### ---------------- CALCULATE CWM AND FD for all traits ------------------------------########
# use FD function from FD package "laliberte and Legendre"
# FD funtion can be used by combining dataset of all species(row) and all traits(columns) 
# and a dataset with species abundance per different plots -> plots in rows, species in columns
# NOTE the species have to be identical 
###################################################################################################
names(df.traits)
names(df.design)

######### FIRST
# create the species-trait matrix --> NOTE the trait data used
names(df.traits)
df.traits.sel<-df.traits[,c("species_name",TRAITS)]
SPECIES.traits<-as.vector(df.traits.sel$species_name)
SPECIES.traits
rownames(df.traits.sel)<-df.traits.sel$species_name
head(df.traits.sel)
sp.tr.m<-as.matrix(df.traits.sel[,TRAITS])

######## SECOND
# create the PTAG - species abundance matrix 
abun.sp<-dcast(df.design[,c("PTAG","species_name","rel.abun.des")],PTAG~species_name)
length(abun.sp)
length(abun.sp$PTAG)
## save PTAg in rownames
rownames(abun.sp)<-abun.sp$PTAG
PTAG.list<-as.vector(abun.sp$PTAG)
PTAG.sp.m<-as.matrix(abun.sp[,SPECIES.traits])

##################################################################################
###### CALCULATE THE FD and CWM VALUES
##      require(FD)
##################################################################################
### create dataset to save this data
FD.data<-as.data.frame(cbind(PTAG=PTAG.list))
head(FD.data)
##### apply function of FD
t=38;TRAITS[t]
for(t in 1:length(TRAITS)){
  sp.tr.m2<-as.matrix(sp.tr.m[,TRAITS[t]])
  head(sp.tr.m2)
  colnames(sp.tr.m2)<-TRAITS[t]
  fd<-dbFD(x=sp.tr.m2, a = PTAG.sp.m, CWM.type="all") ## note the dataset to include
  names(fd)
  a<-as.data.frame(fd$FRic)
  b<-as.data.frame(fd$FDis)
  c<-as.data.frame(fd$FEve)
  d<-as.data.frame(fd$RaoQ)
  e<-as.data.frame(fd$CWM)
  set<-cbind(a,b,c,d,e)
  head(set)
  ## create a data file with CWM/FD per trait 
  if(t<=2)(colnames(set)<-gsub("\\s+","",paste(TRAITS[t],c(".FRic",".FDis",".FEve",".RaoQ","0.CWM",".CWM")))
           )else(colnames(set)<-gsub("\\s+","",paste(TRAITS[t],c(".FRic",".FDis",".FEve",".RaoQ",".CWM"))))
  FD.data<-cbind(FD.data, set)
}
head(FD.data)
write.csv(FD.data,"CWM FD data 38 single traits.csv")

