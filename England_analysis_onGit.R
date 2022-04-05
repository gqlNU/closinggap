#############################################################
####   this script fits the full model to the England data
#############################################################
rm(list=ls());gc()
library(foreign)
library(INLA)
library(sf)
library(ggplot2)
library(viridis)
library(gridExtra)
workdir <- '/Volumes/WorkSpace/ClosingGap/onGit/geo_data'
setwd(workdir)

######################################################################
###   load the shapefile into R
######################################################################
load("formatted_closing_gap_england_09Dec2019.RData")
shp <- england
shp.sf <- st_as_sf(shp)
shp.sf$log.SMI_mean <- log(shp.sf$SMI_mean)

data <- slot(shp,'data')#read.dbf('dataset_20191018.dbf')
data$lsoa11cd <- as.character(data$lsoa11cd)
data$lsoa11nm <- as.character(data$lsoa11nm)
data$log.SMI_mean <- shp.sf$log.SMI_mean
data$Nr_GPPHML1 <- as.numeric(as.character(data$Nr_GPPHML1))

######################################################################
###   load the updated data with FRZ3_dist (flood risk 3*)
######################################################################
updated <- read.csv('10Jun2020/dataset_20200610.csv',header=TRUE)
updated$LSOA11NM <- as.character(updated$LSOA11NM)
ids <- sapply(data$lsoa11nm,function(x){which(updated$LSOA11NM==x)})
data$FRZ3_dist <- updated$FRZ3_dist[ids]

######################################################################
###   load the updated data with FRZ3_dist (flood risk 3*) no defence
######################################################################
updated <- read.csv('25Jun2020/dataset_20200624.csv',header=TRUE)
updated$LSOA11NM <- as.character(updated$LSOA11NM)
ids <- sapply(data$lsoa11nm,function(x){which(updated$LSOA11NM==x)})
data$FRZ3_dist_no_defence <- updated$FR3_no.defences_dist[ids]

######################################################################
###   get %non-white
######################################################################
data$nonWhite_pc <- 100-data$W_pc

######################################################################
###  constructing random effect IDs
######################################################################
###  LSOA
lsoa.icar <- 1:nrow(data)

###  LAD
data$LAD17NM <- as.character(data$LAD17NM)
uniq <- unique(data$LAD17NM)
nlads <- length(uniq)
id.lad <- nrow(data)
for (i in 1:nlads) {
  id.lad[which(data$LAD17NM==uniq[i])] <- i
}
lad.iid <- lad.icar <- id.lad

###  CCG
lookup.file <- 'lookup/LSOA_2011_to_Clinical_Commissioning_Groups_to_Sustainability_and_Transformation_Partnerships_April_2017_Lookup_in_England.csv'
lookup <- read.csv(lookup.file,header=TRUE)
lookup$LSOA11CD <- as.character(lookup$LSOA11CD)
lookup$CCG17CD <- as.character(lookup$CCG17CD)
lsoa11cd <- data$lsoa11cd
ids <- match(lsoa11cd,lookup$LSOA11CD)
data$CCG17CD <- lookup$CCG17CD[ids]
uniq <- unique(data$CCG17CD)
load('graphs/CCGspadj.RData')
summary(sapply(names(ccg.spadj$num),function(x){which(uniq==x)})-1:length(uniq))
nccgs <- length(uniq)
id.ccg <- nrow(data)
for (i in 1:nccgs) {
	id.ccg[which(data$CCG17CD==uniq[i])] <- i
}
ccg.iid <- ccg.icar <- id.ccg

###  MSOA
data$MSOA11NM <- as.character(data$MSOA11NM)
uniq <- unique(data$MSOA11NM)
nr <- length(uniq)
id.msoa.iid <- nrow(data)
for (i in 1:nr) {
  id.msoa.iid[which(data$MSOA11NM ==uniq[i])] <- i
}
msoa.icar <- msoa.iid <- id.msoa.iid

###  region
data$RGN11NM <- as.character(data$RGN11NM)
uniq <- unique(data$RGN11NM)
nr <- length(uniq)
id.region <- nrow(data)
for (i in 1:nr) {
  id.region[which(data$RGN11NM==uniq[i])] <- i
}

#######################################################################################
#########    this changes the original by
#########    1. replace flood 2 by flood 3
#########    2. both W and BL removed but %non-white added
#########    3. Nr_GPPHML1 removed
#########    4. remove the 10 outliers
########################################################################################
###   randome effects (MSOA-BYM/CCG-ICAR/LAD-BYM)     NO LSOA
###     + IMD + AGE + ethnicity
###     + GP + urbanicity + flood + noise + greenspace
###     + PM2.5
###     + bluespace
###     + CCG
###     + region becomes categorical (only 10 regions)
if (TRUE) {
outliers <- c('Kensington and Chelsea 020A','Tower Hamlets 015A','Westminster 019B','Westminster 020C'
             ,'Sefton 021D','Sheffield 036A','North Devon 002D','East Staffordshire 006B','Sheffield 073E','Nottingham 026G')
ids <- sapply(outliers,function(x){which(data$lsoa11nm==x)})
data$log.SMI_mean[ids] <- NA
formula <- log.SMI_mean ~ 1 + 
	PGS_pc + PGS_dist + wood_ha +
    nonWhite_pc +
    PPM25 +
    A18_24 +
    A25_44 +
    A45_64 +
    A65over +
    Noise_5 +
    FRZ3_dist +
    CRIquintile +
    INCquintile +
    BHSquintile +
    ASquintile +
    EMPquintile +
    INDquintile +  
    Type_urban +   
    LakePGS_di + RiverPGS_d +
    RGN11NM +
    f(ccg.iid,model='iid') +
    f(ccg.icar,model='besag',graph='graphs/CCGspadj.graph') +
    f(lad.iid,model='iid') +
    f(lad.icar,model='besag',graph='graphs/LADspadj.graph') +
    f(msoa.iid,model='iid') +
    f(msoa.icar,model='besag',graph='graphs/MSOAspadj.graph')
start.time <- Sys.time()
m <- inla(formula,family='gaussian',data=data,control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE,config=FALSE))
end.time <- Sys.time()
print(paste0('model fitting took ',end.time-start.time,' minutes'))
summary(m,digits=6)
save(file='results/England/LAD-BYM_MSOA-BYM_CCG-BYM_IMD_age_nonWhitePC_urban_noise_flood3_PGS_PM25_bluespace_RegionCAT_outliersRemoved.RData',m,data)
}

if (FALSE) {
###     random effect only model	
outliers <- c('Kensington and Chelsea 020A','Tower Hamlets 015A','Westminster 019B','Westminster 020C'
             ,'Sefton 021D','Sheffield 036A','North Devon 002D','East Staffordshire 006B','Sheffield 073E','Nottingham 026G')
ids <- sapply(outliers,function(x){which(data$lsoa11nm==x)})
data$log.SMI_mean[ids] <- NA
formula <- log.SMI_mean ~ 1 + 
    RGN11NM +
    f(ccg.iid,model='iid') +
    f(ccg.icar,model='besag',graph='England/CCGspadj.graph') +
    f(lad.iid,model='iid') +
    #f(lad.icar,model='besag',graph='England/LADspadj.graph') +
    f(msoa.iid,model='iid') +
    f(msoa.icar,model='besag',graph='England/MSOAspadj.graph')
start.time <- Sys.time()
m <- inla(formula,family='gaussian',data=data,control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE,config=FALSE))
end.time <- Sys.time()
print(paste0('model fitting took ',end.time-start.time,' minutes'))
summary(m,digits=6)
save(file='results/England/LAD-IID_MSOA-BYM_CCG-BYM_only_outliersRemoved.RData',m,data)
}
