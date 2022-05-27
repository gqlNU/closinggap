# Association of environmental and socioeconomic indicators with serious mental illness diagnoses identified from general practitioner practice data in England: A spatial Bayesian modelling study


## Data
- The file `smi_closingthegap.RData` contains the SMI prevalence and the associated covariates at the LSOA level;
- The folder `England` contains the graph files for MSOAs, Districts and CCGs. Each graph file defines the neighbourhood structure at that geographical level. These graph files are INLA inputs.
 
##  R code for model fitting

The code below assumes all the above files are stored under the following folder structure:
- `C:/closinggap/smi_closingthegap.RData`
- `C:/closinggap/England/CCGspadj.graph`
- `C:/closinggap/England/LADspadj.graph`
- `C:/closinggap/England/MSOAspadj.graph`


```R
library(foreign)
library(INLA)
library(sf)
library(ggplot2)
library(viridis)
library(gridExtra)
library(maptools)
####################
#  user input
####################
exclude.outliers <- TRUE

###  end user input

######################################################################
###   set working directory
###   This is the full path to the folder where the data file
###   and the graph files are stored.
######################################################################
workdir <- 'C:/closinggap/'
setwd(workdir)

######################################################################
###   load the dataset
######################################################################
load('smi_closingthegap.RData')
ccg.iid <- indices$ccg.iid
lad.iid <- indices$lad.iid
msoa.iid <- indices$msoa.iid
ccg.icar <- indices$ccg.icar
lad.icar <- indices$lad.icar
msoa.icar <- indices$msoa.icar

######################################################################
###   LSOAs with outlying values (see thge main paper for detail)
######################################################################
outliers <- c('Kensington and Chelsea 020A','Tower Hamlets 015A','Westminster 019B','Westminster 020C'
             ,'Sefton 021D','Sheffield 036A','North Devon 002D','East Staffordshire 006B','Sheffield 073E','Nottingham 026G')
if (exclude.outliers) {
	#  replace the log SMI values of the above LSOAs by NA
	ids <- sapply(outliers,function(x){which(fitdata$lsoa11nm==x)})
	fitdata$log.SMI_mean[ids] <- NA
}

######################################################################
###   define the full model
######################################################################
formula <- log.SMI_mean ~ 1 + 
PGS_Area + PGS_dist + wood_ha +
    nonWhite_pc +
    PPM25 +
    A18_24pc +
    A25_44pc +
    A45_64pc +
    A65overpc +
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
    f(ccg.iid,model='iid',graph='England/CCGspadj.graph') +
    f(ccg.icar,model='besag',graph='England/CCGspadj.graph') +
    f(lad.iid,model='iid',graph='England/LADspadj.graph') +
    f(lad.icar,model='besag',graph='England/LADspadj.graph') +
    f(msoa.iid,model='iid',graph='England/MSOAspadj.graph') +
    f(msoa.icar,model='besag',graph='England/MSOAspadj.graph')

######################################################################
###   model fitting
######################################################################
start.time <- Sys.time()
print(start.time)
m <- inla(formula,family='gaussian',data=fitdata,control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE,config=FALSE))
end.time <- Sys.time()
print(paste0('model fitting took ',end.time-start.time,' minutes'))
summary(m,digits=6)

```
