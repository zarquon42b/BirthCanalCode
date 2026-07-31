library(Morpho)
library(geomorph)
library(shapes)
library(sf)
library(lattice)
library(Rvcg)
library(DescTools)
library(devtools)
library(matchingR)

#setwd("~/Desktop/Primate pelvis project/Birth canal shape and neonatal size/OUTLET/R_scripts & data/New analyses July 2026")

#Load data
Coords <- read.morphologika("1_130_Females.txt")
Mean.coord<-read.morphologika("Pelvis_mean_coords.txt")
List <- read.csv("1_IndividList.csv")
List$Species<-as.factor(List$Species)
MeanPelvesID<-read.csv("1_MeanPelvesList.csv")
MeanPelvesID$Species<-as.factor(MeanPelvesID$Species)
Neonate.data<-read.csv("1_Species_neonatal_data.csv")

# Specify paired (L-R) landmarks
left <- c(1,3,5,7,9,11,13,14,15,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243)
right <- c(2,4,6,8,10,12,13,14,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268)
pairedLMs <- cbind(left,right)
data_sym <- symmetrize(Coords,pairedLM=pairedLMs)

#removing redundant semilandmarks that overlap with landmarks

Coords_filtered <- data_sym[-c(54,83,84,113,114,128,129,
                               143,144,158,159,173,174,188,
                               189,203,204,218,219,243,244,
                               268), , ]

#r2morphologika(Coords_filtered, file="130_246lm_females.txt")

Coords_filtered <- procSym(Coords_filtered,CSinit=F,scale=F,pcAlign=T)$rotated

# Traditional inlets

Species.list<-levels(List$Species) #list of species

# Resample sacral landmarks

sacral.lm<-c(13,110:122,15)
resampled.sacral.lm<-list(NULL)

for (i in 1:nlevels(List$Species)){ # loop to resample sacral landmarks so that we have 59
    resampled.sacral.lm[[i]]<-resampleCurve(Mean.coord[sacral.lm,,MeanPelvesID$Species==levels(List$Species)[i]],59) # resample sacral landmarks
  } 

names(resampled.sacral.lm) <- levels(List$Species)

## puboischial landmarks on each side

puboischial.lmL<-list(NULL) # Left side
puboischial.lmR<-list(NULL) # Right side
for (i in 1:nlevels(List$Species)){
    j <- grep(levels(List$Species)[[i]],dimnames(Mean.coord)[[3]])
  puboischial.lmL[[i]] <- rbind(resampleCurve(Mean.coord[c(28,123:135,30),,j],30),resampleCurve(Mean.coord[c(30,175:187,42),,j],30)[-1,]) #resampling pubic symphysis landmarks to 30 and adding resampled ischio-pubic landmarks (excluding the shared fixed landmark at the bottom of the pubis)
  puboischial.lmR[[i]] <- rbind(resampleCurve(Mean.coord[c(29,136:148,31),,j],30),resampleCurve(Mean.coord[c(31,188:200,43),,j],30)[-1,])
    Mean.coord[c(29,136:148,31,188:200,43),,j]
  symm1 <- rbind(puboischial.lmL[[i]],puboischial.lmR[[i]])
  symm1 <- rbind(symm1,resampled.sacral.lm[[i]])
  symm1 <- symmetrize(symm1,pairedLM = cbind(1:59,(1:59)+59))
  puboischial.lmL[[i]] <- symm1[1:59,]
  puboischial.lmR[[i]] <- symm1[59+(1:59),]
  resampled.sacral.lm[[i]] <- symm1[(1:59)+118,]
  
}

names(puboischial.lmR) <- names(puboischial.lmL) <- dimnames(Mean.coord)[[3]]
print(1)

### Calculate landmark coordinates for the end of an ellipse main axis fitting between the ischio-pubic landmarks

find_A2 <- function(A1,
                    B1,
                    B2,
                    ratio){
  
  ## Midpoint of B1 and B2
  M <- (B1+B2)/2
  
  ## Minor-axis radius
  bsin <- sqrt(sum((B1-M)^2))
  
  ## Major-axis direction
  u <- M-A1
  d <- sqrt(sum(u^2))
  u <- u/d
  
  ## Solve for theta
  f <- function(theta){
    
    a <- d/(1+cos(theta))
    
    ratio*a*sin(theta)-bsin
    
  }
  
  theta <- uniroot(f,
                   interval=c(1e-6,
                              pi-1e-6))$root
  
  ## Major radius
  a <- d/(1+cos(theta))
  
  ## Compute A2
   A2 <- M + u*a*(1-cos(theta))
  
  return(A2)
  
}

## apply function to identify point on the ventral pelvis that defined the ellipse AP max axis fitting into the sacral/ischiopubic landmarks
## The point is based on fitting an ellipse with the same ratio as the neonatal face in each species

ventral.lm<-list(NULL) # Left side

for (i in 1:nlevels(List$Species)){ 
  ventral.lm[[i]]<-find_A2(resampled.sacral.lm[[i]],puboischial.lmL[[i]],puboischial.lmR[[i]],ratio=Neonate.data$Ratio.face[i]) # resample sacral landmarks
} 

names(ventral.lm) <- levels(List$Species)

dir.create("Planes_orientations")
source("./1.2.Data_Preparation.planes.r")
