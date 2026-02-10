library(Morpho)
library(geomorph)
library(shapes)
library(sf)
library(lattice)
library(Rvcg)
library(DescTools)
library(devtools)
library(matchingR)

#setwd("~/Desktop/Primate pelvis project/Birth canal shape and neonatal size/OUTLET/Lia script for full canal descending planes")

#Load data
Coords <- read.morphologika("1_142_Females.txt")
Mean.coord<-read.morphologika("Pelvis_mean_coords.txt")
List <- read.csv("1_IndividList.csv")
List$Species<-as.factor(List$Species)
MeanPelvesID<-read.csv("1_MeanPelvesList.csv")
MeanPelvesID$Species<-as.factor(MeanPelvesID$Species)

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

#r2morphologika(Coords_filtered, file="142_246lm_females.txt")

Coords_filtered <- procSym(Coords_filtered,CSinit=F,scale=F,pcAlign=T)$rotated

# Traditional inlets

Species.list<-levels(List$Species) #list of species

# Identifying sacral landmarks defining the posterior side of the canal

# setting up possible interlandmark distances between ends of iliopectineal line and sacrum for determining the sacral point that minimises this (most constraining inlet point)
distance.1 <- matrix(c(26,13, # angle of auricular surface and Promontorium
                       26,110,#Aur-Semilandmark sacrum 2
                       26,111,#Aur-Sml 3
                       26,112,#Aur-Sml 4
                       26,113,#Aur-Sml 5
                       26,114,#Aur-Sml 6
                       26,115,#Aur-Sml 7
                       26,116,#Aur-Sml 8
                       26,117,#Aur-Sml 9
                       26,118,#Aur-Sml 10
                       26,119,#Aur-Sml 11
                       26,120,#Aur-Sml 12  
                       26,121,#Aur-Sml 13
                       26,122,#Aur-Sml 14
                       26,15),#Aur-bottom of sacrum
                     
                     ncol = 2,byrow = T)
distance.2 <- matrix(c(28,13, #top of pubic symphysis and Promontorium
                       28,110,#pubis-Semilandmark sacrum 2
                       28,111,#pubis-Sml 3
                       28,112,#pubis-Sml 4
                       28,113,#pubis-Sml 5
                       28,114,#pubis-Sml 6
                       28,115,#pubis-Sml 7
                       28,116,#pubis-Sml 8
                       28,117,#pubis-Sml 9
                       28,118,#pubis-Sml 10
                       28,119,#pubis-Sml 11
                       28,120,#pubis-Sml 12  
                       28,121,#pubis-Sml 13
                       28,122,#pubis-Sml 14
                       28,15),#pubis-bottom of sacrum
                       ncol=2, byrow=TRUE)                      

# Resample sacral landmarks species

resampled.sacral.lm<-list(NULL)
w<-c(NULL)

for (i in 1:nlevels(List$Species)){  # loop to identify sacral landmark for most constraining inlet and build a list of most constraining inlets
    dist.1<- interlmkdist(Coords_filtered[,,List$Order[List$Species==levels(List$Species)[i]]], distance.1) # matrix of possible interlandmark distances aur-sacrum for all females in the species
    dist.2<- interlmkdist(Coords_filtered[,,List$Order[List$Species==levels(List$Species)[i]]], distance.2) # matrix of possible interlandmark distances pubis-sacrum for all females in the species
    sum.dist<-dist.1+dist.2 #sum of two distances by individual
    w[i]<-distance.1[which.min(apply(sum.dist,2,mean)),2] # identify landmark to use for sacrum for most constraining inlet 
    
}

for (i in 1:nlevels(List$Species)){ # loop to resample sacral landmarks so that we have 15
sacral.lm<-c(w[i]:122,15)
resampled.sacral.lm[[i]]<-resampleCurve(Mean.coord[sacral.lm,,MeanPelvesID$Species==levels(List$Species)[i]],15) # matrix of possible interlandmark distances pubis-sacrum for all females in the species
sum.dist<-dist.1+dist.2 #sum of two distances by individual
    if(levels(List$Species) != "Hsapiens")
        w[i]<-distance.1[which.min(apply(sum.dist,2,mean)),2] # identify landmark to use for sac)],15)
    else
        w[i] <- 13
resampled.sacral.lm[[i]]<- rbind(resampled.sacral.lm[[i]],resampled.sacral.lm[[i]][rep(15,each=14),]) # repeats the last sacral landmark another 14 time, to match lm on pubo-ischial curve.
}

names(resampled.sacral.lm) <- levels(List$Species)
##dimnames(Mean.coord)[[3]]

## puboischial landmarks on each side

puboischial.lmL<-list(NULL) # Left side
puboischial.lmR<-list(NULL) # Right side
for (i in 1:nlevels(List$Species)){
    j <- grep(levels(List$Species)[[i]],dimnames(Mean.coord)[[3]])
  puboischial.lmL[[i]] <- Mean.coord[c(28,123:135,30,175:187,42),,j]
  puboischial.lmR[[i]] <- Mean.coord[c(29,136:148,31,188:200,43),,j]
  symm1 <- rbind(puboischial.lmL[[i]],puboischial.lmR[[i]])
  symm1 <- rbind(symm1,resampled.sacral.lm[[i]])
  symm1 <- symmetrize(symm1,pairedLM = cbind(1:29,(1:29)+29))
  puboischial.lmL[[i]] <- symm1[1:29,]
  puboischial.lmR[[i]] <- symm1[29+(1:29),]
  resampled.sacral.lm[[i]] <- symm1[(1:29)+58,]
  
}

names(puboischial.lmR) <- names(puboischial.lmL) <- dimnames(Mean.coord)[[3]]
print(1)
