loadmystuff()
require(geomorph)
require(rgl)
source("./functions.r")
source("./1.Data_Preparation. 3pointplanes.r")
lms <- read.morphologika("./3D models pelvis/Pelvis_mean_coords.txt")

lmnames <- dimnames(lms)[[3]]
dir.create("./3D models pelvis/DecMeshes/")
for (i in meshfiles) {
    #species <- name2factor(i,which=2,as.factor = F)
    #speccyl <- grep(species,cyls)
    #if (length(speccyl)) {
        tmpmesh <- vcgImport(i)
        tmpmesh <- vcgQEdecim(tmpmesh, tarface=100000)
        
        vcgPlyWrite(tmpmesh,paste0("./3D models pelvis/DecMeshes/",basename(i)))
        
    #}
    
}
##file.exists(paste0("./3D models pelvis/",dimnames(lms)[[3]],".ply"))
##sphere <- vcgSphere()
##mydists <- rep(0,dim(lms)[3])
##names(mydists) <- dimnames(lms)[[3]]
meshfilesDec <- list.files("./3D models pelvis/DecMeshes/",full.names = T,pattern="*.ply")
lmn <- dimnames(lms)[[3]]
lmnspec <- name2factor(lmn,which=2,as.factor = F)
mn <- gsub(".ply","",basename(meshfilesDec))

workhorse <- function(use) {
    message(paste0("processing: ",lmn[use]))
    m.use <- grep(lmn[use],mn)
    x.m <- vcgImport(meshfilesDec[[m.use]])
    
    coord_use <- grep(lmnspec[use],names(resampled.sacral.lm))
    sacrum <- resampled.sacral.lm[[coord_use]]
    pubicL <- puboischial.lmL[[coord_use]]
    pubicR <- puboischial.lmR[[coord_use]]
    
    
    midpoints <- getMidline(sacrum,pubicL,pubicR) ## get path through birthcanal
    cyldir <- sacrum[29,]-sacrum[1,] ## get cylinder dir by direction of sacrum
    ##cyldir <- midpoints[29,]-midpoints[1,]
    p2p <- points2plane(lms[c(13, 15:19, 24:43, 52:246),,use],midpoints[1,],cyldir)
    cA <- computeArea(p2p)
    ## create cylinder and merge with pelvis
    ellicyl <- elliCyl(cA,cyldir,scale=10)
    merged2.m <- mergeMeshes(x.m,ellicyl)
    
    scan1 <- scanitNew(resampled.sacral.lm[[coord_use]],puboischial.lmL[[coord_use]],puboischial.lmR[[coord_use]],merged2.m)
    scan1$mesh <- merged2.m
    scan1$sacrum <- sacrum
    scan1$midpub <- (pubicL+pubicR)/2
    return(scan1)
}

scanlist <- parallel::mclapply(1:length(meshfilesDec), workhorse, mc.cores = 8)
##homoscan <- workhorse(1)
names(scanlist) <- mn
saveRDS(scanlist,"scanlist.RDS")
scanlist <- readRDS("./scanlist.RDS")
## Visualization of scans
clear3d()
myscan <- scanlist[[5]]
wire3d(myscan$mesh,col="white")
spheres3d(myscan$midpoints)

mycol <- rainbow(nrow(sacrum))
for (i in 1:length(myscan$result)) {
    spheres3d(myscan$result[[i]],col=mycol[i])
    pt=points2plane(myscan$midpoints[i,],normal=myscan$dirs[[i]],v1=c(0,0,0))
    d=sqrt(sum((pt-myscan$midpoints[i,])^2))
    difv <- pt-myscan$midpoints[i,]
    if(crossprod(difv,myscan$dirs[[i]]) < 0)
        d <- -d
    planes3d(myscan$dirs[[i]][1],myscan$dirs[[i]][2],myscan$dirs[[i]][3],d=d,alpha=.3,col=mycol[i])
    
}

## Ellipses
require(DescTools)
rot=F
dir.create("Ellipses")
use <- 5
##scanlist <- lapply(c(1,5),workhorse)
##scanlist[[5]] <- scanlist[[2]]
arealist <- list()
for (use in 1:length(scanlist)) {
    specname <- paste0("Ellipses/",lmn[[use]])
    
    mydims <- vcgKDtree(lms[,,use],lms[,,use],k=32)
    mymax <- max(mydims$distance)*3
    ielistRot <- list()
    myscan <- scanlist[[use]]
        
    for (i in 1:nrow(myscan$midpoints)) {
        
        ca <- computeArea(myscan$result[[i]])
        trafo <- computeTransform(cbind(ca$xpro2D,0),ca$xpro3D)
        sacrum2ca <- applyTransform(myscan$sacrum[i,],trafo)[,1:2]
        pubic2ca <- applyTransform(myscan$midpub[i,],trafo)[,1:2]
        myangle <- cangle(c(0,10),sacrum2ca-pubic2ca)
        if (myangle > pi/2)
            myangle <- -(pi-myangle)
        myrot <- create2Drot(myangle)
        ca$xpro2D <- applyTransform(ca$xpro2D,myrot)
        if (rot)
            ie <- inscribeEllipseRot(ca$xpro2D,iters = 1000,maxpi=pi/8,rotsteps = 16,threads=8)
        else
            ie <- inscribeEllipse(ca$xpro2D,iters = 2000,maxratio = 1.5)
        
        if (rot) {
            specnamerot <- paste0(specname,"_rot")
            dir.create(specnamerot,showWarnings=FALSE)
            png(paste0(specnamerot,"/",sprintf("%02d_rot.png",i)),width = 1000,height = 1000)
            plot(ie$polyRot,asp=1,xlim=c(-mymax,mymax),ylim=c(-mymax,mymax))
            
            lines(rbind(ie$polyRot,ie$polyRot[1,]))
        }
        else {
            dir.create(specname,showWarnings = FALSE)
            png(paste0(specname,"/",sprintf("%02d.png",i)),width = 1000,height = 1000)
            plot(ca$xpro2D,asp=1,xlim=c(-mymax,mymax),ylim=c(-mymax,mymax))
            
            lines(rbind(ca$xpro2D,ca$xpro2D[1,]))
        }
        DrawEllipse(x=ie$center[1],y=ie$center[2],radius.x = ie$radius.x,radius.y = ie$radius.y,col=mycol[i])
        points(applyTransform(rbind(sacrum2ca,pubic2ca),myrot),col=2,pch=19)
        lines(c(0,0),c(-100,100))
        dev.off()
        ielistRot[[i]] <- ie
        
    }
    areas <- sapply(ielistRot,function(x) x <- x$maxarea)
    area.f <- data.frame(slice=1:nrow(sacrum),area=areas)
    arealist[[use]] <- area.f
    which.min(areas)
    if (rot)
        openxlsx::write.xlsx(area.f,paste0(specnamerot,"/",lmn[[use]],"_areaoutlet.xlsx"))
    else
        openxlsx::write.xlsx(area.f,paste0(specname,"/",lmn[[use]],"_areaoutlet.xlsx"))
}
names(arealist) <- names(scanlist)
saveRDS(arealist,"arealist.RDS")
area.m <- (sapply(arealist,function(x) x <- x$area))
area.m <- scale((area.m),scale=T)
plot(x=1:29,y=area.m[,1],ylim=range(area.m),cex=0,ylab="standardized area")
lines(x=1:29,y=area.m[,1],lwd=2)
for (i in 2:32)
    lines(x=1:29,y=area.m[,i],col=i,lwd=2)

mylm <- lm(t(area.m) ~ matrix(1:29,32,29,byrow = T))
car::sp(t(area.m),1:29)
area.mx <- cbind((area.m),path=1:29)
require(reshape2)
areamelt <- melt(area.mx,id.vars = "path")
areamelt <- areamelt[-which(areamelt$Var2=="path"),]
plot(areamelt$Var1,areamelt$value)
mylm <- lm(areamelt$value ~ poly(areamelt$Var1,7))


car::sp(areamelt$value ~ areamelt$Var1,ylab="normalized available space",xlab="slice number")
library(ggplot2)
ggplot(areamelt, aes(x = Var1, y = value)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ poly(x,7), level = 0.95)

require(betareg)
mybeta <- betareg(areamelt$value ~ areamelt$Var1)
### EOC

mesh <- a$mesh
split=100
filename=lmnames[1]
folder="./3D models pelvis/"


for (i in 1:dim(lms)[3]) {
    mylm <- lms[,,i]
    myfile <- paste0("./3D models pelvis/",dimnames(lms)[[3]][i],".ply")
    mymesh <- vcgImport(myfile)
    mymesh <- vcgQEdecim(mymesh,edgeLength = 1)
    mean1 <- colMeans(mylm[22:23,])
    mean2a <- colMeans(mylm[30:31,])
    mean2b <- mylm[15,]
    mean2 <- colMeans(rbind(mean2a,mean2b))
    myseq <- seq(from=.3,to=.7,length.out = 25)
    dir <- mean2b-mean2a
    dist <- 0
    cent <- 0
    for(j in myseq) {
        tmpcent <- mean2a+j*dir
        spheres3d(tmpcent)
        clost <- vcgClostKD(t(tmpcent),mymesh,sign=F,threads = 0)
        if (clost$quality > dist){
            dist <- clost$quality
            myclost <- clost
            cent <- tmpcent
        }
    }
    clostV <- vert2points(myclost)
    ##spheres3d(vert2points(myclost),col=4)
    mysphereS <- scalemesh(sphere,size = myclost$quality)
    mysphereT <- translate3d(mysphereS,cent[1],cent[2],cent[3])
    ##    return(list(sphere=mysphereT, clost=myclost))
    vcgPlyWrite(mysphereT,gsub(".ply","_sphere.ply",myfile))
    mydists[i] <- myclost$quality
}
openxlsx::write.xlsx(data.frame(radius=mydists),"dists.xlsx",colNames=T,rowNames=T)

