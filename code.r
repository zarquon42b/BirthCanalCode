require(Morpho)
require(Rvcg)
require(mesheR)
require(geomorph)
require(rgl)
source("./functions.r")
source("./1.Data_Preparation. 3pointplanes.r")
lms <- read.morphologika("./3D models pelvis/Pelvis_mean_coords.txt")
lmnames <- dimnames(lms)[[3]]

## load and decimate meshes (only do once!)
meshfiles <- list.files("./3D models pelvis",full.names=T,pattern=".ply")
dir.create("./3D models pelvis/DecMeshes/")

for (i in meshfiles) {
   
        tmpmesh <- vcgImport(i)
        tmpmesh <- vcgQEdecim(tmpmesh, tarface=100000)
    vcgPlyWrite(tmpmesh,paste0("./3D models pelvis/DecMeshes/",basename(i)))
}



meshfilesDec <- list.files("./3D models pelvis/DecMeshes/",full.names = T,pattern="*.ply")
lmn <- dimnames(lms)[[3]]
lmnspec <- name2factor(lmn,which=2,as.factor = F)
mn <- gsub(".ply","",basename(meshfilesDec))

### run scan with auto-generated ellipses
scanlist <- parallel::mclapply(1:length(meshfilesDec), workhorse, mc.cores = 8)

## save to disk
names(scanlist) <- mn
saveRDS(scanlist,"scanlist.RDS")
scanlist <- readRDS("./scanlist.RDS")
scanlist1.1 <- parallel::mclapply(1:length(meshfilesDec), workhorse,ellscale=1.1, mc.cores = 8)
names(scanlist1.1) <- mn
saveRDS(scanlist1.1,"scanlist1.1.RDS")



## Create Ellipses for each slice
require(DescTools)
dir.create("Ellipses",showWarnings=FALSE)

for (k in c("scanlist","scanlist1.1")) {
    arealist <- ellipselist <- list()
    myscans <- get(k) ## replace scanlist with scanlist1.1 for plotting with 10% increased cylinder
    dir.create(paste0("Ellipses/",k,"/"),showWarnings = FALSE)
    for (use in 1:length(myscans)) {
       
        specname <- paste0("Ellipses/",k,"/",lmn[[use]])
        mydims <- vcgKDtree(lms[,,use],lms[,,use],k=32)
        mymax <- max(mydims$distance)*3
        ielistRot <- list()
        myscan <- myscans[[use]]
        mycol <- rainbow(nrow(myscan$midpoints))
        dir.create(specname,showWarnings = FALSE)
        for (i in 1:nrow(myscan$midpoints)) {
            
            ca <- computeArea(myscan$result[[i]])
            ## generate trafo matrix from 3D to 2D
            trafo <- computeTransform(cbind(ca$xpro2D,0),ca$xpro3D)
            sacrum2ca <- applyTransform(myscan$sacrum[i,],trafo)[,1:2]
            pubic2ca <- applyTransform(myscan$midpub[i,],trafo)[,1:2]
            ## align to sacrum/midpubic axis
            myangle <- cangle(c(0,10),sacrum2ca-pubic2ca)
            if (myangle > pi/2)
                myangle <- -(pi-myangle)
            ## generate rotation matrix
            myrot <- create2Drot(myangle)

            ca$xpro2D <- applyTransform(ca$xpro2D,myrot)
            ie <- inscribeEllipse(ca$xpro2D,iters = 2000,maxratio = 1.5)
            png(paste0(specname,"/",sprintf("%02d.png",i)),width = 1000,height = 1000)
            plot(ca$xpro2D,asp=1,xlim=c(-mymax,mymax),ylim=c(-mymax,mymax))
            lines(rbind(ca$xpro2D,ca$xpro2D[1,]))
            DrawEllipse(x=ie$center[1],y=ie$center[2],radius.x = ie$radius.x,radius.y = ie$radius.y,col=mycol[i])
            points(applyTransform(rbind(sacrum2ca,pubic2ca),myrot),col=2,pch=19)
            lines(c(0,0),c(-100,100))
            dev.off()
            ie$trafo1 <- trafo
        ie$trafo2 <- myrot
            ielistRot[[i]] <- ie
            
        }
        areas <- sapply(ielistRot,function(x) x <- x$maxarea)
        area.f <- data.frame(slice=1:nrow(myscan$sacrum),area=areas)
        arealist[[use]] <- area.f
        ellipselist[[use]] <-ielistRot 
        which.min(areas)
        openxlsx::write.xlsx(area.f,paste0(specname,"/",lmn[[use]],"_areaoutlet.xlsx"))
    }
    names(arealist) <- names(ellipselist) <- names(scanlist)
    saveRDS(arealist,"arealist.RDS")
    saveRDS(ellipselist,"ellipselist.RDS")
    ellipselist <- readRDS("ellipselist.RDS")
    centerlist <- list()
    for(i in 1:length(ellipselist))
        centerlist[[i]] <- t(sapply(ellipselist[[i]], function(x) x <- applyTransform(cbind(applyTransform(x$center,x$trafo2,T),0),x$trafo1,T)))
    names(centerlist) <- names(scanlist)
    saveRDS(centerlist,"centerlist.RDS")
    
    ellipsedimlist <- list()
    for(i in 1:length(ellipselist))
        ellipsedimlist[[i]] <- t(sapply(ellipselist[[i]], function(x) x <- cbind(x$radius.x,x$radius.y)))
    
    for (i in 1:length(ellipselist)) {
        temp <- cbind(centerlist[[i]],ellipsedimlist[[i]])
        colnames(temp) <- c("center.x","center.y","center.z","radius.x","radius.y")
        specname <- paste0("Ellipses/",k,"/",names(ellipselist)[[i]])
        openxlsx::write.xlsx(temp,paste0(specname,"/",names(ellipselist)[[i]],"_ellipses.xlsx"))
    }
}

### END OF PRODUCTION CODE
## area.m <- (sapply(arealist,function(x) x <- x$area))
## area.m <- scale((area.m),scale=T)
## plot(x=1:29,y=area.m[,1],ylim=range(area.m),cex=0,ylab="standardized area")
## lines(x=1:29,y=area.m[,1],lwd=2)
## for (i in 2:32)
##     lines(x=1:29,y=area.m[,i],col=i,lwd=2)

## area.mx <- cbind((area.m),path=1:29)
## require(reshape2)
## areamelt <- melt(area.mx,id.vars = "path")
## areamelt <- areamelt[-which(areamelt$Var2=="path"),]
## plot(areamelt$Var1,areamelt$value)
## mylm <- lm(areamelt$value ~ poly(areamelt$Var1,7))


## car::sp(areamelt$value ~ areamelt$Var1,ylab="normalized available space",xlab="slice number",regLine=F)
## library(ggplot2)
## ggplot(areamelt, aes(x = Var1, y = value)) +
##   geom_point() +
##   stat_smooth(method = "lm", formula = y ~ poly(x,7), level = 0.999999999)

