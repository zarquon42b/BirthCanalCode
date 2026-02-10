scanitNew <- function(sacrum,pubicL,pubicR,mesh,split=100) {
    result <- dirs <- list()
    midpoints <- NULL

    for (i in 1:nrow(sacrum)) {
        mp <- colMeans(rbind(colMeans(rbind(pubicL[i,],pubicR[i,])),sacrum[i,]))
        midpoints <- rbind(midpoints,mp)
        
       # if (i != nrow(x))
        dir <- crossProduct(sacrum[i,]-pubicL[i,],pubicR[i,]-pubicL[i,])
       ## else
       ##     dir <- x[i,]-x[i-1,]
        
       ## dir <- x[nrow(x),]-x[1,]
        
        dirs[[i]] <- dir
        tP <- tangentPlane(dir)
        startDir <- mp+tP$y
        degrees <- seq(from=0,to=2*pi,length.out=split)
        out <- t(sapply(1:split,function(y) x <- rotaxis3d(t(startDir),mp,mp+dir,theta = degrees[y])))
        out <- t(out)-mp
        tmpmesh <- list(vb=t(matrix(mp,split,3,byrow = T)))
        tmpmesh$normals <- (out)
        class(tmpmesh) <- "mesh3d"
        b <- vcgRaySearch(tmpmesh,mesh)
        bad <- which(b$quality == 0)
        hits <- vert2points(b)
        if (length(bad))
            hits <- hits[-bad,]
            
        result[[i]] <- hits
    }
    return(list(result=result,midpoints=midpoints,dirs=dirs))
}

getMidline <- function(sacrum,pubicL,pubicR) {
    ## midpointInd <- cbind(c(28,123:135,30),c(29,136:148,31))
    ## midpoint1 <- t(sapply(1:nrow(midpointInd),function(y) colMeans(x[midpointInd[y,],])))
    ## ind2 <- c(13,110:122,15)
    ## mp2 <- x[ind2,]
    ## midpoint2 <- (mp2+midpoint1)/2
    ## return(list(midpoints=midpoint2,midpoint1=midpoint1,midpoint2=mp2))
    
    midpoints <- NULL
    for (i in 1:nrow(sacrum)) {
        mp <- colMeans(rbind(colMeans(rbind(pubicL[i,],pubicR[i,])),sacrum[i,]))
        midpoints <- rbind(midpoints,mp)
    }
    return(midpoints)
    
    
}


getIt <- function(lc=100) {
    it <- NULL
    for (i in 1:(lc - 1)) {
        face0 <- c(i, i + 1, i + lc)
        face1 <- c(i + 1, i + lc + 1, i + lc)
        it <- rbind(it, face0, face1)
    }
    it <- rbind(it, c(lc, 1, lc + 1))
    it <- rbind(it, c(2 * lc, lc, lc + 1))
   
    return(it)
}

require(shotGroups)
elliCyl <- function(x,dir,scale=1,ellscale=1) {
    cAtrafo <- computeTransform(x$xpro3D,cbind(x$xpro2D,0))
    me <- getMinEllipse(x$xpro2D)
    me$cov <- me$cov*ellscale
    png()
    plot(x$xpro2D)
    de <- drawEllipse(me, fg='blue')
    dev.off()
    de3d <- applyTransform(cbind(de,0),cAtrafo)
    cyldir <- dir*scale
    de13d <- translate3d(de3d,-cyldir[1],-cyldir[2],-cyldir[3])
    de23d <- translate3d(de3d,cyldir[1],cyldir[2],cyldir[3])
    ellcyl <- list(vb=t(cbind(rbind(de13d,de23d),1)))
    ellcyl$it <- t(getIt(100))
    class(ellcyl) <- "mesh3d"
    return(ellcyl)
}

checkMeshes <- function(lms,meshfilesDec,start=NULL) {
    mn <- gsub(".ply","",basename(meshfilesDec))
    lmn <- dimnames(lms)[[3]]
    lmnspec <- name2factor(lmn,which=2,as.factor = F)
    if(is.null(start))
        i <- 1
    else
        i <- start
    shitlist <- NULL
    while (i <= dim(lms)[3]) {
        use <- i
        m.use <- grep(lmn[use],mn)
        if(length(m.use)) {
        x.m <- vcgImport(meshfilesDec[[m.use]])
        clear3d()
        wire3d(x.m)
        spheres3d(lms[,,use],col=3)
        answer <- readline(paste("viewing #", i, "(return=next | f=false)\n"))
        if (answer=="f") {
            shitlist <- c(shitlist,lmn[use])
            i <- i+1
            }
        else
            i <- i+1
        } else {
            print(lmn[use])
            i <- i+1
        }
        
    }
    return(shitlist)
}

cangle <- function(x,y) {
    dot <- crossprod(x,y)
    det1 <- det(cbind(x,y))
    angle <- atan2(det1,dot)
    return(angle)
}

create2Drot <- function(angle) {
    diag1 <- cos(angle)
    diag2 <- sin(angle)
    out <- matrix(c(diag1,diag2,-diag2,diag1),2,2)
    out <- Morpho:::mat2homg(rbind(out,0))
    return(out)
}

visualizePath <- function(myscan,planes=T) {
    if (rgl.cur()!= 0)
        clear3d()
    wire3d(myscan$meshOrig,col="white")
    spheres3d(myscan$midpoints)
    
    mycol <- rainbow(nrow(myscan$midpoints))
    for (i in 1:length(myscan$result)) {
        spheres3d(myscan$result[[i]][myscan$hits[[i]],],col=mycol[i])
        spheres3d(myscan$result[[i]][-myscan$hits[[i]],],col=mycol[i],radius=.2)
        if (planes) {
            pt=points2plane(myscan$midpoints[i,],normal=myscan$dirs[[i]],v1=c(0,0,0))
            d=sqrt(sum((pt-myscan$midpoints[i,])^2))
            difv <- pt-myscan$midpoints[i,]
            if(crossprod(difv,myscan$dirs[[i]]) < 0)
                d <- -d
            planes3d(myscan$dirs[[i]][1],myscan$dirs[[i]][2],myscan$dirs[[i]][3],d=d,alpha=.3,col=mycol[i])
        }
    }
}


workhorse <- function(use,ellscale=1) {
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
    ellicyl <- elliCyl(cA,cyldir,scale=10,ellscale=ellscale)
    
    merged2.m <- mergeMeshes(x.m,ellicyl)
    
    
    scan1 <- scanitNew(resampled.sacral.lm[[coord_use]],puboischial.lmL[[coord_use]],puboischial.lmR[[coord_use]],merged2.m)
    scan2 <- scanitNew(resampled.sacral.lm[[coord_use]],puboischial.lmL[[coord_use]],puboischial.lmR[[coord_use]],x.m)
    for (j in 1:length(scan1$result)) {
        clost <- vcgKDtree(scan2$result[[j]],scan1$result[[j]],k=1)
        good <- which(clost$distance < 1e-2)
        scan1$hits[[j]] <- good
    }

    scan1$meshOrig <- x.m
    scan1$mesh <- merged2.m
    scan1$sacrum <- sacrum
    scan1$midpub <- (pubicL+pubicR)/2
    return(scan1)
}
