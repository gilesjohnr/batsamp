r <- 30 # radius in meters
center <- SpatialPoints(coords=data.frame(x=0,y=0))
canopy <- gBuffer(center, width=r)
win <- disc(r=r, centre=c(0,0))

ntrees <- 50
treedens <- ntrees/gArea(canopy)
avgtreecan <- 3
avgclustsize <- 100

# model the cluster process
cl <- rThomas(kappa=treedens,
              scale=avgtreecan,
              mu=avgclustsize,
              win=win,
              saveLambda=T)
# estimate kernel density, which models movement around individual bats
cldens <- density(ppp(cl$x,
                      cl$y,
                      window=win),
                  kernel='gaussian',
                  sigma=0.5) # the kernel bandwith here intrinsically models movement on the individual level because a gaussian kernel with radius sigma is placed around each point, the density is the sum of all these


rast <- raster(nrow=nrow(cldens),
               ncol=ncol(cldens),
               ext=extent(-r, r,-r, r),
               crs=NA)
dists <- distanceFromPoints(rast, center)/1e5
dists <- mask(dists, canopy)


# model overall roost level movement with spatial PDF kernel centered on whole roost
gomp <- function(val, range) {
     require(flexsurv)
     p <- val/(range+1)
     gom <- dgompertz(qgompertz(p, shape=1, rate=1),
                      shape=1,
                      rate=1)
     return(gom)
}

probs <- as.data.frame(cbind(v <- getValues(dists),
                             pd <- gomp(val=v, range=r)))

probs[,3] <- probs[,2]/sum(probs[,2], na.rm=T)
rastpdf <- reclassify(dists, as.matrix(probs[-2]))
mvmt <- as.im(rastpdf)
bd <- cldens*mvmt/integral.im(cldens*mvmt)*cl$n


path <- 'figs/fig3_density.pdf'
pdf(path, width=10, height=6)

par(mfrow=c(2,3), mar=c(1,2,1,2), oma=c(0,0,0,0))
plot(cldens, main="Individual position and movement")
plot(mvmt, main="Roost level movement")
plot(bd, main="Final kernel of bat density")

persp(cldens,
      col='lightblue',
      theta=220,
      phi=30,
      zlab='density',
      main="")

persp(mvmt,
      col='lightblue',
      theta=220,
      phi=30,
      zlab='movement',
      main="")

persp(bd,
      col='lightblue',
      theta=220,
      phi=30,
      zlab='density',
      main="")

dev.off()
