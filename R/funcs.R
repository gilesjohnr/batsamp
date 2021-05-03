batdens1D <- function(x, mean=0, sd=2, lower=-r, upper=r) {
     numer <- dnorm(x, mean, sd)
     denom <- pnorm(upper, mean, sd) - pnorm(lower, mean, sd)
     x <- (numer/denom)
     return(x)
}

##############################################################################
sheetsamp <- function(r=30, # radius of roost
                      s=1.2, # sheet area
                      nsheets=100, # number of sheets
                      dsheets=2, # distance between sheets (stratified only)
                      type='stratified', # sampling design
                      plot=TRUE,
                      main=""
) {
     require(sp)
     require(spatstat)
     require(rgeos)
     require(plotrix)
        require(raster)

     canopy <- gBuffer(SpatialPoints(coords=data.frame(x=0,y=0)), width=r)
     hexpts <- spsample(canopy, type="hexagonal", cellsize=s)
     hexpol <- HexPoints2SpatialPolygons(hexpts)

     if (type=='quadrant') {

          nsheets <- 10
          pts <- rSSI(r=1.9*3, # inhibitory radius (dist between large sheet centers)
                      n=nsheets, # number of points = number of sheets
                      win=disc(r=r-s, centre=c(0,0)))

          dists <- gDistance(SpatialPoints(cbind(pts$x, pts$y)), hexpts, byid=T)
          m <- matrix(nrow=0, ncol=2) # the sheet that each quadrant belongs to

          for (sheet in 1:nsheets) {
               quad <- as.numeric(names(sort(dists[,sheet])[1:4]))
               m <- rbind(m, cbind(sheet, quad))
          }

          sheets <- hexpol[as.vector(m[,2])]
          coords <- coordinates(sheets)
          sheets <- SpatialPolygonsDataFrame(sheets,
                                             data=data.frame(sheet=m[,1],
                                                             quad=m[,2],
                                                             x=coords[,1],
                                                             y=coords[,2],
                                                             row.names=row.names(sheets)))
     }

     if (type=='random') {

          sheets <- hexpol[as.integer(sample(1:length(hexpol),
                                             size=nsheets,
                                             replace=F))]
          coords <- coordinates(sheets)
          sheets <- SpatialPolygonsDataFrame(sheets,
                                             data=data.frame(sheet=str_sub(names(sheets), start=3),
                                                             x=coords[,1],
                                                             y=coords[,2],
                                                             row.names=row.names(sheets)))

     }

     if (type=='uniform') {

          pts <- spsample(canopy, n=nsheets, type='regular')
          dists <- gDistance(pts, hexpts, byid=T)
          v <- vector()

          for (i in 1:ncol(dists)) {
               # get names of hexagons nearest point pattern
               v <- c(v, as.numeric(names(sort(dists[,i])[1])))
          }

          sheets <- hexpol[v] # subset hexagonal grid to give sheet pattern
          coords <- coordinates(sheets)
          sheets <- SpatialPolygonsDataFrame(sheets,
                                             data=data.frame(sheet=str_sub(names(sheets), start=3),
                                                             x=coords[,1],
                                                             y=coords[,2],
                                                             row.names=row.names(sheets)))
     }

     if (type=='stratified') {

          pts <- rSSI(r=dsheets, # distance between sheets
                      n=nsheets, # number of sheets
                      win=disc(r=r, centre=c(0,0)))
          dists <- gDistance(SpatialPoints(cbind(pts$x, pts$y)), hexpts, byid=T)
          v <- vector()

          for (i in 1:ncol(dists)) {
               # get names of hexagons nearest point pattern
               v <- c(v, as.numeric(names(sort(dists[,i])[1])))
          }

          sheets <- hexpol[v] # subset hexagonal grid to give sheet pattern
          coords <- coordinates(sheets)
          sheets <- SpatialPolygonsDataFrame(sheets,
                                             data=data.frame(sheet=str_sub(names(sheets), start=3),
                                                             x=coords[,1],
                                                             y=coords[,2],
                                                             row.names=row.names(sheets)))
     }

     crs(sheets) <- crs("+proj=longlat +datum=WGS84")

     if (plot == TRUE) {

          par(mar=c(3,3,3,3))
          plot(hexpol, col='grey90', border='white', main=main)
          plot(sheets, col='red', border='white', add=T)

          for (i in seq(5,r,5)) {draw.circle(0,0, radius=i, lwd=0.8, border='grey70')}
          segments(x0=-r-(r*0.1), y0=0, x1=r+(r*0.1), y1=0, col='grey40')
          segments(x0=0, y0=-r-(r*0.1), x1=0, y1=r+(r*0.1), col='grey40')
          text(x=seq(-r-0.6, r+0.4, 5),
               y=-1.4,
               seq(-r, r, 5),
               cex=0.65,
               col='grey40')
     }

     return(sheets)
}

##############################################################################

hexarea <- function(r) {r^2*6 * sin(60*(pi/180)) / 2} # area of sampling sheet = hexagon area, radius of hexagonal tiles is half the cellsize

##############################################################################

batdens2D <- function(r = 60, # radius of roost in meters
                    nbats=5000,
                    ntrees = 20, # mean number of trees used in a roost
                    treedist = 5, # minimum distance between roosting trees
                    terrsize = 2, # size of male territories in meters
                    haremsize = 8, # mean size of harem. Fem + juv
                    mfactor = 300, # factor to inflate density of males
                    simwin = 30, # cluster process MUST be simulated in a smaller window. This is why the initial roost radius is 60m. simqin is subtracted from the initial r value, So that the final simulation window will be the desired radius of 30m
                    plot=T
) {
     # set simulation area
     win <- disc(r=r, centre=c(0,0)) # object of class owin
     roostarea <- gBuffer(SpatialPoints(coords=data.frame(x=0, y=0)), width=r)

     # locations of roost trees
     trees <- rSSI(r=treedist, # inhibitory radius = expected distance between two mature trees large enough to be used for roosting
                   n=rnorm(1, mean=ntrees, sd=ntrees*0.25), # average number of roosting trees
                   win=disc(r=r-simwin-treedist, centre=c(0,0))) # trees at least 'treedist' meters from roost edge. Roost edge defined by simulation window

     treedens <- density(ppp(trees$x, trees$y, window=disc(r=r, centre=c(0,0))),
                         sigma=3,
                         diggle=T,
                         positive=T)

     # male positions
     males <- rMaternII(kappa=as.im(treedens)*mfactor, # tree density used as intensity function for male positions. Because there are so few trees, the kernel density function needs to be inflated to give the desired intensity of simulated male positions
                        r=terrsize, # inhibitory radius between males (territory size)
                        win=win)

     # harem size and distribution: an inhomogenous Poisson process with male territories supplied as parent point process
     harem <- rMatClust(kappa=density(ppp(males$x, males$y, window=disc(r=r, centre=c(0,0))), #
                                      sigma=1, # bandwidth
                                      diggle=T, # edge correction
                                      positive=T), # force positive density values
                        scale=terrsize, # radius of clusters (territory size)
                        mu=haremsize, # mean number of individuals in cluster
                        win=disc(r=r-simwin, centre=c(0,0)), # reduced simulation window
                        saveLambda=T) # save parent points
     males <- attributes(harem)$parents

     # kernel function of harem density
     haremdens <- density(ppp(harem$x,
                              harem$y,
                              window=disc(r=r-simwin, centre=c(0,0))),
                          sigma=1,
                          diggle=T,
                          positive=T)

     # movement within roost modeled as kernel density of roost trees with large bandwidth
     mvmt <- density(ppp(trees$x, trees$y, window=disc(r=r-simwin, centre=c(0,0))),
                     sigma=8,
                     diggle=T,
                     positive=T)

     # multiply mvmt and harmdens to get final kernel density
     kern <- haremdens*mvmt/integral.im(haremdens*mvmt) # insure combined kernel densities integrate to 1
     kern <- kern*nbats
     return(kern)

     if (plot==T) {
          par(mfrow=c(1,2))
          # Point positions of simulated trees, male territories, and associated harems
          plot(harem, pch=2, cols='grey40', cex=1.25)
          points(males$x, males$y, pch=15, cex=0.8, col='blue')
          points(trees, pch=10, col='red', cex=3)
          # Kernel representing movement among roost trees
          plot(mvmt)
          points(trees, pch=10, cex=1.5)

          par(mfrow=c(2,3))
          # top row shows the density of individual roosting locations, movement kernel, and the combined kernel estimator including both. Bottom row shows the same thing just as perspective plots
          plot(haremdens)
          plot(mvmt)
          plot(kern)
          persp(haremdens,
                col='lightblue',
                theta=30,
                phi=30,
                zlab='Roosting density',
                main="")
          persp(mvmt,
                col='lightblue',
                theta=30,
                phi=30,
                zlab='Movement density',
                main="")
          persp(kern,
                col='lightblue',
                theta=30,
                phi=30,
                zlab='Final bat density',
                main="")

          par(mfrow=c(1,2))
          # simulation of sheet samlping
          sheets <- sheetsamp(r=r-simwin, s=1.9, nsheets=10, type='quadrant', plot=T)
          # point locations of roosting individuals
          points(harem, pch=1, cex=0.8)
          # density of the final combined kernel function along with simulated sheet placement
          plot(kern, main="")
          plot(sheets, border='green', add=T)
     }
}

##############################################################################
# Spatial PDF of Gompertz distribution
gomp <- function(val, range, shape, rate) {
     require(flexsurv)
     p <- val/(range+1)
     gom <- dgompertz(qgompertz(p,
                                shape=shape,
                                rate=rate),
                      shape=1, rate=1)
     return(gom)
}

##################################################################################################
##################################################################################################
# non-parallel version of prevsim
# this version creates one scenario of bat density and performs X number of simulations for all 4 sheet sampling designs
# put into parallel foreach loop to generate different scenarios

prevsim <- function(r=30, # radius of roost
                    p=0.1, # true prevalence
                    pu=0.5, # probability urine contributed & collected
                    s=1.9, # sheet area
                    nsheets=100, # number of sheets (total sheets = x4 for quadrant design)
                    dsheets=2, # distance between sheets (stratified only)
                    type='quadrant', # random or quadrant based design
                    ntrees=30, # number of trees, tree density = kappa
                    scale=1.5, # radius around parent points
                    mu=11, # mean number of cluster pts within radius
                    shape=0.6, # shape parameter for movement model (Gompertz pdf)
                    rate=1, # rate parameter for movement model (Gompertz pdf)
                    d=0.05, # desired precision when estimating prevalence (for minimum individual sample size claculation)
                    sim=1 # unique identifier for the simulation
) {

     require(raster)
     require(stringr)
     require(sp)
     require(rgeos)
     require(maptools)
     require(spatstat)
     require(doParallel)
     require(plotrix)
     require(ape)
     require(flexsurv)

     # set simulation window
     center <- SpatialPoints(coords=data.frame(x=0,y=0))
     canopy <- gBuffer(center, width=r)
     win <- disc(r=r, centre=c(0,0))

     ##################################
     ### Build scenario of bat density
     #################################
     # model the individual cluster process

     treedens <- ntrees/gArea(canopy) # mean tree density = kappa
     cl <- rThomas(kappa=treedens, # parent points (trees) are random homogeneous Poisson process
                   scale=scale, # radius around parent points (tree canopies)
                   mu=mu, # mean cluster size = bats in tree, drawn from Poisson with mean mu
                   win=win,
                   saveLambda=T)

     cldens <- density(ppp(cl$x,
                           cl$y,
                           window=win),
                       kernel='gaussian',
                       sigma=0.5) # sigma models indiviual movement with a gaussian kernel with a radius of 0.5 meters

     # model roost-wide movement with a spatial Gompertz pdf
     # make raster with cell values of distance to center
     rast <- raster(nrow=nrow(cldens),
                    ncol=ncol(cldens),
                    ext=extent(-r, r,-r, r),
                    crs=NA)
     dists <- suppressWarnings(distanceFromPoints(rast, center)/1e5)
     crs(dists) <- CRS("+proj=longlat")
     dists <- mask(dists, canopy)

     # reclassify raster values with the Gompertz PDF values
     probs <- as.data.frame(cbind(v <- getValues(dists),
                                  pd <- gomp(val=v,
                                             range=r,
                                             shape=shape,
                                             rate=rate)))
     probs[,3] <- probs[,2]/sum(probs[,2], na.rm=T)
     rastpdf <- reclassify(dists, as.matrix(probs[-2]))
     mvmt <- as.im(rastpdf)

     # ensure B(x) integrates to one and multiply by Nb
     bd <- cldens*mvmt/integral.im(cldens*mvmt)*cl$n

     ####################################
     ### Generate sheet sampling designs
     ###################################

     msim <- foreach(t=type, .combine='rbind') %do% {

          if (t=='quadrant') {

               sheets <- sheetsamp(r=r, s=1.9, nsheets=10, type=t, plot=F)

          } else {

               sheets <- sheetsamp(r=r, s=s, nsheets=nsheets,
                                   dsheets=dsheets, type=t, plot=F)
               }

          m <- matrix(nrow=0, ncol=7)
          for (i in 1:length(sheets)) {

               Cb <- integral(bd, domain=as.owin(SpatialPolygons(sheets@polygons[i])))
               puCb <- as.integer(round(Cb*rnorm(1, pu, 0.1)))
               prneg <- round((1-pu)^(puCb*p), digits=2)
               sheetbin <- as.integer(sum(rbinom(n=puCb, size=1, prob=p)) >= 1)
               coords <- coordinates(sheets[i,])
               m <- rbind(m,
                          c(sheets$sheet[i],
                            coords[1],
                            coords[2],
                            Cb,
                            puCb,
                            prneg,
                            sheetbin))
          }

          nsamps <- sum(as.numeric(m[,5]) > 0.5) # number of sheets returning samples
          estp <- round(sum(as.integer(m[,7]))/nsamps, 3) # estimated prevalence
          upper.se <- round(estp+(1/sqrt(nsamps)), 3) # upper limit of sampling error
          lower.se <- round(estp-(1/sqrt(nsamps)), 3) # lower limit of sampling error
          ceR <- round(clarkevans(cl, correction='cdf'), 3) # crude measure of clustering
          Nb <- cl$n # simulated roost population size
          # Calculate Morans I for prneg between sheets
          dst <- as.matrix(dist(cbind(as.numeric(as.character(m[,2])),
                                      as.numeric(as.character(m[,3])))))
          dst <- 1/dst # inverse distance matrix
          diag(dst) <- 0 # diagonal set to zero
          mI <- round(Moran.I(as.numeric(as.character(m[,6])), dst)$observed, 3)
          bias <- round(estp - p, 3)
          tot.indiv <- as.integer(sum(as.numeric(m[,5]))) # total number individuals sampled
          # minimum individual sample size required to detect presence
          reqsampsize <- as.integer(minsampsize(p=p,
                                                d=d,
                                                N=Nb))
          minsamp <-  as.integer(tot.indiv >= reqsampsize)
          falseneg <- as.integer(estp == 0 && p > 0)

          if (t == 'quadrant') {
               cbind(sim, t, m, p, estp, nsamps, upper.se, lower.se, ceR, Nb,
                     mI, r, pu, 1.9, 10, dsheets, ntrees, scale, mu, shape,
                     rate, bias, tot.indiv, reqsampsize, minsamp, falseneg)
          } else {
               cbind(sim, t, m, p, estp, nsamps, upper.se, lower.se, ceR, Nb,
                     mI, r, pu, s, nsheets, dsheets, ntrees, scale, mu, shape,
                     rate, bias, tot.indiv, reqsampsize, minsamp, falseneg)
          }
     }

     # clean up matrix
     return(data.frame(sim=as.integer(as.character(msim[,1])),
                       type=as.character(msim[,2]),
                       sheetID=as.integer(as.character(msim[,3])),
                       x=as.numeric(as.character(msim[,4])),
                       y=as.numeric(as.character(msim[,5])),
                       Cb=as.integer(as.character(msim[,6])),
                       puCb=as.integer(as.character(msim[,7])),
                       prneg=as.numeric(as.character(msim[,8])),
                       sheetbin=as.integer(as.character(msim[,9])),
                       trueprev=as.numeric(as.character(msim[,10])),
                       estprev=as.numeric(as.character(msim[,11])),
                       nsamps=as.integer(as.character(msim[,12])),
                       upper.se=as.numeric(as.character(msim[,13])),
                       lower.se=as.numeric(as.character(msim[,14])),
                       ceR=as.numeric(as.character(msim[,15])),
                       Nb=as.integer(as.character(msim[,16])),
                       mI=as.numeric(as.character(msim[,17])),
                       r=as.numeric(as.character(msim[,18])),
                       pu=as.numeric(as.character(msim[,19])),
                       s=as.numeric(as.character(msim[,20])),
                       nsheets=as.integer(as.character(msim[,21])),
                       dsheets=as.numeric(as.character(msim[,22])),
                       ntrees=as.integer(as.character(msim[,23])),
                       scale=as.numeric(as.character(msim[,24])),
                       mu=as.numeric(as.character(msim[,25])),
                       shape=as.numeric(as.character(msim[,26])),
                       rate=as.numeric(as.character(msim[,27])),
                       bias=as.numeric(as.character(msim[,28])),
                       tot.indiv=as.integer(as.character(msim[,29])),
                       reqsampsize=as.integer(as.character(msim[,30])),
                       minsamp=as.integer(as.character(msim[,31])),
                       falseneg=as.integer(as.character(msim[,32]))))
}

############################################################################
# Function to extract the first row of a simulation. This removes all but the measurements for the first sheet. Giving a shorter dataframe containing unique values of metrics calculated for the entire roost (i.e. trueprev, estprev, Nb, ceR)
uniqueprevs <- function (x, split.by='sim') {
     x2 <- x[1,]
     if (split.by=='trueprev') {xsplit <- split(x, x$trueprev)}
     if (split.by=='sim') {xsplit <- split(x, x$sim)}
     for(i in 2:length(xsplit)) {
          x2 <- rbind(x2, xsplit[[i]][1,])
     }
     return(x2)
}

############################################################################
# Function to make quick plots of simulations
plotsim <- function(x, title) {
     # get vectors of estimated true and estimated prevalences for each simulation
     estp <- summarize(group_by(x, sim), estp=estprev[1])

     par(mfrow=c(1,2), mar=c(5,5,2,2), oma = c(0, 0, 2, 0))
     plot(estp$sim, estp$estp, ylim=c(0,1),
          xlab='Simulation',
          ylab='Estimated prevalence')
     abline(h=x$trueprev[1], lty=2, col='red') # true prevalence
     abline(h=mean(estp$estp), col='blue', cex=1.25) # mean estimated prevalence
     abline(h=mean(x$upper.se), col='lightblue', cex=1.25) # mean upper bound of sampling error
     abline(h=mean(x$lower.se),  col='lightblue', cex=1.25) # mean lower bound of sampling error

     plot(x$Cb, x$prneg,
          ylim=c(0,1),
          xlab='Bat density (Cb)',
          ylab='Probability of obtaining a negative sheet')
     abline(h=1-x$trueprev[1], col='green')
     mtext(paste(title), outer = TRUE, cex = 1.5)
}

plotsim.estp <- function(x, main) {
     estp <- summarize(group_by(x, sim), estp=estprev[1])
     plot(estp$sim, estp$estp, ylim=c(0,1),
          pch=16, col='grey45',
          xlab='Simulation',
          ylab='Estimated prevalence',
          main=main)
     abline(h=x$trueprev[1], lty=2, col='red') # true prevalence
     abline(h=mean(estp$estp), col='blue', cex=2) # mean estimated prevalence
     abline(h=mean(x$upper.se), lty=3, col='blue', cex=2) # mean upper bound of sampling error
     x$lower.se[x$lower.se < 0] <- 0
     abline(h=mean(x$lower.se), lty=3, col='blue', cex=2) # mean lower bound of sampling error
}

plotsim.prneg <- function(x, maxCb) {
     estp <- summarize(group_by(x, sim), estp=estprev[1])
     plot(x$Cb, x$prneg,
          xlim=c(0,maxCb),
          ylim=c(0,1),
          pch=16, col='grey45',
          xlab='Bat density (Cb)',
          ylab='Probability of obtaining a negative sheet',
          main="")
     abline(h=1-x$trueprev[1], col='green', cex=2)
}

######################################################################################
minsampsize <- function(p=0.05,
                        d=0.05,
                        N=5000) {
     n <- 1.96^2*p*(1-p) / d^2
     nadj <- (N*n) / (N+n)
     return(round(nadj))
}

colnum <- function(x) as.data.frame(colnames(x))
