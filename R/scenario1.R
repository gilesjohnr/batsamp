source('R/load.R')
source('R/funcs.R')

# Make cluster
clust <- makeCluster(4, setup_timeout=0.5)
registerDoParallel(clust) # register as foreach back end
getDoParRegistered() # check
getDoParWorkers()

# Export functions to cluster
clusterExport(clust, c("foreach", "sheetsamp", "raster", "extent", "distanceFromPoints", "mask",
                       "getValues", "reclassify", "qgompertz", "dgompertz", "gomp", "Moran.I", "minsampsize"))

######################################################################################################
### Simulation 1: sheet sampling designs applied to simulated bat density, true prevalence is fixed
######################################################################################################

# Set 'fixed' parameters for both simulation 1 and 2
r <- 30
p <- 0.1
pu <- 0.5
s <- 0.65
nsheets <- 100
ntrees <- 50
scale <- 3
mu <- 100
shape <- 0.8
rate <- 1
nsims <- 1000

dir.create('output/scenario1', showWarnings=FALSE)

system.time(
     foreach (i=1:nsims,
              .errorhandling="remove",
              .packages=c('sp',
                          'spatstat',
                          'stringr',
                          'plotrix',
                          'rgeos',
                          'maptools',
                          'flexsurv')) %dopar% {

                               write.csv(prevsim(r=r,
                                                 p=p,
                                                 pu=pu,
                                                 s=s,
                                                 nsheets=nsheets,
                                                 type=c('quadrant', 'uniform', 'stratified', 'random'),
                                                 ntrees=ntrees,
                                                 scale=scale,
                                                 mu=mu,
                                                 shape=shape,
                                                 rate=rate,
                                                 sim=i),
                                         file=paste('output/scenario1/sim', i, '.csv', sep=""))
                          }
)

# Collate simulation output
scen1 <- bind_rows(lapply(list.files("output/scenario1",
                                    full.names=T),
                         read.csv))
scen1 <- scen1[,-1]

# Separate simulations by sheet sampling design
sq <- scen1[scen1$type=='quadrant',]
su <- scen1[scen1$type=='uniform',]
ss <- scen1[scen1$type=='stratified',]
sr <- scen1[scen1$type=='random',]

# Plot mean estimated prevalence and prob of neg sheet for each sampling design
par(mfrow=c(2,4), mar=c(5,5,2,2), oma = c(0, 0, 0, 0))
plotsim.estp(sq, main='Quadrant')
plotsim.estp(su, main='Uniform')
plotsim.estp(ss, main='Stratified')
plotsim.estp(sr, main='Random')
plotsim.prneg(sq, maxCb=max(sq$Cb))
plotsim.prneg(su, maxCb=max(sq$Cb))
plotsim.prneg(ss, maxCb=max(sq$Cb))
plotsim.prneg(sr, maxCb=max(sq$Cb))
