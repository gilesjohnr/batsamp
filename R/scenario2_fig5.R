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

#######################################################################
### Simulation 2: Simulations across multiple values of true prevalence
#######################################################################

# Same 'fixed' parameters parameters scenario 1
r <- 30
pu <- 0.5
s <- 0.65
nsheets <- 100
ntrees <- 50
scale <- 3
mu <- 100
shape <- 0.8
rate <- 1
nsims <- 1000 # number of simulations

dir.create('output/scenario2', showWarnings=FALSE)

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

                               p <- runif(1, min=0, max=1) # true prevalence allowed to vary randomly

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
                                         file=paste('output/scenario2/sim_', i, '.csv', sep=""))
                          }
)

# Collate simulation output
scen2 <- bind_rows(lapply(list.files("output/scenario2",
                                    full.names=T),
                         read.csv))
scen2 <- scen2[,-1]

# Separate simulations by sheet sampling design
sq.prevs <- scen2[scen2$type=='quadrant',]
su.prevs <- scen2[scen2$type=='uniform',]
ss.prevs <- scen2[scen2$type=='stratified',]
sr.prevs <- scen2[scen2$type=='random',]

# Just need first row of each simulation
sq.prevs2 <- uniqueprevs(sq.prevs, split.by='trueprev')
su.prevs2 <- uniqueprevs(su.prevs, split.by='trueprev')
ss.prevs2 <- uniqueprevs(ss.prevs, split.by='trueprev')
sr.prevs2 <- uniqueprevs(sr.prevs, split.by='trueprev')



path <- 'figs/fig5_scenario2.pdf'
pdf(path, width=6, height=6)

# Plot true vs estimated prevalance
par(mfrow=c(2,2), mar=c(4,4,2,2))
plot(sq.prevs2$trueprev, sq.prevs2$estprev,
     xlim=c(0,1), ylim=c(0,1),
     pch=1, col='royalblue4',
     ylab='estimated prevalence',
     xlab='',
     main='Quadrant')
lines(c(0,1), c(0,1), lty=2, lwd=2, col='red')
text(0.8, 0.1, paste0('bias = ', round(sum(sq.prevs$estprev - sq.prevs$trueprev)/nrow(sq.prevs), 2)))

plot(su.prevs2$trueprev, su.prevs2$estprev,
     xlim=c(0,1), ylim=c(0,1),
     pch=1, col='royalblue4',
     ylab='',
     xlab='',
     main='Uniform')
lines(c(0,1), c(0,1), lty=2, lwd=2, col='red')
text(0.8, 0.1, paste('bias = ', round(sum(su.prevs$estprev - su.prevs$trueprev)/nrow(su.prevs), 2)), sep="")

plot(ss.prevs2$trueprev, ss.prevs2$estprev,
     xlim=c(0,1), ylim=c(0,1),
     pch=1, col='royalblue4',
     ylab='estimated prevalence',
     xlab='true prevalence',
     main='Stratified')
lines(c(0,1), c(0,1), lty=2, lwd=2, col='red')
text(0.8, 0.1, paste('bias = ', round(sum(ss.prevs$estprev - ss.prevs$trueprev)/nrow(ss.prevs), 2)), sep="")

plot(sr.prevs2$trueprev, sr.prevs2$estprev,
     xlim=c(0,1), ylim=c(0,1),
     pch=1, col='royalblue4',
     ylab='',
     xlab='true prevalence',
     main='Random')
lines(c(0,1), c(0,1), lty=2, lwd=2, col='red')
text(0.8, 0.1, paste('bias = ', round(sum(sr.prevs$estprev - sr.prevs$trueprev)/nrow(sr.prevs), 2)), sep="")

dev.off()