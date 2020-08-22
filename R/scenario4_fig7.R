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


##############################################################################
### Scenario 4: How small should sheets be and how many are needed?
### Sensitivity of the stratified design to changes in sheet size and number
#############################################################################

# Here, I run the same scenario as scenario 3 with only the stratified design.
# All parameter ranges/distributions stay the same except that sheet area, sheet number,
# and distance between sheets (previously fixed at 2m) are varied over a larger intervals.

nsims <- 10000

# Construct latin hypercube
lh <- randomLHS(nsims, 11)

# transform hypercube to parameter values
lhp <- matrix(0, nrow=nsims, ncol=11)

# populate matrix with parameter values
lhp[,1] <- qunif(lh[,1], min=25, max=50) # r = radius of roost
lhp[,2] <- qunif(lh[,2], min=0, max=1) # p = true prevalence
lhp[,3] <- qunif(lh[,3], min=0.25, max=0.75) # pu = probability urine is contributed and collected
lhp[,4] <- qunif(lh[,4], min=0.625, max=1.75) # s = sheet area for tratified design. (0.25 sq m up to 2 sq m area)
lhp[,5] <- qunif(lh[,5], min=25, max=150) # nsheets = number of sheets placed for uniform, stratified, and random
lhp[,6] <- qunif(lh[,6], min=25, max=75) # ntrees = number of trees used in the roost
lhp[,7] <- qunif(lh[,7], min=2, max=6) # scale = radius of cluster process, or the average tree canopy size of a tree
lhp[,8] <- qunif(lh[,8], min=25, max=150) # mu = mean number of individuals in a cluster (bats in a roost tree)
lhp[,9] <- qunif(lh[,9], min=0.5, max=2) # shape of gompertz, how evenly distributed roost level movement is
lhp[,10] <- qunif(lh[,10], min=1, max=2) # rate of gompertz, how quickly movement decays towards roost edge
lhp[,11] <- qunif(lh[,11], min=0.5, max=5) # dsheets = distance between sheets

dir.create('output/scenario4', showWarnings=FALSE)

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

                                     write.csv(prevsim(r=lhp[i,1],
                                                       p=lhp[i,2],
                                                       pu=lhp[i,3],
                                                       s=lhp[i,4],
                                                       nsheets=lhp[i,5],
                                                       dsheets=lhp[i,11],
                                                       type='stratified',
                                                       ntrees=lhp[i,6],
                                                       scale=lhp[i,7],
                                                       mu=lhp[i,8],
                                                       shape=lhp[i,9],
                                                       rate=lhp[i,10],
                                                       d=0.1,
                                                       sim=i),
                                               file=paste('output/scenario4/sim_', i, '.csv', sep=""))
                             }
)

# Collate simulation output
scen4 <- bind_rows(lapply(list.files("output/scenario4",
                                    full.names=T),
                         read.csv))

scen4 <- scen4[,-1]
scen4 <- uniqueprevs(scen4, split.by='sim')

# Estimation bias
BRTscen4 <- gbm.step(data=scen4,
                     gbm.x=c(12,20:22), # type removed
                     gbm.y=28,
                     family='gaussian',
                     tree.complexity=4, # allow third order interactions
                     learning.rate=0.005,
                     bag.fraction=0.7,
                     n.folds=10)

# Prob of false negative
BRTscen4FN <- gbm.step(data=scen4,
                       gbm.x=c(12,20:22), # type removed
                       gbm.y=32,
                       family='bernoulli',
                       tree.complexity=4, # allow third order interactions
                       learning.rate=0.001,
                       bag.fraction=0.7,
                       n.folds=10)



path <- 'figs/fig7_scenario4.pdf'
pdf(path, width=9, height=7)

par(mfrow=c(2,4), mar=c(4,4,1,1))

### s
varname <- 's'
v <- scen4[,colnames(scen4)==varname]
plot(v, BRTscen4$fitted,
     pch=16, col='royalblue4',
     ylab='Fitted estimation bias',
     xlab='',
     xaxt='n')
axis(1, at=seq(0.625,1.75,0.25), labels=c(0.25, 0.5, 0.75, 1.25, 1.75))
abline(h=0, lty=2)
fit <- sreg(v, BRTscen4$fitted, lambda=1e-5)
x <- seq(min(v), max(v), (max(v)-min(v))*0.01)
sp <- predict(fit, x)
lines(x, sp, col='red', cex=1.2, type='l')
tmp <- as.data.frame(cbind(x, sp))
names(tmp) <- c(paste(varname), 'sreg')
text('(a)', x=max(v)*0.4, y=max(BRTscen4$fitted)*0.98, cex=1.2)

### nsheets
varname <- 'nsheets'
v <- scen4[,colnames(scen4)==varname]
plot(v, BRTscen4$fitted,
     pch=16, col='royalblue4',
     ylab='',
     xlab='')
abline(h=0, lty=2)
fit <- sreg(v, BRTscen4$fitted, lambda=1e4)
x <- seq(min(v), max(v), (max(v)-min(v))*0.01)
sp <- predict(fit, x)
lines(x, sp, col='red', cex=1.2, type='l')
tmp <- as.data.frame(cbind(x, sp))
names(tmp) <- c(paste(varname), 'sreg')
text('(b)', x=max(v)*0.2, y=max(BRTscen4$fitted)*0.98, cex=1.2)

### dsheets
varname <- 'dsheets'
v <- scen4[,colnames(scen4)==varname]
plot(v, BRTscen4$fitted,
     pch=16, col='royalblue4',
     ylab='',
     xlab='')
abline(h=0, lty=2)
fit <- sreg(v, BRTscen4$fitted, lambda=1e-3)
x <- seq(min(v), max(v), (max(v)-min(v))*0.01)
sp <- predict(fit, x)
lines(x, sp, col='red', cex=1.2, type='l')
tmp <- as.data.frame(cbind(x, sp))
names(tmp) <- c(paste(varname), 'sreg')
text('(c)', x=max(v)*0.15, y=max(BRTscen4$fitted)*0.98, cex=1.2)

### namps
varname <- 'nsamps'
v <- scen4[,colnames(scen4)==varname]
plot(v, BRTscen4$fitted,
     pch=16, col='royalblue4',
     ylab='',
     xlab='')
abline(h=0, lty=2)
fit <- sreg(v, BRTscen4$fitted, lambda=1e3)
x <- seq(min(v), max(v), (max(v)-min(v))*0.01)
sp <- predict(fit, x)
lines(x, sp, col='red', cex=1.2, type='l')
tmp <- as.data.frame(cbind(x, sp))
names(tmp) <- c(paste(varname), 'sreg')
text('(d)', x=max(v)*0.05, y=max(BRTscen4$fitted)*0.98, cex=1.2)

###############
# Prob of false negative

### s
varname <- 's'
v <- scen4[,colnames(scen4)==varname]
plot(v, BRTscen4FN$fitted,
     pch=16, col='royalblue4',
     ylab='Probability of false negative',
     xlab=expression(paste('Sheet area ', (m^2))),
     xaxt='n')
#round(hexarea((seq(0.625,1.75,0.25)/2)), 2)
axis(1, at=seq(0.625,1.75,0.25), labels=c(0.25, 0.5, 0.75, 1.25, 1.75))

fit <- sreg(v, BRTscen4FN$fitted, lambda=1e-4)
x <- seq(min(v), max(v), (max(v)-min(v))*0.01)
sp <- predict(fit, x)
lines(x, sp, col='red', cex=1.2, type='l')
tmp <- as.data.frame(cbind(x, sp))
names(tmp) <- c(paste(varname), 'sreg')
text('(e)', x=max(v)*0.95, y=max(BRTscen4FN$fitted)*0.98, cex=1.2)

### nsheets
varname <- 'nsheets'
v <- scen4[,colnames(scen4)==varname]
plot(v, BRTscen4FN$fitted,
     pch=16, col='royalblue4',
     ylab='',
     xlab='Number of sheets')
fit <- sreg(v, BRTscen4FN$fitted, lambda=1e4)
x <- seq(min(v), max(v), (max(v)-min(v))*0.01)
sp <- predict(fit, x)
lines(x, sp, col='red', cex=1.2, type='l')
tmp <- as.data.frame(cbind(x, sp))
names(tmp) <- c(paste(varname), 'sreg')
text('(f)', x=max(v)*0.95, y=max(BRTscen4FN$fitted)*0.98, cex=1.2)

### dsheets
varname <- 'dsheets'
v <- scen4[,colnames(scen4)==varname]
plot(v, BRTscen4FN$fitted,
     pch=16, col='royalblue4',
     ylab='',
     xlab='Distance between sheets (m)')
fit <- sreg(v, BRTscen4FN$fitted, lambda=1e-3)
x <- seq(min(v), max(v), (max(v)-min(v))*0.01)
sp <- predict(fit, x)
lines(x, sp, col='red', cex=1.2, type='l')
tmp <- as.data.frame(cbind(x, sp))
names(tmp) <- c(paste(varname), 'sreg')
text('(g)', x=max(v)*0.95, y=max(BRTscen4FN$fitted)*0.98, cex=1.2)

### nsamps
varname <- 'nsamps'
v <- scen4[,colnames(scen4)==varname]
plot(v, BRTscen4FN$fitted,
     pch=16, col='royalblue4',
     ylab='',
     xlab='Number of samples')
fit <- sreg(v, BRTscen4FN$fitted, lambda=1e4)
x <- seq(min(v), max(v), (max(v)-min(v))*0.01)
sp <- predict(fit, x)
lines(x, sp, col='red', cex=1.2, type='l')
tmp <- as.data.frame(cbind(x, sp))
names(tmp) <- c(paste(varname), 'sreg')
text('(h)', x=max(v)*0.95, y=max(BRTscen4FN$fitted)*0.98, cex=1.2)

dev.off()

