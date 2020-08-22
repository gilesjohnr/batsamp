# Plot sheet sample designs
par(mfrow=c(2,2), mar=c(2,2,2,2), oma = c(0, 0, 0, 0))
sheetsamp(r=30, s=1.9, nsheets=10, type='quadrant', plot=T, main="Quadrant")
sheetsamp(r=30, s=1, nsheets=100, type='uniform', plot=T, main="Uniform")
sheetsamp(r=30, s=1, nsheets=100, type='stratified', plot=T, main="Stratified")
sheetsamp(r=30, s=1, nsheets=100, type='random', plot=T, main="Random")

# Make cluster
clust <- makePSOCKcluster(8)
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

system.time(
     foreach (i=101:nsims,
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
                                         file=paste('output/sim1/sim1_', i, '.csv', sep=""))
                          }
)

# Collate simulation output
sim1 <- bind_rows(lapply(list.files("output/sim1",
                                    full.names=T),
                         read.csv))
sim1 <- sim1[,-1]

write.csv(sim1, file='output/sim1_s065_tp01.csv')

sim1 <- read.csv('output/sim1_s065_tp01.csv')
sim1 <- sim1[,-1]
str(sim1)

# Separate simulations by sheet sampling design
sq <- sim1[sim1$type=='quadrant',]
su <- sim1[sim1$type=='uniform',]
ss <- sim1[sim1$type=='stratified',]
sr <- sim1[sim1$type=='random',]

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

# Get unique estimated prevalence values
sq <- uniqueprevs(sq, split.by='sim')
su <- uniqueprevs(su, split.by='sim')
ss <- uniqueprevs(ss, split.by='sim')
sr <- uniqueprevs(sr, split.by='sim')

# Relationship between estimated prevalence and population size
par(mfrow=c(2,2))
xmn <- min(c(sq$Nb, su$Nb, ss$Nb, sr$Nb))
xmx <- max(c(sq$Nb, su$Nb, ss$Nb, sr$Nb))
ymn <- min(c(sq$estprev, su$estprev, ss$estprev, sr$estprev))
ymx <- max(c(sq$estprev, su$estprev, ss$estprev, sr$estprev))

plot(sq$Nb, sq$estprev,
     xlim=c(xmn, xmx), ylim=c(ymn, ymx),
     xlab='Nb (number of bats)',
     ylab='Estimated prevalence',
     main='Quadrant')
abline(lm(sq$estprev ~ sq$Nb), col='red')
abline(h=p, lty=2, col='grey80')

plot(su$Nb, su$estprev,
     xlim=c(xmn, xmx), ylim=c(ymn, ymx),
     xlab='Nb (number of bats)',
     ylab='Estimated prevalence',
     main='Uniform')
abline(lm(su$estprev ~ su$Nb), col='red')
abline(h=p, lty=2, col='grey80')

plot(ss$Nb, ss$estprev,
     xlim=c(xmn, xmx), ylim=c(ymn, ymx),
     xlab='Nb (number of bats)',
     ylab='Estimated prevalence',
     main='Stratified')
abline(lm(ss$estprev ~ ss$Nb), col='red')
abline(h=p, lty=2, col='grey80')

plot(sr$Nb, sr$estprev,
     xlim=c(xmn, xmx), ylim=c(ymn, ymx),
     xlab='Nb (number of bats)',
     ylab='Estimated prevalence',
     main='Random')
abline(lm(sr$estprev ~ sr$Nb), col='red')
abline(h=p, lty=2, col='grey80')

# Relationship between estimated prevalence and clustering within the roost
par(mfrow=c(2,2))
xmn <- min(c(sq$ceR, su$ceR, ss$ceR, sr$ceR))
xmx <- max(c(sq$ceR, su$ceR, ss$ceR, sr$ceR))
ymn <- min(c(sq$estprev, su$estprev, ss$estprev, sr$estprev))
ymx <- max(c(sq$estprev, su$estprev, ss$estprev, sr$estprev))

plot(sq$ceR, sq$estprev,
     xlim=c(xmn, xmx), ylim=c(ymn, ymx),
     xlab='ceR (Clark-Evans clustering coefficient)',
     ylab='Estimated prevalence',
     main='Quadrant')
abline(lm(sq$estprev ~ sq$ceR), col='red')
abline(h=p, lty=2, col='grey80')
plot(su$ceR, su$estprev, xlim=c(xmn, xmx), ylim=c(ymn, ymx),
     xlab='ceR (Clark-Evans clustering coefficient)',
     ylab='Estimated prevalence',
     main='Uniform')
abline(lm(su$estprev ~ su$ceR), col='red')
abline(h=p, lty=2, col='grey80')
plot(ss$ceR, ss$estprev, xlim=c(xmn, xmx), ylim=c(ymn, ymx),
     xlab='ceR (Clark-Evans clustering coefficient)',
     ylab='Estimated prevalence',
     main='Stratified')
abline(lm(ss$estprev ~ ss$ceR), col='red')
abline(h=p, lty=2, col='grey80')
plot(sr$ceR, sr$estprev, xlim=c(xmn, xmx), ylim=c(ymn, ymx),
     xlab='ceR (Clark-Evans clustering coefficient)',
     ylab='Estimated prevalence',
     main='Random')
abline(lm(sr$estprev ~ sr$ceR), col='red')
abline(h=p, lty=2, col='grey80')

# Relationship between spatial autocorrelation and estimated prev
# Remember spatial autocorrelation is calculated with prneg
bp <- rbind(sq[,c('type','mI')],
            su[,c('type','mI')],
            ss[,c('type','mI')],
            sr[,c('type','mI')])

par(mfrow=c(1,1))
boxplot(mI ~ type,
        data=bp,
        ylab='Morans I for probility of obtaining negative sheet',
        col='lightblue')

# Number of samples vs number of bats
par(mfrow=c(2,2))
xmn <- min(c(sq$Nb, su$Nb, ss$Nb, sr$Nb))
xmx <- max(c(sq$Nb, su$Nb, ss$Nb, sr$Nb))
ymn <- min(c(sq$nsmaps, su$nsamps, ss$nsamps, sr$nsamps))
ymx <- max(c(sq$nsmaps, su$nsamps, ss$nsamps, sr$nsamps))

plot(sq$Nb, sq$nsamps,
     xlim=c(xmn, xmx), ylim=c(ymn, ymx),
     xlab='Nb (Number of bats)', ylab='Number of samples',
     main='Quadrant')
abline(lm(sq$nsamps ~ sq$Nb), col='red')


plot(su$Nb, su$nsamps,
     xlim=c(xmn, xmx), ylim=c(ymn, ymx),
     xlab='Nb (Number of bats)', ylab='Number of samples',
     main='Uniform')
abline(lm(su$nsamps ~ su$Nb), col='red')

plot(ss$Nb, ss$nsamps,
     xlim=c(xmn, xmx), ylim=c(ymn, ymx),
     xlab='Nb (Number of bats)', ylab='Number of samples',
     main='Stratified')
abline(lm(ss$nsamps ~ ss$Nb), col='red')

plot(sr$Nb, sr$nsamps,
     xlim=c(xmn, xmx), ylim=c(ymn, ymx),
     xlab='Nb (Number of bats)', ylab='Number of samples',
     main='Random')
abline(lm(sr$nsamps ~ sr$Nb), col='red')

# Estimated prevalence vs Number of samples
par(mfrow=c(2,2))
xmn <- min(c(sq$nsmaps, su$nsamps, ss$nsamps, sr$nsamps))
xmx <- max(c(sq$nsmaps, su$nsamps, ss$nsamps, sr$nsamps))
ymn <- min(c(sq$estprev, su$estprev, ss$estprev, sr$estprev))
ymx <- max(c(sq$estprev, su$estprev, ss$estprev, sr$estprev))

plot(sq$nsamps, sq$estprev,
     xlim=c(xmn, xmx), ylim=c(ymn, ymx),
     xlab='Number of samples', ylab='Estimated prevalence',
     main='Quadrant')
abline(lm(sq$estprev ~ sq$nsamps), col='red')
abline(h=p, lty=2, col='grey80')

plot(su$nsamps, su$estprev,
     xlim=c(xmn, xmx), ylim=c(ymn, ymx),
     xlab='Number of samples', ylab='Estimated prevalence',
     main='Uniform')
abline(lm(su$estprev ~ su$nsamps), col='red')
abline(h=p, lty=2, col='grey80')

plot(ss$nsamps, ss$estprev,
     xlim=c(xmn, xmx), ylim=c(ymn, ymx),
     xlab='Number of samples', ylab='Estimated prevalence',
     main='Stratified')
abline(lm(ss$ estprev ~ ss$nsamps), col='red')
abline(h=p, lty=2, col='grey80')

plot(sr$nsamps, sr$estprev,
     xlim=c(xmn, xmx), ylim=c(ymn, ymx),
     xlab='Number of samples', ylab='Estimated prevalence',
     main='Random')
abline(lm(sr$estprev ~ sr$nsamps), col='red')
abline(h=p, lty=2, col='grey80')

#######################################################################
### Simulation 2: Simulations across multiple values of true prevalence
#######################################################################

# Same 'fixed' parameters parameters simulation 1 above
nsims <- 1000 # number of simulations

system.time(
     foreach (i=1:nsims,
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
                                         file=paste('output/sim2/sim2_', i, '.csv', sep=""))
                          }
)

# Collate simulation output
sim2 <- bind_rows(lapply(list.files("output/sim2",
                                    full.names=T),
                         read.csv))
sim2 <- sim2[,-1]
#write.csv(sim2, file='output/sim2_all.csv')

sim2 <- read.csv('output/sim2_all.csv')
sim2 <- sim2[,-1]

# Separate simulations by sheet sampling design
sq.prevs <- sim2[sim2$type=='quadrant',]
su.prevs <- sim2[sim2$type=='uniform',]
ss.prevs <- sim2[sim2$type=='stratified',]
sr.prevs <- sim2[sim2$type=='random',]

# Just need first row of each simulation
sq.prevs2 <- uniqueprevs(sq.prevs, split.by='trueprev')
su.prevs2 <- uniqueprevs(su.prevs, split.by='trueprev')
ss.prevs2 <- uniqueprevs(ss.prevs, split.by='trueprev')
sr.prevs2 <- uniqueprevs(sr.prevs, split.by='trueprev')

# Plot true vs estimated prevalance
par(mfrow=c(2,2), mar=c(4,4,2,2))
plot(sq.prevs2$trueprev, sq.prevs2$estprev,
     xlim=c(0,1), ylim=c(0,1),
     pch=16, col='royalblue4',
     ylab='estimated prevalence',
     xlab='',
     main='Quadrant')
lines(c(0,1), c(0,1), lty=2, col='red')
text(0.8, 0.1, paste('bias = ', round(sum(sq.prevs$estprev - sq.prevs$trueprev)/nrow(sq.prevs), 2)), sep="")

plot(su.prevs2$trueprev, su.prevs2$estprev,
     xlim=c(0,1), ylim=c(0,1),
     pch=16, col='royalblue4',
     ylab='',
     xlab='',
     main='Uniform')
lines(c(0,1), c(0,1), lty=2, col='red')
text(0.8, 0.1, paste('bias = ', round(sum(su.prevs$estprev - su.prevs$trueprev)/nrow(su.prevs), 2)), sep="")

plot(ss.prevs2$trueprev, ss.prevs2$estprev,
     xlim=c(0,1), ylim=c(0,1),
     pch=16, col='royalblue4',
     ylab='estimated prevalence',
     xlab='true prevalence',
     main='Stratified')
lines(c(0,1), c(0,1), lty=2, col='red')
text(0.8, 0.1, paste('bias = ', round(sum(ss.prevs$estprev - ss.prevs$trueprev)/nrow(ss.prevs), 2)), sep="")

plot(sr.prevs2$trueprev, sr.prevs2$estprev,
     xlim=c(0,1), ylim=c(0,1),
     pch=16, col='royalblue4',
     ylab='',
     xlab='true prevalence',
     main='Random')
lines(c(0,1), c(0,1), lty=2, col='red')
text(0.8, 0.1, paste('bias = ', round(sum(sr.prevs$estprev - sr.prevs$trueprev)/nrow(sr.prevs), 2)), sep="")

##############################################################################
### Simulation 3: Global sensitivity analysis using latin hyper cube sampling
#############################################################################
# Following Prowse et al 2016
nsims <- 10000

# Construct latin hypercube
lh <- randomLHS(nsims, 8)

# Empty matrix to catch parameter values
lhp <- matrix(0, nrow=nsims, ncol=8)

# populate matrix with parameter values
lhp[,1] <- qunif(lh[,1], min=25, max=50) # r = radius of roost
lhp[,2] <- qunif(lh[,2], min=0, max=1) # p = true prevalence
lhp[,3] <- qunif(lh[,3], min=0.25, max=0.75) # pu = probability urine is contributed and collected
lhp[,4] <- qunif(lh[,4], min=25, max=75) # ntrees = number of trees used in the roost
lhp[,5] <- qunif(lh[,5], min=2, max=6) # scale = radius of cluster process, or the average tree canopy size of a tree
lhp[,6] <- qunif(lh[,6], min=25, max=150) # mu = mean number of individuals in a cluster (bats in a roost tree)
lhp[,7] <- qunif(lh[,7], min=0.5, max=2) # shape of gompertz, how evenly distributed roost level movement is
lhp[,8] <- qunif(lh[,8], min=1, max=2) # rate of gompertz, how quickly movement decays towards roost edge

system.time(
     foreach (i=1:nsims,
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
                                                 s=0.65,
                                                 nsheets=100,
                                                 type=c('quadrant', 'uniform', 'stratified', 'random'),
                                                 ntrees=lhp[i,4],
                                                 scale=lhp[i,5],
                                                 mu=lhp[i,6],
                                                 shape=lhp[i,7],
                                                 rate=lhp[i,8],
                                                 sim=i),
                                         file=paste('output/sim3_run3/sim3_', i, '.csv', sep=""))
                          }
)
# Collate simulation output
sim3 <- bind_rows(lapply(list.files("output/sim3_run3",
                                    full.names=T),
                         read.csv))
str(sim3)
sim3 <- sim3[,-1]

# Separate simulations by sheet sampling design
sim3.sq <- sim3[sim3$type=='quadrant',]
#sim3.su <- sim3[sim3$type=='uniform',]
sim3.ss <- sim3[sim3$type=='stratified',]
#sim3.sr <- sim3[sim3$type=='random',]

# Just need first row of each simulation
sim3.sq <- uniqueprevs(sim3.sq, split.by='sim')
#sim3.su <- uniqueprevs(sim3.su, split.by='sim')
sim3.ss <- uniqueprevs(sim3.ss, split.by='sim')
#sim3.sr <- uniqueprevs(sim3.sr, split.by='sim')

sim3 <- rbind(sim3.sq, sim3.ss)
#sim3 <- rbind(sim3.sq, sim3.su, sim3.ss, sim3.sr)
#rm(sim3.sq, sim3.su, sim3.ss, sim3.sr)
# write datafram to file: contains unique simulations (does not include all sheets, just the first) for all four sampling designs
write.csv(sim3, file="output/sim3_alluniquesims_run3.csv")

sim3 <- read.csv('output/sim3_alluniquesims_run3.csv')
sim3 <- sim3[,-1]

BRTsim3 <- gbm.step(data=sim3,
                    gbm.x=c(2,10,12,15:21,23:27,29),
                    gbm.y=28,
                    family='gaussian',
                    tree.complexity=4, # allow third order interactions
                    learning.rate=0.005,
                    bag.fraction=0.7,
                    n.folds=10)
# save to file
saveRDS(BRTsim3, file="output/BRT_sim3_run3.rds")
BRTsim3 <- readRDS("output/BRT_sim3_run3.rds")

summary(BRTsim3, las=1, xlim=c(0,40))
par(mar=c(4,4,2,2))
gbm.plot(BRTsim3, n.plots=12)
gbm.plot.fits(BRTsim3)

# Plot selected fits
tmp <- as.character(factor(sim3$type))
tmp[tmp=="quadrant"] <- 'royalblue4'
tmp[tmp=="stratified"] <- 'darkorange'

par(mfrow=c(2,3), mar=c(5,5,1,1))
# Custom relative influence plot
ri <- sort(summary(BRTsim3, plot=F)$rel.inf, decreasing=F)
vars <- as.character(rev(summary(BRTsim3, plot=F)$var))
vars[vars=='tot.indiv'] <- 'Sum Cb'
vars[vars=='trueprev'] <- 'true prev'
vars[vars=='mI'] <- 'Morans I'
vars[vars=='s'] <- 'sheet area'
vars[vars=='r'] <- 'radius'

barplot(ri,
        names.arg=vars,
        horiz=T,
        col=rainbow(22, start=0.4, end=0.8),
        xlim=c(0,50),
        xlab='Relative influence (%)',
        las=1)
text('(a)', x=47.5, y=17, cex=1.2)

plot(sim3$tot.indiv, BRTsim3$fitted,
     pch=16, col=tmp,
     ylab='Fitted estimation bias',
     xlab='Sum Cb')
abline(h=0, lty=2, cex=2)
text('(b)', x=max(sim3$tot.indiv)*0.95, y=max(BRTsim3$fitted)*0.975, cex=1.2)

plot(sim3$trueprev, BRTsim3$fitted,
     pch=16, col=tmp,
     ylab='',
     xlab='true prevalence')
abline(h=0, lty=2, cex=2)
text('(c)', x=max(sim3$trueprev)*0.95, y=max(BRTsim3$fitted)*0.975, cex=1.2)

plot(sim3$mI, BRTsim3$fitted,
     pch=16, col=tmp,
     ylab='Fitted estimation bias',
     xlab='Morans I')
abline(h=0, lty=2, cex=2)
text('(d)', x=max(sim3$mI)*0.95, y=max(BRTsim3$fitted)*0.975, cex=1.2)

plot(factor(sim3$type), BRTsim3$fitted,
     pch=16, col=c('royalblue4', 'darkorange'),
     ylab='',
     xlab='Type of sampling design')
abline(h=0, lty=2, cex=2)
text('(e)', x=2.4, y=max(BRTsim3$fitted)*0.975, cex=1.2)

plot(sim3$Nb, BRTsim3$fitted,
     pch=16, col=tmp,
     ylab='',
     xlab='Number of bat in roost (Nb)')
abline(h=0, lty=2, cex=2)
text('(f)', x=max(sim3$Nb)*0.975, y=max(BRTsim3$fitted)*0.975, cex=1.2)

# Calculate false negative rate
round(sum(sim3.sq$falseneg)/nrow(sim3.sq), 2)
round(sum(sim3.ss$falseneg)/nrow(sim3.ss), 2)


#BRTint <- gbm.interactions(BRTsim3)
#BRTint$rank.list
#par(mfrow=c(1,2), mar=c(1,1,1,1))
#gbm.perspec(BRT1.p4, 31, 27, z.range=c(0,2), theta=30, col='lightblue')
#gbm.perspec(BRT1.p4, 40, 38, z.range=c(0,3.5), theta=120, col='lightblue')

##############################################################################
### Simulation 4: How small should sheets be and how many are needed?
### Sensitivity of the stratified design to changes in sheet size and number
#############################################################################

# Here, I run the same simulation as simulation 3 with only the stratified design. All parameter ranges/distributions stay the same except that sheet area, sheet number, and distance between sheets (previously fixed at 2m) are varied over a larger intervals.

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
system.time(
     foreach (i=1:nsims,
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
                                         file=paste('output/sim4_run3/sim4_', i, '.csv', sep=""))
                          }
)
# Collate simulation output
sim4 <- bind_rows(lapply(list.files("output/sim4_run3",
                                    full.names=T),
                         read.csv))
sim4 <- sim4[,-1]
str(sim4)

sim4 <- uniqueprevs(sim4, split.by='sim')

write.csv(sim4, file="output/sim4_alluniquesims_run2.csv")

sim4 <- read.csv('output/sim4_alluniquesims_run2.csv')
sim4 <- sim4[,-1]
#c(10,12,15:27)
c(12,16,18,20:23)
BRTsim4 <- gbm.step(data=sim4,
                    gbm.x=c(12,20:22), # type removed
                    gbm.y=28,
                    family='gaussian',
                    tree.complexity=4, # allow third order interactions
                    learning.rate=0.005,
                    bag.fraction=0.7,
                    n.folds=10)
saveRDS(BRTsim4, file="output/BRT_sim4.rds")
BRTsim4 <- readRDS("output/BRT_sim4.rds")

par(mfrow=c(1,1))
summary(BRTsim4, las=1)
par(mar=c(4,4,2,2))
gbm.plot(BRTsim4)
gbm.plot.fits(BRTsim4)

# Plot selected fits
plotBRTfit <- function(varname) {
     v <- sim4[,colnames(sim4)==varname]
     plot(v, BRTsim4$fitted,
          pch=16, col='royalblue4',
          ylab='Fitted estimation bias',
          xlab=varname)
     abline(h=0, lty=2)
     fit <- sreg(v, BRTsim4$fitted)
     x <- seq(min(v), max(v), (max(v)-min(v))*0.01)
     sp <- predict(fit, x)
     lines(x, sp, col='red', cex=1.2, type='l')
     out <- as.data.frame(cbind(x, sp))
     names(out) <- c(paste(varname), 'sreg')
     return(out)
}

par(mfrow=c(2,3), mar=c(5,4,1,1))
# Custom relative influence plot
ri <- sort(summary(BRTsim4, plot=F)$rel.inf, decreasing=F)
vars <- as.character(rev(summary(BRTsim4, plot=F)$var))
vars[vars=='tot.indiv'] <- 'Sum Cb'
vars[vars=='mI'] <- 'Morans I'

barplot(ri,
        names.arg=vars,
        horiz=T,
        col=rainbow(22, start=0.4, end=0.8),
        xlab='Relative influence (%)',
        las=1)

plotBRTfit('trueprev')
plotBRTfit('s')
plotBRTfit('Nb')

plotBRTfit('nsamps')



tmp <- tmp[tmp$sreg < 0.05,]
tmp2 <- tmp[order(tmp$sreg, decreasing=T),][1,]
points(tmp2, cex=2, col='green')
abline(v=tmp2$s, col='green')
abline(h=tmp2$sreg, col='green')


plotBRTfit('nsheets')
plotBRTfit('dsheets')

par(mfrow=c(2,4), mar=c(5,4,1,1))

##############################################################################
### Simulation 4 - BRT analysis 2: How well does a stratified design with 100 1x1m sheets perform?
# One drawback of the small-sheet design is that it less area will collect fewer
# urine samples and the total number of individuals samples will not be enough to detect low values of prevalence. Does the stratified small-sheet design consistently sample enough individuals to detect presence?
# Add minimum number of individuals that must be sampled to detect low viral prevalence from Thrusfield 2007

BRTsim4FN <- gbm.step(data=sim4,
                      gbm.x=c(12,20:22), # type removed
                      gbm.y=32,
                      family='bernoulli',
                      tree.complexity=4, # allow third order interactions
                      learning.rate=0.001,
                      bag.fraction=0.7,
                      n.folds=10)
saveRDS(BRTsim4FN, file="output/BRT_sim4run3FN.rds")
BRTsim4FN <- readRDS("output/BRT_sim4run3FN.rds")
par(mfrow=c(1,1))
summary(BRTsim4FN)
par(mar=c(4,4,2,2))
gbm.plot(BRTsim4FN)
gbm.plot.fits(BRTsim4FN)


par(mfrow=c(2,5), mar=c(5,4,1,1))
plotBRTfit <- function(varname) {
     v <- sim4[,colnames(sim4)==varname]
     plot(v, BRTsim4$fitted,
          pch=16, col='royalblue4',
          ylab='Fitted estimation bias',
          xlab=varname)
     abline(h=0, lty=2)
     fit <- sreg(v, BRTsim4$fitted)
     x <- seq(min(v), max(v), (max(v)-min(v))*0.01)
     sp <- predict(fit, x)
     lines(x, sp, col='red', cex=1.2, type='l')
     out <- as.data.frame(cbind(x, sp))
     names(out) <- c(paste(varname), 'sreg')
     return(out)
}
tmp <- plotBRTfit('s')
tmp <- tmp[order(tmp$sreg),][1,]
points(tmp, cex=2, col='green')
abline(v=tmp$s, lwd=0.5, col='green')
#abline(h=tmp$sreg, lty=2, lwd=0.5, col='green')
plotBRTfit('nsheets')
plotBRTfit('dsheets')
plotBRTfit('nsamps')
plotBRTfit('Nb')

plotBRTfit <- function(varname) {
     v <- sim4[,colnames(sim4)==varname]
     plot(v, BRTsim4FN$fitted,
          pch=16, col='royalblue4',
          ylab='Probability of false negative',
          xlab=varname)
     fit <- sreg(v, BRTsim4FN$fitted)
     x <- seq(min(v), max(v), (max(v)-min(v))*0.01)
     sp <- predict(fit, x)
     lines(x, sp, col='red', cex=1.2, type='l')
     out <- as.data.frame(cbind(x, sp))
     names(out) <- c(paste(varname), 'sreg')
     return(out)
}
plotBRTfit('s')
plotBRTfit('nsheets')
plotBRTfit('dsheets')
plotBRTfit('nsamps')
plotBRTfit('Nb')