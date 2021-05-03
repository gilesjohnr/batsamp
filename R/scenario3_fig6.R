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
### Simulation 3: Global sensitivity analysis using latin hyper cube sampling
#############################################################################

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

dir.create('output/scenario3', showWarnings=FALSE)

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
                                                       s=0.65,
                                                       nsheets=100,
                                                       type=c('quadrant', 'stratified'),
                                                       ntrees=lhp[i,4],
                                                       scale=lhp[i,5],
                                                       mu=lhp[i,6],
                                                       shape=lhp[i,7],
                                                       rate=lhp[i,8],
                                                       sim=i),
                                               file=paste('output/scenario3/sim_', i, '.csv', sep=""))
                             }
)

# Collate simulation output
scen3 <- bind_rows(lapply(list.files("output/scenario3",
                                    full.names=T),
                         read.csv))
scen3 <- scen3[,-1]

# Separate simulations by sheet sampling design
scen3.sq <- scen3[scen3$type=='quadrant',]
scen3.ss <- scen3[scen3$type=='stratified',]

# Just need first row of each simulation
scen3.sq <- uniqueprevs(scen3.sq, split.by='sim')
scen3.ss <- uniqueprevs(scen3.ss, split.by='sim')
scen3 <- rbind(scen3.sq, scen3.ss)
scen3$type <- factor(scen3$type)

BRTscen3 <- gbm.step(data=scen3,
                    gbm.x=c(2,10,12,15:21,23:27,29),
                    gbm.y=28,
                    family='gaussian',
                    tree.complexity=4, # allow third order interactions
                    learning.rate=0.005,
                    bag.fraction=0.7,
                    n.folds=10)

summary(BRTscen3, las=1, xlim=c(0,40))


# Plot selected fits
tmp <- as.character(factor(scen3$type))
tmp[tmp=="quadrant"] <- 'royalblue4'
tmp[tmp=="stratified"] <- 'darkorange'


path <- 'figs/fig6_scenario3.pdf'
pdf(path, width=9, height=7)

par(mfrow=c(2,3), mar=c(5,5,1,1))
# Custom relative influence plot
ri <- sort(summary(BRTscen3, plot=F)$rel.inf, decreasing=F)
vars <- as.character(rev(summary(BRTscen3, plot=F)$var))
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

plot(scen3$tot.indiv, BRTscen3$fitted,
     pch=16, col=tmp,
     ylab='Fitted estimation bias',
     xlab='Sum Cb',
     yaxt='n')
axis(side=2, las=2)
abline(h=0, lty=2, cex=2)
text('(b)', x=max(scen3$tot.indiv)*0.95, y=max(BRTscen3$fitted)*0.975, cex=1.2)

plot(scen3$trueprev, BRTscen3$fitted,
     pch=16, col=tmp,
     ylab='',
     xlab='true prevalence',
     yaxt='n')
axis(side=2, las=2)
abline(h=0, lty=2, cex=2)
text('(c)', x=max(scen3$trueprev)*0.95, y=max(BRTscen3$fitted)*0.975, cex=1.2)

plot(scen3$mI, BRTscen3$fitted,
     pch=16, col=tmp,
     ylab='Fitted estimation bias',
     xlab='Morans I',
     yaxt='n')
axis(side=2, las=2)
abline(h=0, lty=2, cex=2)
text('(d)', x=max(scen3$mI)*0.95, y=max(BRTscen3$fitted)*0.975, cex=1.2)

plot(factor(scen3$type), BRTscen3$fitted,
     pch=16, col=c('royalblue4', 'darkorange'),
     ylab='',
     xlab='Type of sampling design',
     yaxt='n')
axis(side=2, las=2)
abline(h=0, lty=2, cex=2)
text('(e)', x=2.4, y=max(BRTscen3$fitted)*0.975, cex=1.2)

plot(scen3$Nb, BRTscen3$fitted,
     pch=16, col=tmp,
     ylab='',
     xlab='Number of bat in roost (Nb)',
     yaxt='n')
axis(side=2, las=2)
abline(h=0, lty=2, cex=2)
text('(f)', x=max(scen3$Nb)*0.975, y=max(BRTscen3$fitted)*0.975, cex=1.2)


dev.off()