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


hev_indiv <- read.csv('./data/hev_indiv_data.csv', stringsAsFactors=F)
hev_roost <- read.csv('./data/hev_roost_data.csv', stringsAsFactors=F)

#------------------------------------------------------------
# Individual-level prevalence
#------------------------------------------------------------

hev_indiv <- hev_indiv[hev_indiv$Species == 'BFF',]
hev_indiv$Date <- as.Date(hev_indiv$Date, format='%m/%d/%y')
hev_indiv$yearmon <- format(hev_indiv$Date, '%Y-%m')
hev_indiv$mon <- format(hev_indiv$Date, '%m')
hev_indiv <- hev_indiv[order(hev_indiv$Date),]

hev_indiv$res <- NA
for (i in 1:nrow(hev_indiv)) hev_indiv$res[i] <- as.integer(any(hev_indiv[i,c(19:25)] != '>40'))
hev_indiv$res[is.na(hev_indiv$res)] <- 0

tmp <- split(hev_indiv$res, factor(hev_indiv$Date))
tmp <- data.frame(date=as.Date(names(tmp)),
                  pos=unlist(lapply(tmp, sum)),
                  n=unlist(lapply(tmp, length)),
                  row.names=NULL)
tmp$p <- tmp$pos/tmp$n

mod <- mgcv::gam(res ~ s(as.integer(Date), bs='tp', k=7), data=hev_indiv, family='quasibinomial')
summary(mod)

date_seq <- min(hev_indiv$Date):max(hev_indiv$Date)
preds <- predict(mod, newdata=data.frame(Date=date_seq), type='response', se.fit=T)
df <- data.frame(date=as.Date(date_seq, origin="1970-01-01"), Sample='Individual', as.data.frame(preds))
df_indiv <- merge(df, tmp, all=T)


#------------------------------------------------------------
# Pooled prevalence
#------------------------------------------------------------

hev_roost <- hev_roost[hev_roost$loc == 'Boonah',]
hev_roost$date <- as.Date(hev_roost$date, format='%d-%b-%y')
hev_roost <- hev_roost[hev_roost$date >= as.Date('2013-05-02'),]

date_seq <- min(hev_roost$date):max(hev_roost$date)

tmp <- split(hev_roost$hev, factor(hev_roost$date))
tmp <- data.frame(date=as.Date(names(tmp)),
                  pos=unlist(lapply(tmp, sum)),
                  n=unlist(lapply(tmp, length)),
                  row.names=NULL)
tmp$p <- tmp$pos/tmp$n

mod <- mgcv::gam(hev ~ s(as.integer(date), bs='tp', k=8), data=hev_roost, family='quasibinomial')
summary(mod)
preds <- predict(mod, newdata=data.frame(date=date_seq), type='response', se.fit=T)
df <- data.frame(date=as.Date(date_seq, origin="1970-01-01"), Sample='Pooled quadrant', as.data.frame(preds))
df_roost <- merge(df, tmp, all=T)

#------------------------------------------------------------
# Sheet prevalence
#------------------------------------------------------------

hev_sheet <- hev_roost %>%
     group_by(date, sheet) %>%
     group_modify(~ data.frame(hev=as.integer(any(.x$hev == 1))))

tmp <- hev_sheet %>%
     group_by(date) %>%
     group_modify(~ data.frame(pos=sum(.x$hev),
                               n=length(.x$hev),
                               p=sum(.x$hev)/length(.x$hev)))


mod <- mgcv::gam(hev ~ s(as.integer(date), bs='tp', k=8), data=hev_sheet, family='quasibinomial')
summary(mod)

preds <- predict(mod, newdata=data.frame(date=date_seq), type='response', se.fit=T)
df <- data.frame(date=as.Date(date_seq, origin="1970-01-01"), Sample='Pooled sheet', as.data.frame(preds))
df_sheet <- merge(df, tmp, all=T)

# Combine all models
df <- rbind(df_indiv, df_roost, df_sheet)



#---------------------------------------------------------------------------------
# Validation of theoretical model using estimated viral prevalence from GAM models
#---------------------------------------------------------------------------------

# Here all parameters are allowed to vary, but true prevalence is bounded by the values of viral prevalence that were estimated
# by the GAM models gitted to observed field data (Boonah, QLD, June 2013--June 2014)


nsims <- 1000

lh <- randomLHS(nsims, 11) # Construct latin hypercube
lhp <- matrix(0, nrow=nsims, ncol=11) # transform hypercube to parameter values

# populate matrix with parameter values
lhp[,1] <- qunif(lh[,1], min=25, max=50) # r = radius of roost

tmp <- MASS::fitdistr(df_indiv$p[!is.na(df_indiv$p)] + 1e-3, 'beta', start=list(shape1=2, shape2=2)) # fit Beta distribution to observed individual-level prevalence data
lhp[,2] <- qbeta(lh[,2], shape1=tmp$estimate[1], shape2=tmp$estimate[2]) # p = true prevalence

lhp[,3] <- qunif(lh[,3], min=0.25, max=0.75) # pu = probability urine is contributed and collected
lhp[,4] <- qunif(lh[,4], min=0.625, max=1.75) # s = sheet area for tratified design. (0.25 sq m up to 2 sq m area)
lhp[,5] <- qunif(lh[,5], min=25, max=150) # nsheets = number of sheets placed for uniform, stratified, and random
lhp[,6] <- qunif(lh[,6], min=25, max=75) # ntrees = number of trees used in the roost
lhp[,7] <- qunif(lh[,7], min=2, max=6) # scale = radius of cluster process, or the average tree canopy size of a tree
lhp[,8] <- qunif(lh[,8], min=25, max=150) # mu = mean number of individuals in a cluster (bats in a roost tree)
lhp[,9] <- qunif(lh[,9], min=0.5, max=2) # shape of gompertz, how evenly distributed roost level movement is
lhp[,10] <- qunif(lh[,10], min=1, max=2) # rate of gompertz, how quickly movement decays towards roost edge
lhp[,11] <- qunif(lh[,11], min=0.5, max=5) # dsheets = distance between sheets

dir.create('output/scenario5', showWarnings=FALSE)

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
                                                 type='quadrant',
                                                 ntrees=lhp[i,6],
                                                 scale=lhp[i,7],
                                                 mu=lhp[i,8],
                                                 shape=lhp[i,9],
                                                 rate=lhp[i,10],
                                                 d=0.1,
                                                 sim=i),
                                         file=paste('output/scenario5/sim_', i, '.csv', sep=""))
                          }
)

# Collate simulation output
scen5 <- bind_rows(lapply(list.files("output/scenario5",
                                    full.names=T),
                         read.csv))

head(scen5)

sim_sheet <- scen5 %>%
     group_by(sim, sheetID) %>%
     group_modify(~ data.frame(trueprev=.x$trueprev,
                               poolprev=.x$estprev,
                               nsheets=sum(.x$puCb) >= 1,                  # sheet provided a pooled urine sample
                               pos=as.integer(any(.x$sheetbin == 1)))) %>%  # sheet was positive
     group_by(sim) %>%
     group_modify(~ data.frame(trueprev=.x$trueprev,
                               poolprev=.x$poolprev,
                               sheetprev=sum(.x$pos)/sum(.x$nsheets)))

sim_sheet <- uniqueprevs(sim_sheet, split.by='sim')

sim$bias1 <- sim$poolprev - sim$trueprev
sim$bias2 <- sim$sheetprev - sim$trueprev
sim$mag1 <- sim$poolprev/sim$trueprev
sim$mag2 <- sim$sheetprev/sim$trueprev
