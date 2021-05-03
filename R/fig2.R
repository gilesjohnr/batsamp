library(ggplot2)
library(dplyr)
library(lemon)

#------------------------------------------------------------
# Data
#------------------------------------------------------------

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




#------------------------------------------------------------
# Plot GAM estimates of viral prevalence
#------------------------------------------------------------

# All models together
p0 <- ggplot(df) +
     geom_ribbon(aes(x=date, ymin=fit - se.fit, ymax=fit + se.fit, fill=Sample), alpha=0.4) +
     geom_line(aes(x=date, y=fit, color=Sample), size=1.5) +
     ylab('HeV Prevalence') + xlab('Date') +
     ylim(min(df$fit - df$se.fit), max(df$fit + df$se.fit)) +
     scale_color_brewer(palette='Set1', type='qual') +
     scale_fill_brewer(palette='Set1', type='qual') +
     theme_classic() +
     theme(axis.title.x=element_text(size=12, margin=margin(t=15)),
           axis.title.y=element_text(size=12, margin=margin(r=15)),
           axis.text=element_text(size=11))

path <- 'figs/fig2_prevalence.pdf'
pdf(path, width=8, height=5, onefile=FALSE)

reposition_legend(p0 + theme_classic(), 'top right', x=1-0.1, y=1-0.1)

dev.off()
