# Plot sheet sample designs

path <- 'figs/fig4_sheetsamp.pdf'
pdf(path, width=6, height=6)

par(mfrow=c(2,2), mar=c(2,2,2,2), oma = c(0, 0, 0, 0))
sheetsamp(r=30, s=1.9, nsheets=10, type='quadrant', plot=T, main="Quadrant")
sheetsamp(r=30, s=1, nsheets=100, type='uniform', plot=T, main="Uniform")
sheetsamp(r=30, s=1, nsheets=100, type='stratified', plot=T, main="Stratified")
sheetsamp(r=30, s=1, nsheets=100, type='random', plot=T, main="Random")

dev.off()