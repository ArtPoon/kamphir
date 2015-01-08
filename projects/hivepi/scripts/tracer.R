# default skip lines starting with '#' (comment char)
#log <- read.table('~/git/kamphir/projects/hivepi/beast/MASTER1.log', sep='\t', header=TRUE)
log <- read.table('~/git/kamphir/projects/hivepi/beast/SIRTree.n100.log', sep='\t', header=TRUE)

# true R0 is 3.33 (0.001*999/0.3)
log$beta <- log$R0Es * log$becomeUninfectiousRateEs / log$S0Es

par(family='sans', mfrow=c(2,2))
#plot(log$Sample, log$likelihood, type='b', cex=1.2, ylim=c(log$likelihood[20], max(log$likelihood)))
plot(log$Sample, log$beta, type='l'); abline(h=0.001, lty=2, col='red')
plot(log$Sample, log$S0, type='l', log='y'); abline(h=999, lty=2, col='red')
plot(log$Sample, log$become, type='l'); abline(h=0.3, lty=2, col='red')
plot(log$Sample, log$samplingProportionEs, type='l'); abline(h=0.15, lty=2, col='red')

require(coda)
mc <- mcmc(data=log[100:nrow(log), which(is.element(names(log), c('Sample', 'likelihood', 'clockRate', 'kappa', 'beta', 'S0Es', 'becomeUninfectiousRateEs', 'samplingProportionEs')))])

# using MASTER to directly fit simulated tree (274 tips)
kam <- read.table('~/git/kamphir/projects/hivepi/data/MASTER.SIR.beta1E-3.gamma03.phi015.N1000.n100.log.4', header=TRUE, sep='\t')
kam$R0 <- kam$beta * (kam$N-1) / kam$gamma
# true value of R0 is 3.33

# using rcolgem to directly fit simulated tree (1000 tips)
exact <- read.table('~/git/kamphir/projects/hivepi/data/SI.beta001.N10E5.log.11', header=TRUE, sep='\t')