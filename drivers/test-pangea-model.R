
# evaluation period is for 4 years (2015-2018) 
# this period ends 1 year before end of simulation (2019)
# there are some simulations that end Dec 2017, so evaluation period
# ends Dec 2016 (2 years)
# HIV first introduced into population in 1980

demes <- c('I0', 'J0', 'K0', 'I1', 'J1', 'K1', 'I2', 'J2', 'K2')

# I = acute, J = asymptomatic, K = chronic; S = susceptible
# 0 = source population, 1 = main regional population, 2 = high risk minority population


# mixing rates - NOTE rho is now the proportion reserved for contacts OUTSIDE of group
p11 <- '((1-parms$rho1) * parms$c1*(S1+I1+J1+K1)*(1-parms$rho1) / (parms$c1*(S1+I1+J1+K1)*(1-parms$rho1) + parms$c2*(S2+I2+J2+K2)*(1-parms$rho2)))'
p12 <- '(parms$rho1 + (1-parms$rho1) * parms$c2*(S2+I2+J2+K2)*(1-parms$rho2) / (parms$c1*(S1+I1+J1+K1)*(1-parms$rho1) + parms$c2*(S2+I2+J2+K2)*(1-parms$rho2)))'
p21 <- '(parms$rho2 + (1-parms$rho2) * parms$c1*(S1+I1+J1+K1)*(1-parms$rho1) / (parms$c1*(S1+I1+J1+K1)*(1-parms$rho1) + parms$c2*(S2+I2+J2+K2)*(1-parms$rho2)))'
p22 <- '((1-parms$rho2) * parms$c2*(S2+I2+J2+K2)*(1-parms$rho2) / (parms$c1*(S1+I1+J1+K1)*(1-parms$rho1) + parms$c2*(S2+I2+J2+K2)*(1-parms$rho2)))'

# birth rates
births <- rbind(
	c('parms$z*parms$beta1*parms$c1*I0/(S0+I0+J0+K0)*S0', '0', '0',  
	  '0', '0', '0',  
	  '0', '0', '0'),  # I0 ->
	c('parms$z*parms$beta2*parms$c1*J0/(S0+I0+J0+K0)*S0', '0', '0',  
	  '0', '0', '0',  
	  '0', '0', '0'),  # J0 ->
	c('parms$z*parms$beta3*parms$c1*K0/(S0+I0+J0+K0)*S0', '0', '0',  
	  '0', '0', '0',  
	  '0', '0', '0'),  # K0 ->
	c('0', '0', '0',
	  paste(sep='*', 'parms$z*parms$beta1*parms$c1', p11, 'I1/(S1+I1+J1+K1)*S1'), '0', '0', # I1 -> I1
	  paste(sep='*', 'parms$z*parms$beta1*parms$c2', p21, 'I1/(S1+I1+J1+K1)*S2'), '0', '0'  # I1 -> I2
	  ),
	c('0', '0', '0',
	  paste(sep='*', 'parms$z*parms$beta2*parms$c1', p11, 'J1/(S1+I1+J1+K1)*S1'), '0', '0', # J1 -> I1
	  paste(sep='*', 'parms$z*parms$beta2*parms$c2', p21, 'J1/(S1+I1+J1+K1)*S2'), '0', '0'  # J1 -> I2
	  ),
	c('0', '0', '0',
	  paste(sep='*', 'parms$z*parms$beta3*parms$c1', p11, 'K1/(S1+I1+J1+K1)*S1'), '0', '0', # K1 -> I1
	  paste(sep='*', 'parms$z*parms$beta3*parms$c2', p21, 'K1/(S1+I1+J1+K1)*S2'), '0', '0' # K1 -> I2
	  ),
	c('0', '0', '0',
	  paste(sep='*', 'parms$z*parms$beta1*parms$c1', p12, 'I2/(S2+I2+J2+K2)*S1'), '0', '0', # I2 -> I1
	  paste(sep='*', 'parms$z*parms$beta1*parms$c2', p22, 'I2/(S2+I2+J2+K2)*S2'), '0', '0' # I2 -> I2
	  ),
	c('0', '0', '0', 
	  paste(sep='*', 'parms$z*parms$beta2*parms$c1', p12, 'J2/(S2+I2+J2+K2)*S1'), '0', '0', # J2 -> J1
	  paste(sep='*', 'parms$z*parms$beta2*parms$c2', p22, 'J2/(S2+I2+J2+K2)*S2'), '0', '0' # J2 -> J2
	  ),
	c('0', '0', '0',
	  paste(sep='*', 'parms$z*parms$beta3*parms$c1', p12, 'K2/(S2+I2+J2+K2)*S1'), '0', '0', # K2 -> K1
	  paste(sep='*', 'parms$z*parms$beta3*parms$c2', p22, 'K2/(S2+I2+J2+K2)*S2'), '0', '0' # K2 -> K2
	  )
)

rownames(births)=colnames(births) <- demes

# migration rates (transition from acute to chronic, not spatial migration)
# note constant term is used to allow for early migration from source to sink while
#  infected population in sink is zero
migrations <- rbind(
           c('0', 'parms$alpha1*I0', '0', '(parms$mig*I1 + 1e-4) * I0/(S0+I0+J0+K0)', '0', '0', '0', '0', '0'), # I0->J0
           c('0', '0', 'parms$alpha2*J0', '0', 'parms$mig*J1*J0/(S0+I0+J0+K0)', '0', '0', '0', '0'), # J0->K0
           c('0', '0', '0', '0', '0', '0', '0', '0', '0'),  
           c('0', '0', '0', '0', 'parms$alpha1*I1', '0', '0', '0', '0'), # I1->J1
           c('0', '0', '0', '0', '0', 'parms$alpha2*J1', '0', '0', '0'), # J1->K1
           c('0', '0', '0', '0', '0', '0', '0', '0', '0'),
           c('0', '0', '0', '0', '0', '0', '0', 'parms$alpha1*I2', '0'), # I2->J2
           c('0', '0', '0', '0', '0', '0', '0', '0', 'parms$alpha2*J2'), # J2->K2
           c('0', '0', '0', '0', '0', '0', '0', '0', '0')
)
rownames(migrations)=colnames(migrations) <- demes

# death rates
deaths <- c(
	'parms$mu*I0', 'parms$mu*J0', '(parms$mu+parms$gamma)*K0', 
	'parms$mu*I1', 'parms$mu*J1', '(parms$mu+parms$gamma)*K1', 
	'parms$mu*I2', 'parms$mu*J2', '(parms$mu+parms$gamma)*K2'
)
names(deaths) <- demes



# dynamics of susceptible class
nonDemeDynamics <- c(
	'-parms$mu*S0 + parms$mu*(S0+I0+J0+K0) + parms$gamma*K0 - S0*parms$c1*parms$z*(parms$beta1*I0+parms$beta2*J0+parms$beta3*K0)/(S0+I0+J0+K0) + parms$mig*(I1*I0+J1*J0)/(S0+I0+J0+K0)', 
	paste(sep='', '-parms$mu*S1 + parms$mu*(S1+I1+J1+K1) + parms$gamma*K1 - S1*parms$c1*parms$z*(', p11, 
	  '*(parms$beta1*I1 + parms$beta2*J1 + parms$beta3*K1)/(S1+I1+J1+K1) + ', p12, 
	  '*(parms$beta1*I2 + parms$beta2*J2 + parms$beta3*K2)/(S2+I2+J2+K2))'),
	paste(sep='', '-parms$mu*S2 + parms$mu*(S2+I2+J2+K2) + parms$gamma*K2 - S2*parms$c2*parms$z*(', p21, 
	  '*(parms$beta1*I1 + parms$beta2*J1 + parms$beta3*K1)/(S1+I1+J1+K1) + ', p22, 
	  '*(parms$beta1*I2 + parms$beta2*J2 + parms$beta3*K2)/(S2+I2+J2+K2))')
)

names(nonDemeDynamics) <- c('S0', 'S1', 'S2')


# initialize population sizes
N0=1E6
S0=N0-1; I0=1; J0=0; K0=0
N1=1E4; p=0.05
S1=N1*(1-p); I1=0; J1=0; K1=0
S2=N1*p; I2=0; J2=0; K2=0
x0 = c(I0=I0, J0=J0, K0=K0, 
	I1=I1, J1=J1, K1=K1, 
	I2=I2, J2=J2, K2=K2, S0=S0, S1=S1, S2=S2)

## default settings
nreps = 1
fgyResolution = 500.  # large value gives smaller time step
integrationMethod = 'rk4'
t.end = 40.*52  # weeks
n.tips = 200

sampleTimes <- rep(t.end, times=n.tips)
maxSampleTime <- max(sampleTimes)

# numerical solution of ODE system
require(rcolgem)

eval.period <- 5*52.  # weeks
t0 <- 0
t1 <- (maxSampleTime-eval.period) / 2
fgyRes.1 <- round(fgyResolution * (t1-t0) / maxSampleTime)

parms <- list(alpha1=0.01, alpha2=0.001, gamma=1/520., mu=1/3640., beta1=0.03, beta2=0.00084, beta3=0.02, rho1=0.1, rho2=0.9, c1=1.0, c2=5.0, mig=0.05/52, z=1)

tfgy.1 <- make.fgy(t0, t1, births, deaths, nonDemeDynamics, x0, migrations=migrations, parms=parms, fgyResolution=fgyRes.1, integrationMethod=integrationMethod)

# use state of system at last time point to initialize next interval
x1 <- tfgy.1[[5]][fgyRes.1, 2:ncol(tfgy.1[[5]])]
t2 <- maxSampleTime-eval.period
fgyRes.2 <- round(fgyResolution * (t2-t1) / maxSampleTime)
parms$z <- 1.3

tfgy.2 <- make.fgy(t1, t2, births, deaths, nonDemeDynamics, x1, migrations=migrations, parms=parms, fgyResolution=fgyRes.2, integrationMethod=integrationMethod)

x2 <- tfgy.2[[5]][fgyRes.2, 2:ncol(tfgy.2[[5]])]
t3 <- maxSampleTime
fgyRes.3 <- fgyResolution - fgyRes.1 - fgyRes.2
parms$z <- 0.5

tfgy.3 <- make.fgy(t2, t3, births, deaths, nonDemeDynamics, x2, migrations=migrations, parms=parms, fgyResolution=fgyRes.3, integrationMethod=integrationMethod)


# sanity checks
par(mfrow=c(3,1), mar=c(5,5,0,2))

mat <- rbind(c(11,2,3,4), c(12,5,6,7), c(13,8,9,10))
for (i in c(1:nrow(mat))) {
	x <- mat[i,]
	df <- tfgy.1[[5]]
	t <- df[,1] / 52.
	plot(t, df[,x[1]], ylim=c(0, max(df[,x[1]])), type='l', col='forestgreen', lwd=2, xlim=c(0, maxSampleTime)/52.)
	lines(t, df[,x[2]], col='red', lwd=2)
	lines(t, df[,x[3]], col='orange', lwd=2)
	lines(t, df[,x[4]], col='black', lwd=2)
	
	df <- tfgy.2[[5]]
	t <- df[,1] / 52.
	lines(t, df[,x[1]], col='forestgreen', lwd=2, lty=2)
	lines(t, df[,x[2]], col='red', lwd=2, lty=2)
	lines(t, df[,x[3]], col='orange', lwd=2, lty=2)
	lines(t, df[,x[4]], col='black', lwd=2, lty=2)
	
	df <- tfgy.3[[5]]
	t <- df[,1] / 52.
	lines(t, df[,x[1]], col='forestgreen', lwd=1, lty=3)
	lines(t, df[,x[2]], col='red', lwd=1, lty=3)
	lines(t, df[,x[3]], col='orange', lwd=1, lty=3)
	lines(t, df[,x[4]], col='black', lwd=1, lty=3)
}


## check for stability of source population
#all(round(apply(tfgy[[5]][,c(2:4, 11)], 1, sum))==1E6)

## growth of sink population
#apply(tfgy[[5]][,c(5:7,12)], 1, sum)

# which deme was each tip sampled from
sampleStates <- matrix(0, nrow=n.tips, ncol=length(demes))
colnames(sampleStates) <- demes
rownames(sampleStates) <- 1:n.tips

demes.t.end <- tfgy.3[[4]][[1]]
sampled.demes <- which(!grepl("0$", names(demes.t.end)))
demes.sample <- sample(rep(sampled.demes, times=round(demes.t.end[sampled.demes])), size=n.tips)
#demes.sample <- sample(rep(1:length(demes), times=round(demes.t.end)), size=n.tips)
sampleStates <- matrix(0, nrow=n.tips, ncol=length(demes))
colnames(sampleStates) <- demes
for (i in 1:n.tips) {
	sampleStates[i, demes.sample[i]] <- 1
}
rownames(sampleStates) <- paste(1:n.tips, demes.sample, sep='_')

# reconstritute the entire ODE
#y.times <- c(tfgy.3[[1]], tfgy.2[[1]], tfgy.1[[1]])
# concatenation fails to get equal increments
y.times <- rev(seq(0, maxSampleTime, length.out=fgyResolution))
y.births <- c(tfgy.3[[2]], tfgy.2[[2]], tfgy.1[[2]])
y.migrations <- c(tfgy.3[[3]], tfgy.2[[3]], tfgy.1[[3]])
y.demeSizes <- c(tfgy.3[[4]], tfgy.2[[4]], tfgy.1[[4]])

trees <- simulate.binary.dated.tree.fgy( y.times, y.births, y.migrations, y.demeSizes, sampleTimes, sampleStates, integrationMethod = integrationMethod, n.reps=nreps)
plot(trees[[1]])
plot(trees[[1]]$height/52.)
'multiPhylo' -> class(trees)
write.tree(trees, file=output.nwk, append=FALSE)

# probabilities of ancestral demes - I'm not sure this works...
#ll.tree <- coalescent.log.likelihood.fgy(trees[[1]], rev(y.times), rev(y.births), rev(y.migrations), rev(y.demeSizes), returnTree=TRUE)
#plot(apply(ll.tree[[2]]$lstate[,1:3], 1, sum))
#plot(apply(ll.tree[[2]]$lstate[,7:9], 1, sum))
