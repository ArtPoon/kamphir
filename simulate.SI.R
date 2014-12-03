library(rcolgem)

demes <- c('I1', 'I2')

births <- rbind(c('parms$beta1*S*I1 / (S+I1+I2)', '0'), c('parms$beta2*S*I2 / (S+I1+I2)', '0'))
rownames(births)=colnames(births) <- demes

migrations <- rbind(c('0', 'parms$alpha * I1'), c('0', '0'))
rownames(migrations)=colnames(migrations) <- demes

deaths <- c('parms$mu*I1', '(parms$mu+parms$gamma)*I2')
names(deaths) <- demes

# dynamics for susceptible class (S)
nonDemeDynamics <- paste(sep='', '-parms$mu*S + parms$mu*(S+I1+I2)', '-S*(parms$beta1*I1 + parms$beta2*I2) / (S+I1+I2)')
names(nonDemeDynamics) <- 'S'



# sample 100 individuals
n.tips <- 1000
sampleTimes <- rep(30.*52, times=n.tips)

# this is a binary-valued matrix where number of columns equals demes
# and number of rows equals tips in tree
# row names should correspond to tip labels
# for example:
# sampleStates <- cbind( rep(1, length(sampleTimes)), rep(0, length(sampleTimes)) )
sampleStates <- matrix(0, nrow=n.tips, ncol=length(demes))
colnames(sampleStates) <- demes
for (i in 1:n.tips) {
	sampleStates[i, i%%2+1] <- 1
}
rownames(sampleStates) <- paste(1:n.tips, rep(c(1,2), times=n.tips/2), sep='_')


# run this for numerical solution of ODE
m <- nrow(births)
maxSampleTime <- max(sampleTimes)
fgyResolution <- 100.  # large value gives smaller time step
integrationMethod <- 'rk4'
parms <- list(beta1=0.01, beta2=0.001, alpha=1/1352., gamma=1/520., mu=1/3640.)
t0 <- 0

# simulate population of 1000 individuals
x0 <- c(I1=1, I2=0, S=9999)

tfgy <- make.fgy( t0, maxSampleTime, births, deaths, nonDemeDynamics,  x0,  migrations=migrations,  parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod )

#trees <- simulate.binary.dated.tree(births=births, deaths=deaths, nonDemeDynamics=nonDemeDynamics, t0=0, x0=x0, sampleTimes=sampleTimes, sampleStates=sampleStates, migrations=migrations, parms=parms, n.reps=10)

trees <- simulate.binary.dated.tree.fgy( tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, integrationMethod = integrationMethod, n.reps=10)
'multiPhylo' -> class(trees)



