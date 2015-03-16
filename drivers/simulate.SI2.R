#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
#args <- c('/tmp/input.csv', '/tmp/tips.csv', '/tmp/output.nwk')

if (length(args) != 3) { stop('Usage: simulate.SI.R <input CSV> <tips CSV> <output NWK>') }
input.csv = args[1]  # simulation and model parameter settings
tips.csv = args[2]  # tip dates and states
output.nwk = args[3]

if (!file.exists(input.csv)) {
	stop('input file does not exist')
}
if (!file.exists(tips.csv)) {
	stop('tip label file does not exist')
}

n.cores <- 6  # for simulation in parallel

## default settings
n.tips = 100
nreps = 10
fgyResolution = 500.  # large value gives smaller time step
integrationMethod = 'adams'
t0 = 0
t_end = 30.*52  # weeks
t_break = 0.5  # relative time where transmission rate changes

N = 1000  # total population size

# model parameters
beta1 = 0.01
beta2 = 0.005
gamma = 1/520.
mu = 1/3640.



# parse settings from control file
inputs <- read.csv(input.csv, header=FALSE, na.strings='')
for (i in 1:nrow(inputs)) {
	eval(parse(text=paste(sep='', inputs[i,1], '<-', inputs[i,2])))
}

# initial population frequencies
S = N-1
I = 1
x0 <- c(I=I, S=S)
if (any(x0 < 0)) {
	stop('Population sizes cannot be less than 0.')
}

parms <- list(beta=beta1, gamma=gamma, mu=mu)
if (any(parms<0)) {
	stop ('No negative values permitted for model parameters.')
}


# define ODE system
demes <- c('I')

births <- rbind(c('parms$beta*S*I / (S+I)'))
rownames(births)=colnames(births) <- demes

migrations <- rbind(c('0'))
rownames(migrations)=colnames(migrations) <- demes

deaths <- c('(parms$mu+parms$gamma)*I')
names(deaths) <- demes

# dynamics for susceptible class (S)
nonDemeDynamics <- paste(sep='', '-parms$mu*S + parms$mu*S + (parms$mu+parms$gamma)*I', '-S*(parms$beta*I) / (S+I)')
names(nonDemeDynamics) <- 'S'

# read tip labels from file
tip.labels <- read.csv(tips.csv, header=FALSE, na.strings='')
names(tip.labels) <- c('tip.label', 'tip.height')
n.tips <- nrow(tip.labels)
if (sum(x0) < n.tips) {
	stop ('Population size is smaller than requested number of tips.')
}

# interpret missing tip dates as 0 (sampled at t.end)
tip.labels$tip.height[is.na(tip.labels$tip.height)] <- 0

if (max(tip.labels$tip.height) > t_end) {
	stop('Max tip height in labels file exceeds t.end setting.')
}
if (any(tip.labels$tip.height < 0)) {
	stop('Negative tip heights not allowed.')
}

# a vector indicating when each tip was sampled
#sampleTimes <- rep(t_end, times=n.tips)
sampleTimes <- t_end - tip.labels$tip.height

sampleStates <- matrix(1, nrow=n.tips, ncol=length(demes))
colnames(sampleStates) <- demes
rownames(sampleStates) <- 1:n.tips


# run this for numerical solution of ODE
m <- nrow(births)
maxSampleTime <- max(sampleTimes)


require(rcolgem, quietly=TRUE)

#tfgy <- make.fgy( t0, maxSampleTime, births, deaths, nonDemeDynamics,  x0,  migrations=migrations,  parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod )

# adjust fgyResolution for t_break
times <- seq(t0, t_end, length.out=fgyResolution)
fgyRes.1 <- round(fgyResolution * t_break)
fgyRes.2 <- fgyResolution - fgyRes.1


tfgy.1 <- make.fgy( t0, times[fgyRes.1], births, deaths, nonDemeDynamics,  x0,  migrations=migrations,  parms=parms, fgyResolution = fgyRes.1, integrationMethod = integrationMethod )

x1 <- tfgy.1[[5]][fgyRes.1, 2:3]
parms$beta <- beta2
tfgy.2 <- make.fgy( times[fgyRes.1+1], maxSampleTime, births, deaths, nonDemeDynamics,  x1,  migrations=migrations,  parms=parms, fgyResolution = fgyRes.2, integrationMethod = integrationMethod )

# are these the same?
#plot(tfgy[[5]][,1], tfgy[[5]][,2], type='l', ylim=c(0,1000))
#lines(tfgy[[5]][,1], tfgy[[5]][,3], col='red')
#plot(tfgy.1[[5]][,1], tfgy.1[[5]][,2], lty=2, ylim=c(0, 1000), xlim=c(0, maxSampleTime))
#points(tfgy.1[[5]][,1], tfgy.1[[5]][,3], lty=2, col='red')
#points(tfgy.2[[5]][,1], tfgy.2[[5]][,2], lty=2)
#points(tfgy.2[[5]][,1], tfgy.2[[5]][,3], lty=2, col='red')


#trees <- simulate.binary.dated.tree(births=births, deaths=deaths, nonDemeDynamics=nonDemeDynamics, t0=0, x0=x0, sampleTimes=sampleTimes, sampleStates=sampleStates, migrations=migrations, parms=parms, n.reps=10)

# reconstitute entire tfgy
y.times <- c(tfgy.2[[1]], tfgy.1[[1]])
y.births <- c(tfgy.2[[2]], tfgy.1[[2]])
y.migrations <- c(tfgy.2[[3]], tfgy.1[[3]])
y.demeSizes <- c(tfgy.2[[4]], tfgy.1[[4]])

# times, births, migrations, demeSizes
#trees <- simulate.binary.dated.tree.fgy( tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, integrationMethod = integrationMethod, n.reps=nreps)
trees <- simulate.binary.dated.tree.fgy(y.times, y.births, y.migrations, y.demeSizes, sampleTimes, sampleStates, integrationMethod, nreps)

'multiPhylo' -> class(trees)

write.tree(trees, file=output.nwk, append=FALSE)
