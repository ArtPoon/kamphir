#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
#args <- c('/tmp/input.csv', '/tmp/tips.csv', '/tmp/output.nwk')

if (length(args) != 3) { stop('Usage: mapCasesToFSA.R <input CSV> <tips CSV> <output NWK>') }
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
nreps = 10
fgyResolution = 500.  # large value gives smaller time step
integrationMethod = 'rk4'
t0 = 0
t_end = 30.*52  # weeks

# model parameters (per week)
beta1 = 0.01
beta2 = 0.001
alpha = 0.01  # transition rate from acute to chronic
gamma = 1/520. # excess mortality 
mu = 1/3640.  # baseline mortality

# parse settings from control file
inputs <- read.csv(input.csv, header=FALSE, na.strings='')
for (i in 1:nrow(inputs)) {
	eval(parse(text=paste(sep='', inputs[i,1], '<-', inputs[i,2])))
}

parms <- list(alpha=alpha, beta=beta, gamma=gamma, mu=mu)
if (any(parms<0)) {
	stop ('No negative values permitted for model parameters.')
}


# initial population frequencies
N = 1000  # total population size
S = N-1
I1 = 1
I2 = 0
x0 <- c(I1=I1, I2=I2, S=S)
if (any(x0 < 0)) {
	stop('Population sizes cannot be less than 0.')
}



# define the ODE system
demes <- c('I1', 'I2')

births <- rbind(c('parms$beta1*S*I1 / (S+I1+I2)', '0'), c('parms$beta2*S*I2 / (S+I1+I2)', '0'))
rownames(births)=colnames(births) <- demes

migrations <- rbind(c('0', 'parms$alpha * I1'), c('0', '0'))
rownames(migrations)=colnames(migrations) <- demes

deaths <- c('parms$mu*I1', '(parms$mu+parms$gamma)*I2')
names(deaths) <- demes

# dynamics for susceptible class (S) - complete replacement
nonDemeDynamics <- paste(sep='', '-parms$mu*S + parms$mu*(S+I1+I2)', '-S*(parms$beta1*I1 + parms$beta2*I2) / (S+I1+I2)')
names(nonDemeDynamics) <- 'S'


# parse tip labels
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
sampleTimes <- t_end - tip.labels$tip.height

# binary-valued matrix indicating state of each tip
sampleStates <- matrix(0, nrow=n.tips, ncol=length(demes))
colnames(sampleStates) <- demes
for (i in 1:n.tips) {
	sampleStates[i, tip.labels$tip.label[i]] <- 1
}
rownames(sampleStates) <- paste(1:n.tips, tip.labels$tip.label, sep='_')


# numerical solution of ODE
library(rcolgem, quietly=TRUE)

m <- nrow(births)
maxSampleTime <- max(sampleTimes)
tfgy <- make.fgy( t0, maxSampleTime, births, deaths, nonDemeDynamics,  x0,  migrations=migrations,  parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod )

# simulate trees
trees <- simulate.binary.dated.tree.fgy( tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, integrationMethod = integrationMethod, n.reps=n.reps, n.cores=ncores)

'multiPhylo' -> class(trees)

write.tree(trees, file=output.nwk, append=FALSE)


