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

N = 3000  # total population size
n.tips <- 300
p = 0.5  # frequency of risk group 1

# model parameters
beta = 0.01
gamma = 1/520.
mu = 1/3640.
c1 = 2.0
c2 = 1.0
rho = 0.9


# parse settings from control file
inputs <- read.csv(input.csv, header=FALSE, na.strings='')
for (i in 1:nrow(inputs)) {
	eval(parse(text=paste(sep='', inputs[i,1], '<-', inputs[i,2])))
}

# initial population frequencies
S1 = p*N-1
S2 = (1-p)*N
I1 = 1
I2 = 0
x0 <- c(I1=I1, I2=I2, S1=S1, S2=S2)
if (any(x0 < 0)) {
	stop('Population sizes cannot be less than 0.')
}

parms <- list(beta=beta, gamma=gamma, mu=mu, c1=c1, c2=c2, rho=rho)
if (any(parms<0)) {
	stop ('No negative values permitted for model parameters.')
}


# define ODE system

demes <- c('I1', 'I2')  # two risk groups

p11 <- '(parms$rho + (1-parms$rho) * parms$c1*(S1+I1) / (parms$c1*(S1+I1) + parms$c2*(S2+I2)))'
p12 <- '(1-parms$rho) * parms$c2*(S2+I2) / (parms$c1*(S1+I1) + parms$c2*(S2+I2))'
p21 <- '(1-parms$rho) * parms$c1*(S1+I1) / (parms$c1*(S1+I1) + parms$c2*(S2+I2))'
p22 <- '(parms$rho + (1-parms$rho) * parms$c2*(S2+I2) / (parms$c1*(S1+I1) + parms$c2*(S2+I2)))'

births <- rbind(c(paste(sep='*', 'parms$beta*parms$c1', p11, 'I1/(S1+I1)*S1'), 
				paste(sep='*', 'parms$beta*parms$c2', p21, 'I1/(S1+I1)*S2')),
				c(paste(sep='*', 'parms$beta*parms$c1', p12, 'I2/(S2+I2)*S1'),
				paste(sep='*', 'parms$beta*parms$c2', p22, 'I2/(S2+I2)*S2')))

rownames(births)=colnames(births) <- demes

migrations <- rbind(c('0', '0'), c('0', '0'))
rownames(migrations)=colnames(migrations) <- demes

deaths <- c('(parms$mu+parms$gamma)*I1', '(parms$mu+parms$gamma)*I2')
names(deaths) <- demes

# dynamics for susceptible classes (S)
nonDemeDynamics <- c(paste(sep='', '-parms$mu*S1 + parms$mu*S1 + (parms$mu+parms$gamma)*I1', paste(sep='*', '-S1*(parms$beta*parms$c1', p11, 'I1/(S1+I1) + parms$beta*parms$c1', p12, 'I2/(S2+I2))')),
	paste(sep='', '-parms$mu*S2 + parms$mu*S2 + (parms$mu+parms$gamma)*I2', paste(sep='*', '-S2*(parms$beta*parms$c2', p21, 'I1/(S1+I1) + parms$beta*parms$c2', p22, 'I2/(S2+I2))')))

names(nonDemeDynamics) <- c('S1', 'S2')


## parse tip labels from file
# tip label should be an integer 1..n where n is number of demes
# tip date should be some numerical value < t.end
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


# numerical solution of ODE
require(rcolgem, quietly=TRUE)

m <- nrow(births)
maxSampleTime <- max(sampleTimes)

tfgy <- make.fgy( t0, maxSampleTime, births, deaths, nonDemeDynamics,  x0,  migrations=migrations,  parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod)
# returns list of [1] times, [2] births, [3] migrations, [4] demeSizes, 
# and [5] ODE solution
# NOTE items 1-4 are in reverse time


# use prevalence of respective infected classes to determine sample states
demes.t.end <- tfgy[[4]][[1]]
if (sum(demes.t.end) < n.tips) {
	stop('Number of infected individuals at t.end less than n.tips')
}

demes.sample <- sample(rep(1:length(demes), times=round(demes.t.end)), size=n.tips)
sampleStates <- matrix(0, nrow=n.tips, ncol=length(demes))
colnames(sampleStates) <- demes
for (i in 1:n.tips) {
	sampleStates[i, demes.sample[i]] <- 1
}
rownames(sampleStates) <- paste(1:n.tips, demes.sample, sep='_')


# simulate trees
trees <- simulate.binary.dated.tree.fgy( tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, integrationMethod = integrationMethod, n.reps=nreps, n.cores=n.cores)
'multiPhylo' -> class(trees)

write.tree(trees, file=output.nwk, append=FALSE)

