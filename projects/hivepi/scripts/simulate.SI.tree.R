#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

if (length(args) != 2) { stop('Usage: Rscript simulate.SI.tree.R <#tips> <outfile>') }
n.tips = as.integer(args[1])
output.nwk = args[2]

# simulate tree under SI model
n.reps = 1
n.cores = 1
fgyResolution = 500.  # large value gives smaller time step
integrationMethod = 'rk4'

t0 = 0
t.end = 30.*52  # weeks

N = 10000  # total population size
S = N-1
I = 1
x0 <- c(I=I, S=S)  # initial population frequencies

# model parameters
beta = 0.01
gamma = 1/520.
mu = 1/3640.
parms <- list(beta=beta, gamma=gamma, mu=mu)

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

# assume collection dates uniformly distributed over last 5 years
sampleTimes <- rep(t.end, times=n.tips)
#sampleTimes <- rep(t.end, times=n.tips) - runif(n.tips, min=0, max=5.*52)
sampleStates <- matrix(1, nrow=n.tips, ncol=length(demes))
colnames(sampleStates) <- demes
rownames(sampleStates) <- paste(sep='_', 1:n.tips, as.integer(sampleTimes))



m <- nrow(births)
maxSampleTime <- max(sampleTimes)

require(rcolgem, quietly=TRUE)
# run this for numerical solution of ODE
tfgy <- make.fgy( t0, maxSampleTime, births, deaths, nonDemeDynamics,  x0,  migrations=migrations,  parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod )

trees <- simulate.binary.dated.tree.fgy( tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, integrationMethod = integrationMethod, n.reps=n.reps)

tree <- trees[[1]]
tree$tip.label <- rownames(sampleStates)[as.integer(tree$tip.label)]

write.tree(tree, file=output.nwk)
