import sys
from rpy2.rinterface import set_readconsole
set_readconsole(None)

import rpy2.robjects as robjects  # R is instantiated upon load module

class Rcolgem ():
    def __init__ (self, ncores, nreps, t0=0, fgy_resolution=500., integration_method='rk4', seed=None):
        # load Rcolgem package
        robjects.r("require(rcolgem, quietly=TRUE)")

        # default settings
        robjects.r('n.cores=%d; nreps=%d; fgyResolution=%d; integrationMethod="%s"; t0=%f' % (
            ncores, nreps, fgy_resolution, integration_method, t0))

        # set up parallelization environment
        robjects.r("require(parallel, quietly=TRUE)")

        if (ncores > 1):
            robjects.r("cl <- makeCluster(%d, 'FORK')" % (ncores,))
            if seed is not None:
                robjects.r("clusterSetRNGStream(cl, {})".format(seed))
        else:
            robjects.r("cl <- NULL")
            if seed is not None:
                robjects.r("set.seed({})".format(seed))

    def init_SI_model (self):
        """
        Defines a susceptible-infected-recovered model in rcolgem.
        :return:
        """
        # define ODE system - as strings, these will be evaluated with new parameters
        robjects.r('demes <- c("I")')

        robjects.r('births <- rbind(c("parms$beta*S*I / (S+I)"))')
        robjects.r('rownames(births) <- colnames(births) <- demes')

        robjects.r('migrations <- rbind(c("0"))')
        robjects.r('rownames(migrations)=colnames(migrations) <- demes')

        robjects.r("deaths <- c('(parms$mu+parms$gamma)*I')")
        robjects.r('names(deaths) <- demes')

        robjects.r("nonDemeDynamics <- paste(sep='', '-parms$mu*S + parms$lambd*S + "
                   "(parms$mu+parms$gamma)*I', '-S*(parms$beta*I) / (S+I)')")
        robjects.r("names(nonDemeDynamics) <- 'S'")


    def simulate_SI_trees (self, params, tree_height, tip_heights, post=False):
        """
        Simulate coalescent trees under the SI model.
        :param tip_heights:
        :return: List of trees; if post=True, then a tuple of ([trees], tfgy)
        """
        # set parameters
        robjects.r('N=%f; beta=%f; gamma=%f' % (params['N'], params['beta'], params['gamma']))
        robjects.r('mu=%f; lambd=%f' % (params['mu'], params.get('lambd', params['mu'])))
        robjects.r('S = N-1')
        robjects.r('I = 1')
        robjects.r('x0 <- c(I=I, S=S)')
        robjects.r('parms <- list(beta=beta, gamma=gamma, mu=mu, lambd=lambd)')

        robjects.r("n.tips <- %d" % len(tip_heights))
        robjects.r("tip.heights <- c(%s)" % ','.join(map(str, tip_heights)))
        robjects.r("t_end <- %f" % (tree_height,))

        robjects.r("sampleTimes <- t_end - tip.heights")
        robjects.r("sampleStates <- matrix(1, nrow=n.tips, ncol=length(demes))")
        robjects.r("colnames(sampleStates) <- demes")
        robjects.r("rownames(sampleStates) <- 1:n.tips")

        robjects.r("m <- nrow(births)")
        robjects.r("maxSampleTime <- max(sampleTimes)")

        # solve ODE
        robjects.r("tfgy <- make.fgy( t0, maxSampleTime, births, deaths, nonDemeDynamics, x0, migrations=migrations, "
                   "parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod)")

        n_inf = robjects.r("tfgy[[4]][[1]]")[0]
        if n_inf < len(tip_heights):
            # number of infected at end of simulation is less than number of tips
            sys.stderr.write("Expected {} infected at end of simulation, but got {}".format(len(tip_heights), n_inf))
            return []

        # simulate trees
        try:
            robjects.r("trees <- simulate.binary.dated.tree.fgy( tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], "
                       "sampleTimes, sampleStates, integrationMethod = integrationMethod, "
                       "n.reps=nreps, cluster=cl)")
        except:
            return []

        robjects.r("'multiPhylo' -> class(trees)")
        try:
            retval = robjects.r("lapply(trees, write.tree)")
        except:
            return []

        trees = map(lambda x: str(x).split()[-1].strip('" '), retval)
        if post:
            return (trees, robjects.r("tfgy[[5]]"))
        else:
            return trees

    def simulate_SI2_trees(self, params, tree_height, tip_heights, post=False):
        """
        Simulate coalescent trees under a two-phase SI model.
        :param params:
        :param tip_heights:
        :return:
        """

        # set parameters
        robjects.r('N=%f; beta1=%f; beta2=%f' % (params['N'], params['beta1'], params['beta2']))
        robjects.r('gamma=%f; mu=%f; lambd=%f' % (params['gamma'], params['mu'], params.get('lambd', params['mu'])))
        robjects.r('t_end=%f; t_break=%f' % (tree_height, params['t_break']))

        # adjust fgyResolution for t_break
        robjects.r("times <- seq(t0, t_end, length.out=fgyResolution)")
        robjects.r("fgyRes.1 <- round(fgyResolution * t_break)")
        robjects.r("fgyRes.2 <- fgyResolution - fgyRes.1")

        # if break is too close to either limit, return single ODE solution
        tp1, tp2 = robjects.r("c(fgyRes.1, fgyRes.2)")
        if tp1 < 3:
            params2 = dict((k, v) for k, v in params.iteritems())  # deep copy
            params2.update({'beta': params['beta2']})
            return self.simulate_SI_trees(params2, tree_height, tip_heights, post)
        if tp2 < 3:
            params2 = dict((k, v) for k, v in params.iteritems())  # deep copy
            params2.update({'beta': params['beta1']})
            return self.simulate_SI_trees(params2, tree_height, tip_heights, post)

        # set model parameters
        robjects.r('S = N-1')
        robjects.r('I = 1')
        robjects.r('x0 <- c(I=I, S=S)')
        robjects.r('parms <- list(beta=beta1, gamma=gamma, mu=mu, lambd=lambd)')

        robjects.r("n.tips <- %d" % len(tip_heights))
        robjects.r("tip.heights <- c(%s)" % ','.join(map(str, tip_heights)))

        robjects.r("sampleTimes <- t_end - tip.heights")
        robjects.r("sampleStates <- matrix(1, nrow=n.tips, ncol=length(demes))")
        robjects.r("colnames(sampleStates) <- demes")
        robjects.r("rownames(sampleStates) <- 1:n.tips")

        robjects.r("m <- nrow(births)")
        robjects.r("maxSampleTime <- max(sampleTimes)")

        # solve first ODE
        robjects.r("tfgy.1 <- make.fgy( t0, times[fgyRes.1], births, deaths, nonDemeDynamics, x0, "
                   "migrations=migrations, parms=parms, fgyResolution = fgyRes.1, "
                   "integrationMethod = integrationMethod )")

        # update model parameter with second beta
        robjects.r("x1 <- tfgy.1[[5]][fgyRes.1, 2:3]")
        robjects.r("parms$beta <- beta2")

        # solve second ODE
        robjects.r("tfgy.2 <- make.fgy( times[fgyRes.1+1], maxSampleTime, births, deaths, nonDemeDynamics, x1, "
                   "migrations=migrations, parms=parms, fgyResolution = fgyRes.2, "
                   "integrationMethod = integrationMethod)")

        n_inf = robjects.r("tfgy.2[[4]][[1]]")[0]
        if n_inf < len(tip_heights):
            # number of infected at end of simulation is less than number of tips
            return []

        # reconstitute the entire tfgy
        robjects.r("y.times <- c(tfgy.2[[1]], tfgy.1[[1]])")
        robjects.r("y.births <- c(tfgy.2[[2]], tfgy.1[[2]])")
        robjects.r("y.migrations <- c(tfgy.2[[3]], tfgy.1[[3]])")
        robjects.r("y.demeSizes <- c(tfgy.2[[4]], tfgy.1[[4]])")

        # simulate trees
        try:
            robjects.r("trees <- simulate.binary.dated.tree.fgy(y.times, y.births, y.migrations, y.demeSizes, "
                       "sampleTimes, sampleStates, integrationMethod=integrationMethod, "
                       "n.reps=nreps, cluster=cl)")
        except:
            return []

        # convert R objects into Python strings in Newick format
        robjects.r("class(trees) <- 'multiPhylo'")
        try:
            retval = robjects.r("lapply(trees, write.tree)")
        except:
            # error converting trees
            return []

        trees = map(lambda x: str(x).split()[-1].strip('" '), retval)
        if post:
            return (trees, robjects.r("rbind(tfgy.1[[5]], tfgy.2[[5]])"))
        else:
            return trees


    def init_DiffRisk_model(self):
        """
        Define ODE system for differential risk SI model.
        """

        robjects.r("demes <- c('I1', 'I2')")

        robjects.r("p11 <- '(parms$rho + (1-parms$rho) * parms$c1*(S1+I1) / (parms$c1*(S1+I1) + parms$c2*(S2+I2)))'")
        robjects.r("p12 <- '(1-parms$rho) * parms$c2*(S2+I2) / (parms$c1*(S1+I1) + parms$c2*(S2+I2))'")
        robjects.r("p21 <- '(1-parms$rho) * parms$c1*(S1+I1) / (parms$c1*(S1+I1) + parms$c2*(S2+I2))'")
        robjects.r("p22 <- '(parms$rho + (1-parms$rho) * parms$c2*(S2+I2) / (parms$c1*(S1+I1) + parms$c2*(S2+I2)))'")
        robjects.r("births <- rbind(c(paste(sep='*', 'parms$beta*parms$c1', p11, 'I1/(S1+I1)*S1'),"
				   "paste(sep='*', 'parms$beta*parms$c2', p21, 'I1/(S1+I1)*S2')),"
		    	   "c(paste(sep='*', 'parms$beta*parms$c1', p12, 'I2/(S2+I2)*S1'),"
			       "paste(sep='*', 'parms$beta*parms$c2', p22, 'I2/(S2+I2)*S2')))")
        robjects.r("rownames(births)=colnames(births) <- demes")

        robjects.r("migrations <- rbind(c('0', '0'), c('0', '0'))")
        robjects.r("rownames(migrations)=colnames(migrations) <- demes")

        robjects.r("deaths <- c('(parms$mu+parms$gamma)*I1', '(parms$mu+parms$gamma)*I2')")
        robjects.r("names(deaths) <- demes")

        robjects.r("nonDemeDynamics <- c(paste(sep='', '-parms$mu*S1 + parms$mu*S1 + (parms$mu+parms$gamma)*I1', "
                   "paste(sep='*', '-S1*(parms$beta*parms$c1', p11, 'I1/(S1+I1) + parms$beta*parms$c1', p12, "
                   "'I2/(S2+I2))')), paste(sep='', '-parms$mu*S2 + parms$mu*S2 + (parms$mu+parms$gamma)*I2', "
                   "paste(sep='*', '-S2*(parms$beta*parms$c2', p21, 'I1/(S1+I1) + parms$beta*parms$c2', p22, "
                   "'I2/(S2+I2))')))")
        robjects.r("names(nonDemeDynamics) <- c('S1', 'S2')")

    def simulate_DiffRisk_trees(self, params, tree_height, tip_heights, post=False):
        """

        :param params:
        :param tip_heights:
        :return:
        """

        # set parameters
        robjects.r('N=%f; beta=%f; c1=%f; c2=%f' % (params['N'], params['beta'], params['c1'], params['c2']))
        robjects.r('rho=%f; p=%f; gamma=%f; mu=%f' % (params['rho'], params['p'], params['gamma'], params['mu']))
        robjects.r('t_end=%f' % (tree_height,))

        # update model parameters
        robjects.r("S1=p*N-1; S2=(1-p)*N; I1=1; I2=0")
        robjects.r("x0 <- c(I1=I1, I2=I2, S1=S1, S2=S2)")
        robjects.r("parms <- list(beta=beta, gamma=gamma, mu=mu, c1=c1, c2=c2, rho=rho)")

        robjects.r("n.tips <- %d" % len(tip_heights))
        robjects.r("tip.heights <- c(%s)" % ','.join(map(str, tip_heights)))

        robjects.r("sampleTimes <- t_end - tip.heights")
        robjects.r("sampleStates <- matrix(1, nrow=n.tips, ncol=length(demes))")
        robjects.r("colnames(sampleStates) <- demes")
        robjects.r("rownames(sampleStates) <- 1:n.tips")

        robjects.r("m <- nrow(births)")
        robjects.r("maxSampleTime <- max(sampleTimes)")

        # solve ODE
        robjects.r("tfgy <- make.fgy( t0, maxSampleTime, births, deaths, nonDemeDynamics, x0, migrations=migrations, "
                   "parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod)")

        # use prevalence of respective infected classes at end of simulation to determine sample states
        robjects.r("demes.t.end <- tfgy[[4]][[1]]")

        if robjects.r("sum(demes.t.end)")[0] < len(tip_heights):
            # number of infected individuals at end of simulation is less than number of tips
            return []

        robjects.r("demes.sample <- sample(rep(1:length(demes), times=round(demes.t.end)), size=n.tips)")
        robjects.r("sampleStates <- matrix(0, nrow=n.tips, ncol=length(demes))")
        robjects.r("colnames(sampleStates) <- demes")
        robjects.r("for (i in 1:n.tips) { sampleStates[i, demes.sample[i]] <- 1 }")
        robjects.r("rownames(sampleStates) <- paste(1:n.tips, demes.sample, sep='_')")

        # simulate trees
        try:

            robjects.r("trees <- simulate.binary.dated.tree.fgy( tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], "
                       "sampleTimes, sampleStates, integrationMethod = integrationMethod, "
                       "n.reps=nreps, cluster=cl)")
        except:
            return []

        # convert R objects into Python strings in Newick format
        robjects.r("class(trees) <- 'multiPhylo'")
        try:
            retval = robjects.r("lapply(trees, write.tree)")
        except:
            # error converting trees
            return []

        trees = map(lambda x: str(x).split()[-1].strip('" '), retval)
        if post:
            return (trees, robjects.r("tfgy[[5]]"))
        else:
            return trees

    def init_stages_model (self):
        """
        Defines an ODE system (by string expressions) for an SIR model where
        the infected class moves through three stages: acute, asymptomatic, and chronic.
        :return:
        """
        robjects.r("demes <- c('I1', 'I2', 'I3')")

        # transition from susceptible by infection (from any stage) to stage one only
        robjects.r("births <- rbind(c('parms$beta1*S*I1 / (S+I1+I2+I3)', '0', '0'), "
                   "c('parms$beta2*S*I2 / (S+I1+I2+I3)', '0', '0'), "
                   "c('parms$beta3*S*I3 / (S+I1+I2+I3)', '0', '0'))")
        robjects.r("rownames(births)=colnames(births) <- demes")

        # transition between stages of infection by "migration"
        robjects.r("migrations <- rbind(c('0', 'parms$alpha1 * I1', '0'), "
                   "c('0', '0', 'parms$alpha2 * I2'),"
                   "c('0', '0', '0'))")
        robjects.r("rownames(migrations)=colnames(migrations) <- demes")

        # assume that increased death rate only at final stage
        robjects.r("deaths <- c('(parms$mu)*I1', '(parms$mu)*I2', '(parms$mu+parms$gamma)*I3')")
        robjects.r("names(deaths) <- demes")

        # dynamics of susceptible class (replacement of deaths, loss to infection)
        robjects.r("nonDemeDynamics <- paste(sep='', '-parms$mu*S + parms$mu*(S+I1+I2) + (parms$mu+parms$gamma)*I3', "
                   "'-S*(parms$beta1*I1 + parms$beta2*I2 + parms$beta3*I3) / (S+I1+I2+I3)')")
        robjects.r("names(nonDemeDynamics) <- 'S'")

    def simulate_stages_trees(self, params, tree_height, tip_heights, post=False):
        """

        :return:
        """
        # set parameters
        robjects.r('N=%f; beta1=%f; beta2=%f; beta3=%f' % (params['N'], params['beta1'], params['beta2'],
                                                           params['beta3']))
        robjects.r('alpha1=%f; alpha2=%f' % (params['alpha1'], params['alpha2']))
        robjects.r('gamma=%f; mu=%f' % (params['gamma'], params['mu']))
        robjects.r('t_end=%f' % (tree_height,))

        # update model parameters
        robjects.r("S=N-1; I1=1; I2=0; I3=0")
        robjects.r("x0 <- c(I1=I1, I2=I2, I3=I3, S=S)")
        robjects.r("parms <- list(beta1=beta1, beta2=beta2, beta3=beta3, alpha1=alpha1, alpha2=alpha2, "
                   "gamma=gamma, mu=mu)")

        robjects.r("n.tips <- %d" % len(tip_heights))
        robjects.r("tip.heights <- c(%s)" % ','.join(map(str, tip_heights)))

        robjects.r("sampleTimes <- t_end - tip.heights")
        robjects.r("sampleStates <- matrix(1, nrow=n.tips, ncol=length(demes))")
        robjects.r("colnames(sampleStates) <- demes")
        robjects.r("rownames(sampleStates) <- 1:n.tips")

        robjects.r("m <- nrow(births)")
        robjects.r("maxSampleTime <- max(sampleTimes)")

        # solve ODE
        robjects.r("tfgy <- make.fgy( t0, maxSampleTime, births, deaths, nonDemeDynamics, x0, migrations=migrations, "
                   "parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod)")

        # use prevalence of respective infected classes at end of simulation to determine sample states
        robjects.r("demes.t.end <- tfgy[[4]][[1]]")
        if robjects.r("sum(demes.t.end)")[0] < len(tip_heights):
            # number of infected individuals at end of simulation is less than number of tips
            return []

        try:
            robjects.r("demes.sample <- sample(rep(1:length(demes), times=round(demes.t.end)), size=n.tips)")
        except:
            print robjects.r('tfgy[[4]]')
            raise

        robjects.r("sampleStates <- matrix(0, nrow=n.tips, ncol=length(demes))")
        robjects.r("colnames(sampleStates) <- demes")
        robjects.r("for (i in 1:n.tips) { sampleStates[i, demes.sample[i]] <- 1 }")
        robjects.r("rownames(sampleStates) <- paste(1:n.tips, demes.sample, sep='_')")

        # simulate trees
        try:
            robjects.r("trees <- simulate.binary.dated.tree.fgy( tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], "
                       "sampleTimes, sampleStates, integrationMethod = integrationMethod, "
                       "n.reps=nreps, cluster=cl)")
        except:
            return []

        # convert R objects into Python strings in Newick format
        robjects.r("class(trees) <- 'multiPhylo'")
        try:
            retval = robjects.r("lapply(trees, write.tree)")
        except:
            # error converting trees
            return []

        trees = map(lambda x: str(x).split()[-1].strip('" '), retval)
        if post:
            return (trees, robjects.r("tfgy[[5]]"))
        else:
            return trees


    def init_pangea (self):
        """
        Bespoke model for PANGEA model comparison exercise.
        :return:
        """

        # I = acute, J = asymptomatic, K = chronic; S = susceptible
        # 0 = source population, 1 = main regional population, 2 = high risk minority population
        robjects.r("demes <- c('I0', 'J0', 'K0', 'I1', 'J1', 'K1', 'I2', 'J2', 'K2')")

        # mixing rates - convenient shorthand for birth rate matrix
        robjects.r("p11 <- '((1-parms$rho1) * parms$c1*(S1+I1+J1+K1)*(1-parms$rho1) / (parms$c1*(S1+I1+J1+K1)*(1-parms$rho1) + parms$c2*(S2+I2+J2+K2)*(1-parms$rho2)))'")
        robjects.r("p12 <- '(parms$rho1 + (1-parms$rho1) * parms$c2*(S2+I2+J2+K2)*(1-parms$rho2) / (parms$c1*(S1+I1+J1+K1)*(1-parms$rho1) + parms$c2*(S2+I2+J2+K2)*(1-parms$rho2)))'")
        robjects.r("p21 <- '(parms$rho2 + (1-parms$rho2) * parms$c1*(S1+I1+J1+K1)*(1-parms$rho1) / (parms$c1*(S1+I1+J1+K1)*(1-parms$rho1) + parms$c2*(S2+I2+J2+K2)*(1-parms$rho2)))'")
        robjects.r("p22 <- '((1-parms$rho2) * parms$c2*(S2+I2+J2+K2)*(1-parms$rho2) / (parms$c1*(S1+I1+J1+K1)*(1-parms$rho1) + parms$c2*(S2+I2+J2+K2)*(1-parms$rho2)))'")

        # birth rates
        robjects.r("births <- rbind(c('parms$beta1*parms$c1*I0/(S0+I0+J0+K0)*S0', '0', '0', "
                   "'0', '0', '0',"
                   "'0', '0', '0'),"
                   "c('parms$beta2*parms$c1*J0/(S0+I0+J0+K0)*S0', '0', '0',"
                   "'0', '0', '0',"
                   "'0', '0', '0'),"
                   "c('parms$beta3*parms$c1*K0/(S0+I0+J0+K0)*S0', '0', '0',"
                   "'0', '0', '0',"
                   "'0', '0', '0'),"
                   "c('0', '0', '0',"
                   "paste(sep='*', 'parms$z*parms$beta1*parms$c1', p11, 'I1/(S1+I1+J1+K1)*S1'), '0', '0',"
                   "paste(sep='*', 'parms$z*parms$beta1*parms$c2', p21, 'I1/(S1+I1+J1+K1)*S2'), '0', '0'),"
                   "c('0', '0', '0',"
                   "paste(sep='*', 'parms$z*parms$beta2*parms$c1', p11, 'J1/(S1+I1+J1+K1)*S1'), '0', '0',"
                   "paste(sep='*', 'parms$z*parms$beta2*parms$c2', p21, 'J1/(S1+I1+J1+K1)*S2'), '0', '0'),"
                   "c('0', '0', '0',"
                   "paste(sep='*', 'parms$z*parms$beta3*parms$c1', p11, 'K1/(S1+I1+J1+K1)*S1'), '0', '0',"
                   "paste(sep='*', 'parms$z*parms$beta3*parms$c2', p21, 'K1/(S1+I1+J1+K1)*S2'), '0', '0'),"
                   "c('0', '0', '0',"
                   "paste(sep='*', 'parms$z*parms$beta1*parms$c1', p12, 'I2/(S2+I2+J2+K2)*S1'), '0', '0',"
                   "paste(sep='*', 'parms$z*parms$beta1*parms$c2', p22, 'I2/(S2+I2+J2+K2)*S2'), '0', '0'),"
                   "c('0', '0', '0',"
                   "paste(sep='*', 'parms$z*parms$beta2*parms$c1', p12, 'J2/(S2+I2+J2+K2)*S1'), '0', '0',"
                   "paste(sep='*', 'parms$z*parms$beta2*parms$c2', p22, 'J2/(S2+I2+J2+K2)*S2'), '0', '0'),"
                   "c('0', '0', '0',"
                   "paste(sep='*', 'parms$z*parms$beta3*parms$c1', p12, 'K2/(S2+I2+J2+K2)*S1'), '0', '0',"
                   "paste(sep='*', 'parms$z*parms$beta3*parms$c2', p22, 'K2/(S2+I2+J2+K2)*S2'), '0', '0'))")
        robjects.r("rownames(births)=colnames(births) <- demes")

        # migration rates (transition from acute to chronic, not spatial migration)
        robjects.r("migrations <- rbind("
                   "c('0', 'parms$alpha1*I0', '0', '(parms$mig*I1 + 1e-4) * I0/(S0+I0+J0+K0)', '0', '0', '0', '0', '0'),"
                   "c('0', '0', 'parms$alpha2*J0', '0', 'parms$mig*J1*J0/(S0+I0+J0+K0)', '0', '0', '0', '0'),"
                   "c('0', '0', '0', '0', '0', '0', '0', '0', '0'),"
                   "c('0', '0', '0', '0', 'parms$alpha1*I1', '0', '0', '0', '0'),"
                   "c('0', '0', '0', '0', '0', 'parms$alpha2*J1', '0', '0', '0'),"
                   "c('0', '0', '0', '0', '0', '0', '0', '0', '0'),"
                   "c('0', '0', '0', '0', '0', '0', '0', 'parms$alpha1*I2', '0'),"
                   "c('0', '0', '0', '0', '0', '0', '0', '0', 'parms$alpha2*J2'),"
                   "c('0', '0', '0', '0', '0', '0', '0', '0', '0'))")
        robjects.r("rownames(migrations)=colnames(migrations) <- demes")

        # death rates (apply migrations from source population here)
        robjects.r("deaths <- c('parms$mu*I0', 'parms$mu*J0', '(parms$mu+parms$gamma)*K0',"
                   "'parms$mu*I1', 'parms$mu*J1', '(parms$mu+parms$gamma)*K1',"
                   "'parms$mu*I2', 'parms$mu*J2', '(parms$mu+parms$gamma)*K2')")
        robjects.r("names(deaths) <- demes")

        # dynamics of susceptible class
        robjects.r("nonDemeDynamics <- c('-parms$mu*S0 + parms$lam*S0 + parms$mu*(S0+I0+J0+K0) + parms$gamma*K0 "
                   "- S0*parms$c1*parms$z*(parms$beta1*I0+parms$beta2*J0+parms$beta3*K0)/(S0+I0+J0+K0) "
                   "+ parms$mig*(I1*I0+J1*J0)/(S0+I0+J0+K0)', "
                   "paste(sep='', '-parms$mu*S1 + parms$lam*S1 + parms$mu*(S1+I1+J1+K1) + parms$gamma*K1 "
                   "- S1*parms$c1*parms$z*(', p11, "
                   "'*(parms$beta1*I1 + parms$beta2*J1 + parms$beta3*K1)/(S1+I1+J1+K1) + ', p12, "
                   "'*(parms$beta1*I2 + parms$beta2*J2 + parms$beta3*K2)/(S2+I2+J2+K2))'),"
                   "paste(sep='', '-parms$mu*S2 + parms$lam*S2 + parms$mu*(S2+I2+J2+K2) + parms$gamma*K2 "
                   "- S2*parms$c2*parms$z*(', p21, "
                   "'*(parms$beta1*I1 + parms$beta2*J1 + parms$beta3*K1)/(S1+I1+J1+K1) + ', p22, "
                   "'*(parms$beta1*I2 + parms$beta2*J2 + parms$beta3*K2)/(S2+I2+J2+K2))'))")
        robjects.r("names(nonDemeDynamics) <- c('S0', 'S1', 'S2')")


    def simulate_pangea(self, params, tree_height, tip_heights, eval_period=5*52.):
        """
        Time scaled to weeks.
        :param params: Dictionary with model parameter settings.
        :param eval_period: For regional simulations, the last five years.
        :param tree_height: Should be about 40 years.
        :param tip_heights:
        :return:
        """

        vars = [
            'N0',  # size of source population
            'N1',  # size of target population
            'p',  # proportion of target population in high-risk group
            'alpha1',  # rate of transition from acute stage
            'alpha2',  # rate of progression to AIDS
            'lam',  # growth rate of populations
            'gamma',  # excess mortality rate due to infection
            'mu',  # baseline mortality rate
            'beta1',  # transmission rate from acute infections
            'beta2',  # transmission rate from asymptomatic
            'beta3',  # transmission rate from AIDS
            'rho1',  # proportion of contacts reserved for outside group
            'rho2',
            'c1',  # contact rate for group 1
            'c2',  # contact rate for group 2
            'mig',  # migration rate
            'z'  # factor for adjusting transmission rates over time
        ]

        # set model parameters
        for v in vars:
            if v == 'z':
                continue
            robjects.r("%s=%f" % (v, params[v]))
        robjects.r("z <- 1")
        robjects.r("parms <- list(%s)" % (','.join(['%s=%s'%(v, v) for v in vars]), ))

        # initialize population sizes
        robjects.r("S0=N0-1; I0=1; J0=0; K0=0")
        robjects.r("S1=N1*(1-p); I1=0; J1=0; K1=0")
        robjects.r("S2=N1*p; I2=0; J2=0; K2=0")
        robjects.r("x0 = c(I0=I0, J0=J0, K0=K0,"
                   "I1=I1, J1=J1, K1=K1,"
                   "I2=I2, J2=J2, K2=K2, S0=S0, S1=S1, S2=S2)")

        # set simulation conditions
        robjects.r("n.tips=%d" % (len(tip_heights), ))
        robjects.r("t.end=%f; eval.period=%f" % (tree_height, eval_period))
        robjects.r("sampleTimes <- rep(t.end, times=n.tips)")
        robjects.r("maxSampleTime <- max(sampleTimes)")

        # solve ODE system for first time interval
        robjects.r("t0 <- 0; t1 <- (maxSampleTime-eval.period) / 2")
        robjects.r("fgyRes.1 <- round(fgyResolution * (t1-t0) / maxSampleTime)")
        try:
            robjects.r("tfgy.1 <- make.fgy(t0, t1, births, deaths, nonDemeDynamics, x0, migrations=migrations, parms=parms,"
                       "fgyResolution=fgyRes.1, integrationMethod=integrationMethod)")
        except:
            sys.stderr.write("Failed to solve first ODE system\n")
            return []


        # use state of system at end of time interval to initialize next time interval
        robjects.r("x1 <- tfgy.1[[5]][fgyRes.1, 2:ncol(tfgy.1[[5]])]")
        robjects.r("t2 <- maxSampleTime-eval.period")
        robjects.r("fgyRes.2 <- round(fgyResolution * (t2-t1) / maxSampleTime)")
        robjects.r("parms$z <- %f" % (params['z1'],))
        try:
            robjects.r("tfgy.2 <- make.fgy(t1, t2, births, deaths, nonDemeDynamics, x1, migrations=migrations, parms=parms,"
                       "fgyResolution=fgyRes.2, integrationMethod=integrationMethod)")
        except:
            sys.stderr.write("Failed to solve second ODE system\n")
            return []


        # solve last time interval
        robjects.r("x2 <- tfgy.2[[5]][fgyRes.2, 2:ncol(tfgy.2[[5]])]")
        robjects.r("t3 <- maxSampleTime")
        robjects.r("fgyRes.3 <- fgyResolution - fgyRes.1 - fgyRes.2")
        robjects.r("parms$z <- %f" % (params['z2'],))
        try:
            robjects.r("tfgy.3 <- make.fgy(t2, t3, births, deaths, nonDemeDynamics, x2, migrations=migrations, parms=parms,"
                       "fgyResolution=fgyRes.3, integrationMethod=integrationMethod)")
        except:
            sys.stderr.write("Failed to solve third ODE system\n")
            return []

        # use prevalence of respective infected classes at end of simulation to determine sample states
        robjects.r("demes.t.end <- tfgy.3[[4]][[1]]")
        robjects.r("sampled.demes <- which(!grepl('0$', names(demes.t.end)))")
        if robjects.r("sum(demes.t.end[sampled.demes])")[0] < len(tip_heights):
            # number of infected individuals at end of simulation is less than number of tips
            sys.stderr.write("At end of simulation, {} individuals infected, but {} tips in tree\n".format(robjects.r("sum(demes.t.end[sampled.demes])")[0], len(tip_heights)))
            return []

        try:
            robjects.r("demes.sample <- sample(rep(sampled.demes, times=round(demes.t.end[sampled.demes])), size=n.tips)")
        except:
            print 'demes.t.end', robjects.r("demes.t.end")
            print 'n.tips', robjects.r('n.tips')
            raise

        # sample demes based on frequencies at last time point
        robjects.r("sampleStates <- matrix(0, nrow=n.tips, ncol=length(demes))")
        robjects.r("colnames(sampleStates) <- demes")
        robjects.r("for (i in 1:n.tips) { sampleStates[i, demes.sample[i]] <- 1 }")
        robjects.r("rownames(sampleStates) <- paste(1:n.tips, demes.sample, sep='_')")

        # reconstritute the entire ODE
        robjects.r("y.times <- rev(seq(0, maxSampleTime, length.out=fgyResolution))")
        robjects.r("y.births <- c(tfgy.3[[2]], tfgy.2[[2]], tfgy.1[[2]])")
        robjects.r("y.migrations <- c(tfgy.3[[3]], tfgy.2[[3]], tfgy.1[[3]])")
        robjects.r("y.demeSizes <- c(tfgy.3[[4]], tfgy.2[[4]], tfgy.1[[4]])")

        # simulate trees
        try:
            robjects.r("trees <- simulate.binary.dated.tree.fgy( y.times, y.births, y.migrations, y.demeSizes, sampleTimes,"
                       " sampleStates, integrationMethod = integrationMethod, n.reps=nreps)")
        except:
            sys.stderr.write("Failed to simulate trees\n")
            return []

        # convert R objects into Python strings in Newick format
        robjects.r("class(trees) <- 'multiPhylo'")
        try:
            retval = robjects.r("lapply(trees, write.tree)")
        except:
            sys.stderr.write("Failed to write trees\n")
            # error converting trees
            return []

        trees = map(lambda x: str(x).split()[-1].strip('" '), retval)
        return trees
