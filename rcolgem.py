import rpy2.robjects as robjects  # R is instantiated upon load module

class Rcolgem ():
    def __init__ (self, ncores, nreps, t0=0, fgy_resolution=500., integration_method='adams'):
        # default settings
        robjects.globalenv['n.cores'] = ncores
        robjects.globalenv['nreps'] = nreps
        robjects.globalenv['fgyResolution'] = fgy_resolution
        robjects.globalenv['integrationMethod'] = integration_method
        robjects.globalenv['t0'] = t0


    def init_SI_model (self, N=1000, beta=0.01, gamma=1/520., mu=1/3640.):
        """
        Defines a susceptible-infected-recovered model in rcolgem.
        :param N: total population size
        :param beta: transmission rate
        :param gamma: excess mortality rate of infected individuals
        :param mu: baseline mortality rate
        :return:
        """
        robjects.r('N=%f; beta=%f; gamma=%f; mu=%f' % (N, beta, gamma, mu))
        robjects.r('S = N-1')
        robjects.r('I = 1')
        robjects.r('x0 <- c(I=I, S=S)')
        robjects.r('parms <- list(beta=beta, gamma=gamma, mu=mu)')

        # define ODE system
        robjects.r('demes <- c("I")')

        robjects.r('births <- rbind(c("parms$beta*S*I / (S+I)"))')
        robjects.r('rownames(births) <- colnames(births) <- demes')

        robjects.r('migrations <- rbind(c("0"))')
        robjects.r('rownames(migrations)=colnames(migrations) <- demes')

        robjects.r("deaths <- c('(parms$mu+parms$gamma)*I')")
        robjects.r('names(deaths) <- demes')

        robjects.r("nonDemeDynamics <- paste(sep='', '-parms$mu*S + parms$mu*S + "
                   "(parms$mu+parms$gamma)*I', '-S*(parms$beta*I) / (S+I)')")
        robjects.r("names(nonDemeDynamics) <- 'S'")


    def simulate_SI_trees (self, params, tip_heights):
        """
        Simulate coalescent trees under the SI model.
        :param tip_heights:
        :return:
        """

        # set parameters
        robjects.r('N=%f; beta=%f; gamma=%f; mu=%f; t_end=%f' % (params['N'], params['beta'], params['gamma'],
                                                                 params['mu'], params['t_end']))
        robjects.r('S = N-1')
        robjects.r('I = 1')
        robjects.r('x0 <- c(I=I, S=S)')
        robjects.r('parms <- list(beta=beta, gamma=gamma, mu=mu)')

        robjects.r("n.tips <- %d" % len(tip_heights))
        robjects.r("tip.heights <- c(%s)" % ','.join(map(str, tip_heights)))

        robjects.r("sampleTimes <- t_end - tip.heights")
        robjects.r("sampleStates <- matrix(1, nrow=n.tips, ncol=length(demes))")
        robjects.r("colnames(sampleStates) <- demes")
        robjects.r("rownames(sampleStates) <- 1:n.tips")

        robjects.r("m <- nrow(births)")
        robjects.r("maxSampleTime <- max(sampleTimes)")
        robjects.r("require(rcolgem, quietly=TRUE)")

        # solve ODE
        robjects.r("tfgy <- make.fgy( t0, maxSampleTime, births, deaths, nonDemeDynamics, x0, migrations=migrations, "
                   "parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod )")

        # simulate trees
        robjects.r("trees <- simulate.binary.dated.tree.fgy( tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, "
                   "sampleStates, integrationMethod = integrationMethod, n.reps=nreps)")
        robjects.r("'multiPhylo' -> class(trees)")
        retval = robjects.r("lapply(trees, write.tree)")
        trees = map(lambda x: str(x).split()[-1].strip('" '), retval)
        return trees


