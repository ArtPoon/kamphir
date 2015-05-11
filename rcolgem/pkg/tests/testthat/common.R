set.seed(0)

# common parameters for several FGY tests
si.t0 <- 0
si.t1 <- 100
n.tips <- 5
sampleTimes <- rep(si.t1, n.tips)
integrationMethod <- 'rk4'

# initial conditions
si.x0 <- c(I=1, S=9)
diffrisk.x0 <- c(I1=1, I2=1, S1=9, S2=9)

# model parameters
si.parms <- list(beta=0.1, gamma=0.005, mu=0.005)
diffrisk.parms <- c(si.parms, c1=0.5, c2=1, rho=0.9)

# common things to be tested for all FGY's
common.fgy.test <- function (fgy) {
    expect_equal(fgy[[1]], seq(100, 0, -10))
    expect_equal(fgy[[5]][,"time"], seq(0, 100, 10))
}

# set up and solve a simple SI model
setup.si.fgy <- function (t0=si.t0, t1=si.t1, x0=si.x0, parms=si.parms, fgyResolution=11) {
    # set up a simple SI model
    births <- matrix("parms$beta*S*I / (S+I)", dimnames=list("I", "I"))
    deaths <- c(I="(parms$mu+parms$gamma)*I")
    nonDemeDynamics <- c(S="(parms$mu+parms$gamma)*I - parms$beta*I*S / (S+I)")

    # no migration
    migrations <- NA

    make.fgy(t0, t1, births, deaths, nonDemeDynamics, x0, migrations, parms, fgyResolution, integrationMethod)
}

# set up and solve a two-phase SI model (infectiousness changes once, all other
# parameters remain constant
setup.si2.fgy <- function () {
    fgy1 <- setup.si.fgy(t1=40, fgyResolution=5)

    parms2 <- si.parms
    parms2[["beta"]] <- 0.2
    x0.2 <- fgy1[[5]][5, 2:3]

    fgy2 <- setup.si.fgy(t0=50, x0=x0.2, parms=parms2, fgyResolution=6)
    list(c(fgy2[[1]], fgy1[[1]]),
         c(fgy2[[2]], fgy1[[2]]),
         c(fgy2[[3]], fgy1[[3]]),
         c(fgy2[[4]], fgy1[[4]]),
         rbind(fgy1[[5]], fgy2[[5]]))
}

# set up and solve a differential risk model (two groups with preferential
# contact and different contact rates)
setup.diffrisk.fgy <- function (t0=si.t0, t1=si.t1, x0=diffrisk.x0, parms=diffrisk.parms, fgyResolution=11) {
    demes <- c("I1", "I2")
    p.denom <- "(parms$c1 * (S1 + I1) + parms$c2 * (S2 + I2))"
    p11 <- paste0("(parms$rho + (1-parms$rho) * parms$c1 * (S1 + I1) / ", p.denom, ")")
    p12 <- paste0("(1-parms$rho) * parms$c2 * (S2 + I2) / ", p.denom)
    p21 <- paste0("(1-parms$rho) * parms$c1 * (S1 + I1) / ", p.denom)
    p22 <- paste0("(parms$rho + (1-parms$rho) * parms$c2 * (S2 + I2) / ", p.denom, ")")

    births <- rbind(c(paste("parms$beta * parms$c1", p11, "I1 * S1 / (S1 + I1)", sep="*"),
                      paste("parms$beta * parms$c2", p21, "I1 * S2 / (S1 + I1)", sep="*")),
                    c(paste("parms$beta * parms$c1", p12, "I2 * S1 / (S2 + I2)", sep="*"),
                      paste("parms$beta * parms$c2", p22, "I2 * S2 / (S2 + I2)", sep="*")))
    dimnames(births) <- list(demes, demes)
    migrations <- matrix(0, nrow=2, ncol=2, dimnames=list(demes, demes))
    deaths <- setNames(paste("(parms$mu + parms$gamma) *", demes), demes)

    nonDemeDynamics <- c(S1=paste0("(parms$mu + parms$gamma) * I1 - ",
                                   "parms$beta * parms$c1 * S1 * ",
                                   "(", p11, "* I1 / (S1 + I1) + ", p12, "* I2 / (S2 + I2))"),
                         S2=paste0("(parms$mu + parms$gamma) * I2 - ",
                                   "parms$beta * parms$c2 * S2 * ",
                                   "(", p21, "* I1 / (S1 + I1) + ", p22, "* I2 / (S2 + I2))"))
    make.fgy(t0, t1, births, deaths, nonDemeDynamics, x0, migrations, parms, fgyResolution, integrationMethod)
}
