context("FGY")

source(file="common.R")

test_that("A simple SI trajectory is correct", { 
    si.fgy <- setup.si.fgy()
    common.fgy.test(si.fgy)

    # number of births at time t, I think
    expect_equal(unlist(si.fgy[[2]]), 
                 c(0.09074935, 0.09181836, 0.09439305, 0.10050123, 0.11447158,
                   0.14372733, 0.19331960, 0.24280038, 0.23725521, 0.16673089,
                   0.09000000))

    # number of migrations should be zero
    expect_equal(unlist(si.fgy[[3]]), rep(0, 11))

    # number of infected at time t
    expect_equal(unname(unlist(si.fgy[[4]])), 
                 c(8.990622140, 8.977205574, 8.944704626, 8.866507095,
                   8.681418495, 8.259948941, 7.380764495, 5.848505461,
                   3.871071724, 2.114361265, 1.000000000))

    # the last result should be an aggregation of the other results
    expect_equal(si.fgy[[5]][,"I"], unname(unlist(rev(si.fgy[[4]]))))
    expect_equal(si.fgy[[5]][,"S"], 10-unname(unlist(rev(si.fgy[[4]]))))
})

test_that("An SI model stays at equilibrium", {
    eq.parms <- list(beta=0.02, gamma = 0.009, mu = 0.009)
    eq.si.fgy <- setup.si.fgy(parms=eq.parms)
    common.fgy.test(eq.si.fgy)

    # at equilibrium, gamma + mu people die and are born at every time step
    expect_equal(unlist(eq.si.fgy[[2]]), rep(0.018, 11))

    # no migration
    expect_equal(unlist(eq.si.fgy[[3]]), rep(0, 11))

    # no change in number of infected people
    expect_equal(unname(unlist(eq.si.fgy[[4]])), rep(1, 11))
})

test_that("A two-phase SI model works", {
    si2.fgy <- setup.si2.fgy()
    common.fgy.test(si2.fgy)

    expect_equal(unlist(si2.fgy[[2]]),
                 c(0.095849145, 0.097782551, 0.104096514, 0.124533779,
                   0.189451594, 0.386639208, 0.193319604, 0.242800385,
                   0.237255209, 0.166730891, 0.090000000))
    expect_equal(unlist(si2.fgy[[3]]), rep(0, 11))
    expect_equal(unname(unlist(si2.fgy[[4]])),
                 c(9.4952801, 9.4845147, 9.4491768, 9.3328179, 8.9404848,
                   7.3807645, 7.3807645, 5.8485055, 3.8710717, 2.1143613,
                   1.0000000))
    expect_equal(si2.fgy[[5]][,"I"], unname(unlist(rev(si2.fgy[[4]]))))
    expect_equal(si2.fgy[[5]][,"S"], 10-unname(unlist(rev(si2.fgy[[4]]))))
})

test_that("A differential risk SI model works", {
    diffrisk.fgy <- setup.diffrisk.fgy()
    common.fgy.test(diffrisk.fgy)

    births <- t(sapply(diffrisk.fgy[[2]], c))
    expect_equal(births,
                 matrix(c(0.092093939, 0.0080993716, 0.0024762417, 0.088417726,
                          0.098532934, 0.0090511164, 0.0024046609, 0.089681035,
                          0.105747318, 0.0103272443, 0.0023341814, 0.092549853,
                          0.112439532, 0.0119296630, 0.0023001664, 0.099081923,
                          0.116405818, 0.0137201967, 0.0023726394, 0.113538569,
                          0.115049662, 0.0152705810, 0.0026507957, 0.142847452,
                          0.106689411, 0.0157173808, 0.0031880130, 0.190679924,
                          0.091992754, 0.0140309747, 0.0038108775, 0.235985457,
                          0.074089564, 0.0102280539, 0.0040739272, 0.228336465,
                          0.056687426, 0.0060224761, 0.0037238632, 0.160623136,
                          0.042000000, 0.0030000000, 0.0030000000, 0.087000000),
                        ncol=4, byrow=TRUE))

    migrations <- c(sapply(diffrisk.fgy[[3]], c))
    expect_equal(migrations, rep(0, 44))

    demes <- as.data.frame(do.call(rbind, diffrisk.fgy[[4]]))
    expected.demes <- data.frame(I1=c(7.29468615, 6.97124322, 6.52966020,
                                      5.95174287, 5.23642379, 4.41135658,
                                      3.53781556, 2.70059420, 1.97945950,
                                      1.41493454, 1.00000000),
                                 I2=c(8.98162513, 8.96517988, 8.92757906,
                                      8.84059186, 8.64069098, 8.19729219,
                                      7.29662592, 5.76662337, 3.82569751,
                                      2.10451842, 1.00000000))
    expect_equal(demes, expected.demes)
})

test_that("We can't make an FGY with bad intial conditions", {
    births <- matrix("parms$beta*S*I / (S+I)", dimnames=list("I", "I"))
    deaths <- c(I="(parms$mu+parms$gamma)*I")
    nonDemeDynamics <- c(S="(parms$mu+parms$gamma)*I - parms$beta*I*S / (S+I)")
    migrations <- NA

    # try the wrong number of demes
    x0 <- c(I1=1, I2=1, S1=49, S2=49)
    expect_error(make.fgy(0, 100, births, deaths, nonDemeDynamics, x0, parms=si.parms))

    # try right number but wrong names
    x0 <- c(I1=1, S1=1)
    expect_error(make.fgy(0, 100, births, deaths, nonDemeDynamics, x0, parms=si.parms))
})
