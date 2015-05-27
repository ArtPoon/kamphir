context("Tree")

# It's important that the first line in any stochastic test be a call to
# "set.seed", otherwise the order in which the tests are run will influence the
# results.

options(scipen=10)

test_that("A tree can be simulated from an SI model", {
    set.seed(0)
    si.fgy <- setup.si.fgy()
    sampleStates <- matrix(1, nrow=n.tips, ncol=1, dimnames=list(1:n.tips, "I"))
    tree <- simulate.binary.dated.tree.fgy(si.fgy[[1]], si.fgy[[2]], si.fgy[[3]], si.fgy[[4]],
                                           sampleTimes, sampleStates, integrationMethod)[[1]]

    # this is the tree I got from running it the first time
    # I don't know whether it's "correct" or not, this is just a test for
    # consistency
    treestr <- "(4:95.82816554,(3:92.33793996,(5:78.61132095,(1:37.77986406,2:37.77986406):40.83145689):13.72661901):3.490225578);"
    expected.tree <- read.tree(text=treestr)
    expect_equal(tree, expected.tree)
})

test_that("A tree can be quickly simulated from an SI model", {
    set.seed(0)
    si.fgy <- setup.si.fgy()
    sampleStates <- matrix(1, nrow=n.tips, ncol=1, dimnames=list(1:n.tips, "I"))
    tree <- simulate.bdt.onedeme(si.fgy[[1]], si.fgy[[2]], si.fgy[[4]],
                                 sampleTimes, sampleStates,
                                 integrationMethod)[[1]]

    treestr <- "(4:95.82816554,(3:92.33793996,(5:78.61132095,(1:37.77986406,2:37.77986406):40.83145689):13.72661901):3.490225578);"
    expected.tree <- read.tree(text=treestr)
    expect_equal(tree, expected.tree)
})

test_that("A tree can be simulated with an SI2 model", {
    set.seed(0)
    si2.fgy <- setup.si2.fgy()
    sampleStates <- matrix(1, nrow=n.tips, ncol=1, dimnames=list(1:n.tips, "I"))
    tree <- simulate.binary.dated.tree.fgy(si2.fgy[[1]], si2.fgy[[2]], si2.fgy[[3]], si2.fgy[[4]],
                                           sampleTimes, sampleStates, integrationMethod)[[1]]
    treestr <- "(4:95.62881561,(3:91.9718106,(5:77.15114683,(1:37.10303152,2:37.10303152):40.04811531):14.82066377):3.657005013);"
    expected.tree <- read.tree(text=treestr)
    expect_equal(tree, expected.tree)
})

test_that("A tree can be simulated quickly with an SI2 model", {
    set.seed(0)
    si2.fgy <- setup.si2.fgy()
    sampleStates <- matrix(1, nrow=n.tips, ncol=1, dimnames=list(1:n.tips, "I"))
    tree <- simulate.bdt.onedeme(si2.fgy[[1]], si2.fgy[[2]], si2.fgy[[4]],
                                 sampleTimes, sampleStates,
                                 integrationMethod)[[1]]
    treestr <- "(4:95.62881561,(3:91.9718106,(5:77.15114683,(1:37.10303152,2:37.10303152):40.04811531):14.82066377):3.657005013);"
    expected.tree <- read.tree(text=treestr)
    expect_equal(tree, expected.tree)
})

test_that("A tree can be simulated with a DiffRisk model", {
    set.seed(0)
    diffrisk.fgy <- setup.diffrisk.fgy()
    sampleStates <- matrix(rep(c(1, 0, 0, 1), each=5), ncol=2, dimnames=list(1:10, c("I1", "I2")))
    sampleTimes <- rep(100, 10)
    tree <- simulate.binary.dated.tree.fgy(diffrisk.fgy[[1]], diffrisk.fgy[[2]], diffrisk.fgy[[3]], diffrisk.fgy[[4]],
                                           sampleTimes, sampleStates, integrationMethod)[[1]]
    treestr <- "((3:76.32224589,(4:50.12690607,(2:23.6530747,(1:8.742420893,5:8.742420893):14.91065381):26.47383137):26.19533982):22.8438812,((8:60.08119625,9:60.08119625):19.88816683,(7:74.18552162,(6:27.5078683,10:27.5078683):46.67765332):5.783841456):19.19676402);"
    expected.tree <- read.tree(text=treestr)
    expect_equal(tree, expected.tree)
})

test_that("A tree can be simulated with migrations", {
    set.seed(0)
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
    deaths <- setNames(paste("(parms$mu + parms$gamma) *", demes), demes)

    nonDemeDynamics <- c(S1=paste0("(parms$mu + parms$gamma) * I1 - ",
                                   "parms$beta * parms$c1 * S1 * ",
                                   "(", p11, "* I1 / (S1 + I1) + ", p12, "* I2 / (S2 + I2))"),
                         S2=paste0("(parms$mu + parms$gamma) * I2 - ",
                                   "parms$beta * parms$c2 * S2 * ",
                                   "(", p21, "* I1 / (S1 + I1) + ", p22, "* I2 / (S2 + I2))"))
    x0 <- diffrisk.x0
    parms <- diffrisk.parms

    migrations <- matrix(0, nrow=2, ncol=2, dimnames=list(demes, demes))
    no.mig.fgy <- make.fgy(0, 100, births, deaths, nonDemeDynamics, x0,
                           migrations, parms, 11, 'rk4')

    migrations <- matrix(c("0", "0", "parms$alpha*I1", "0"), nrow=2, ncol=2, dimnames=list(demes, demes))
    parms <- c(parms, alpha=0.01)
    mig.fgy <- make.fgy(0, 100, births, deaths, nonDemeDynamics, x0, migrations, parms,
                        11, 'rk4')

    sampleStates <- matrix(rep(c(1, 0, 0, 1), each=5), ncol=2, dimnames=list(1:10, c("I1", "I2")))
    sampleTimes <- rep(100, 10)

    set.seed(0)
    no.mig.tree <- simulate.binary.dated.tree.fgy(no.mig.fgy[[1]], no.mig.fgy[[2]], no.mig.fgy[[3]], no.mig.fgy[[4]],
                                           sampleTimes, sampleStates, integrationMethod)[[1]]
    expect_equal(write.tree(no.mig.tree),
                 "((3:76.32224589,(4:50.12690607,(2:23.6530747,(1:8.742420893,5:8.742420893):14.91065381):26.47383137):26.19533982):22.8438812,((8:60.08119625,9:60.08119625):19.88816683,(7:74.18552162,(6:27.5078683,10:27.5078683):46.67765332):5.783841456):19.19676402);")

    set.seed(0)
    mig.tree <- simulate.binary.dated.tree.fgy(mig.fgy[[1]], mig.fgy[[2]], mig.fgy[[3]], mig.fgy[[4]],
                                           sampleTimes, sampleStates, integrationMethod)[[1]]
    expect_equal(write.tree(mig.tree), "((7:71.02587233,9:71.02587233):28.07440856,(8:78.02379982,(6:73.61111634,(4:52.60637368,(10:38.91458526,(3:19.13767693,(2:16.41012696,(1:5.751801791,5:5.751801791):10.65832517):2.72754997):19.77690833):13.69178842):21.00474266):4.412683474):21.07648107);")
})

test_that("We are warned about uneven sample times", {
    fgy <- setup.si.fgy()
    times <- fgy[[1]] + rnorm(length(fgy[[1]]))
    sampleTimes <- rep(max(times), 5)
    sampleStates <- matrix(1, ncol=1, nrow=5, dimnames=list(1:5, "I"))
    expect_warning(
        simulate.binary.dated.tree.fgy(times, fgy[[2]], fgy[[3]], fgy[[4]],
                                       sampleTimes, sampleStates))
})

test_that("Approximately even sample times are ok", {
    set.seed(0)
    fgy <- setup.si.fgy()
    times <- fgy[[1]] + rnorm(length(fgy[[1]])) * .Machine$double.eps ^ 0.75
    sampleTimes <- rep(max(times), 5)
    sampleStates <- matrix(1, ncol=1, nrow=5, dimnames=list(1:5, "I"))
    expect_that(simulate.binary.dated.tree.fgy(times, fgy[[2]], fgy[[3]],
                                               fgy[[4]], sampleTimes,
                                               sampleStates),
                not(gives_warning()))
})

test_that("We are warned about missing row names", {
    fgy <- setup.si.fgy()
    sampleTimes <- setNames(rep(max(fgy[[1]]), 5), letters[1:5])
    sampleStates <- matrix(1, ncol=1, nrow=5)
    colnames(sampleStates) <- c("I")
    expect_warning(
        simulate.binary.dated.tree.fgy(fgy[[1]], fgy[[2]], fgy[[3]], fgy[[4]],
                                       sampleTimes, sampleStates))
})

test_that("We can't use a data.frame for sample states", {
    set.seed(0)
    fgy <- setup.si.fgy()
    sampleTimes <- rep(max(fgy[[1]]), 5)
    sampleStates <- matrix(1, ncol=1, nrow=5, dimnames=list(1:5, "I"))
    sampleStates <- as.data.frame(sampleStates)
    expect_error(
        simulate.binary.dated.tree.fgy(fgy[[1]], fgy[[2]], fgy[[3]], fgy[[4]],
                                       sampleTimes, sampleStates))
})

test_that("We can't sample past the end of the time axis", {
    fgy <- setup.si.fgy()
    sampleTimes <- c(rep(max(fgy[[1]]), 4), max(fgy[[1]]) + 1)
    sampleStates <- matrix(1, ncol=1, nrow=5, dimnames=list(1:5, "I"))
    expect_error(
        simulate.binary.dated.tree.fgy(fgy[[1]], fgy[[2]], fgy[[3]], fgy[[4]],
                                       sampleTimes, sampleStates))
})

test_that("Out of order samples are consistent", {
    set.seed(0)
    fgy <- setup.si.fgy()
    sampleTimes <- c(rep(max(fgy[[1]]), 4), median(fgy[[1]]))
    sampleStates <- matrix(1, ncol=1, nrow=5, dimnames=list(1:5, "I"))
    expect_false(all(sampleTimes == sort(sampleTimes)))
    tree <- simulate.binary.dated.tree.fgy(fgy[[1]], fgy[[2]], fgy[[3]], fgy[[4]],
                                           sampleTimes, sampleStates)[[1]]
    treestr <- "(4:96.02066402,(3:92.69148584,(5:29.87938507,(1:51.7406756,2:51.7406756):28.13870947):12.81210077):3.329178179);"
    expected.tree <- read.tree(text=treestr)
    expect_equal(tree, expected.tree)
})

test_that("Replicate samples are consistent", {
    set.seed(0)
    si.fgy <- setup.si.fgy()
    sampleStates <- matrix(1, nrow=n.tips, ncol=1, dimnames=list(1:n.tips, "I"))
    trees <- simulate.binary.dated.tree.fgy(si.fgy[[1]], si.fgy[[2]], si.fgy[[3]], si.fgy[[4]],
                                            sampleTimes, sampleStates, integrationMethod, 
                                            n.reps=2)

    treestrs <- c("(4:95.82816554,(3:92.33793996,(5:78.61132095,(1:37.77986406,2:37.77986406):40.83145689):13.72661901):3.490225578);",
                  "((3:61.05611437,(2:2.463866035,4:2.463866035):58.59224834):29.11117026,(1:61.3252033,5:61.3252033):28.84208134);")
    expected.trees <- mapply(read.tree, text=treestrs, SIMPLIFY=FALSE)
    mapply(expect_equal, trees, expected.trees)
})

test_that("Heterochronous samples are consistent", {
    set.seed(0)
    si.fgy <- setup.si.fgy()
    sampleStates <- matrix(1, nrow=n.tips, ncol=1, dimnames=list(1:n.tips, "I"))
    sampleTimes <- c(rep(si.t1, n.tips-2), round(3*si.t1/4), round(si.t1/2))
    trees <- simulate.binary.dated.tree.fgy(si.fgy[[1]], si.fgy[[2]], si.fgy[[3]], si.fgy[[4]],
                                            sampleTimes, sampleStates, integrationMethod, 
                                            n.reps=2)
    treestrs <- c("(4:71.0800247,(3:92.80050865,(5:30.2529299,(1:54.42753403,2:54.42753403):25.82539587):12.54757875):3.279516052);",  
                  "((4:40.86913106,(1:7.012593634,3:7.012593634):58.85653743):24.81543379,(2:66.09120694,5:16.09120694):24.59335791);")
    expected.trees <- mapply(read.tree, text=treestrs, SIMPLIFY=FALSE)
    mapply(expect_equal, trees, expected.trees)
})

test_that("We are warned about mismatching sample info", {
    fgy <- setup.si.fgy()
    sampleTimes <- setNames(rep(max(fgy[[1]]), 5), letters[1:5])
    sampleStates <- matrix(1, ncol=1, nrow=5)
    rownames(sampleStates) <- letters[2:6]
    colnames(sampleStates) <- c("I")
    expect_warning(
        simulate.binary.dated.tree.fgy(fgy[[1]], fgy[[2]], fgy[[3]], fgy[[4]],
                                       sampleTimes, sampleStates))
})

test_that("A simulation is consistent", {
    load("input.Rdata")
    load("output.Rdata")
    set.seed(0)
    res <- simulate.binary.dated.tree.fgy(times, births, migrations, demeSizes,
                                          sampleTimes, sampleStates,
                                          integrationMethod, n.reps)
    mapply(expect_equal, res, result)
})
