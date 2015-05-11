context("Tree")

test_that("A tree can be simulated from an SI model", {
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

test_that("A tree can be simulated with an SI2 model", {
    si2.fgy <- setup.si2.fgy()
    sampleStates <- matrix(1, nrow=n.tips, ncol=1, dimnames=list(1:n.tips, "I"))
    tree <- simulate.binary.dated.tree.fgy(si2.fgy[[1]], si2.fgy[[2]], si2.fgy[[3]], si2.fgy[[4]],
                                           sampleTimes, sampleStates, integrationMethod)[[1]]
    treestr <- "((3:55.55862338,(2:2.597006421,4:2.597006421):52.96161696):34.21260406,(1:55.91584579,5:55.91584579):33.85538165);"
    expected.tree <- read.tree(text=treestr)
    expect_equal(tree, expected.tree)
})

test_that("A tree can be simulated with a DiffRisk model", {
    diffrisk.fgy <- setup.diffrisk.fgy()
    sampleStates <- matrix(rep(c(1, 0, 0, 1), each=5), ncol=2, dimnames=list(1:10, c("I1", "I2")))
    sampleTimes <- rep(100, 10)
    tree <- simulate.binary.dated.tree.fgy(diffrisk.fgy[[1]], diffrisk.fgy[[2]], diffrisk.fgy[[3]], diffrisk.fgy[[4]],
                                           sampleTimes, sampleStates, integrationMethod)[[1]]
    treestr <- "(((6:10.00063456,7:10.00063456):51.68538645,(1:54.8638819,(2:52.12075637,(4:32.5721444,(3:13.69108906,5:13.69108906):18.88105534):19.54861197):2.743125534):6.822139108):25.16850896,(10:72.43834329,(8:41.77070662,9:41.77070662):30.66763667):14.41618668);"
    expected.tree <- read.tree(text=treestr)
    expect_equal(tree, expected.tree)
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

test_that("Replicate samples are consistent", {
    si.fgy <- setup.si.fgy()
    sampleStates <- matrix(1, nrow=n.tips, ncol=1, dimnames=list(1:n.tips, "I"))
    trees <- simulate.binary.dated.tree.fgy(si.fgy[[1]], si.fgy[[2]], si.fgy[[3]], si.fgy[[4]],
                                            sampleTimes, sampleStates, integrationMethod, 
                                            n.reps=2)

    treestrs <- c(
        "((1:13.74595713,3:13.74595713):57.55451524,(4:50.67885423,(2:39.6305724,5:39.6305724):11.04828183):20.62161814);",
        "(5:92.69174089,(2:83.0433153,(3:75.53197339,(1:52.45128594,4:52.45128594):23.08068745):7.511341914):9.648425583);"
    )
    expected.trees <- mapply(read.tree, text=treestrs, SIMPLIFY=FALSE)
    mapply(expect_equal, trees, expected.trees)
})

test_that("Heterochronous samples are consistent", {
    si.fgy <- setup.si.fgy()
    sampleStates <- matrix(1, nrow=n.tips, ncol=1, dimnames=list(1:n.tips, "I"))
    sampleTimes <- c(rep(si.t1, n.tips-2), round(3*si.t1/4), round(si.t1/2))
    trees <- simulate.binary.dated.tree.fgy(si.fgy[[1]], si.fgy[[2]], si.fgy[[3]], si.fgy[[4]],
                                            sampleTimes, sampleStates, integrationMethod, 
                                            n.reps=2)

    treestrs <- c(
        "(5:36.20938295,(4:37.52339929,(1:59.38273654,(2:55.23945964,3:55.23945964):4.143276902):3.140662751):23.68598366);",
        "(2:64.15499076,((1:51.47312362,4:26.47312362):11.65615768,(3:59.55398855,5:9.553988554):3.575292749):1.025709462);"
    )
    expected.trees <- mapply(read.tree, text=treestrs, SIMPLIFY=FALSE)
    mapply(expect_equal, trees, expected.trees)
})

