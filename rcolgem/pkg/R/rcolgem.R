#' this file contains all R functions of the rcolgem package
#' @import ape
#' @import deSolve

#~ DONE  finite size correction for pair (i,j) of lineages at internal node (see written notes): 
#~     if i transmits and is type k: 
#~          pjk -> pjk * (Yk-1)/Yk
#~          pjl (l\neq k) -> pjl * Yk/(Yk-pjk)
#~ TODO option to correct for direct ancestor sampling if doing serial samples and there is a 0-branch length
#~      add these terms to the likelihood:
#~          if 0 bl at s_i: \sum_k pik Ak / Yk
#~          else: \sum_k pik (1-Ak/Yk)
#~ DONE validate input & raise warnings
#~      check sampleTimes compatible with edge.length; FGY functions defined over length of tree;
#~ DONE adapt coalescent simulator to heterochronous sample
#~ TODO add option to switch to semi-structured coalescent if difference in line states below threshold

#PRELIMINARIES 
.calculate.heights <- function(phylo){
     phylo$maxSampleTime   <- max(phylo$sampleTimes)
     heights <- rep(0, (phylo$Nnode + length(phylo$tip.label)) )
     heights[1:length(phylo$sampleTimes)] <- phylo$maxSampleTime - phylo$sampleTimes
     curgen <- 1:length(phylo$sampleTimes)
     while( length(curgen) > 0) { 
         nextgen <- c()
         icurgenedges <- which(  phylo$edge[,2] %in% curgen  )
         for (i in icurgenedges){
             u<- phylo$edge[i,1]
             v<- phylo$edge[i,2]
             # inspect tree
             if ( heights[u] > 0 & abs(heights[u] - (phylo$edge.length[i] + heights[v]))/heights[u] > 1e-2 )
             { #browser()
               stop( 'Tree is poorly formed. Branch lengths incompatible with sample times.')
               }
             heights[u] <- phylo$edge.length[i] + heights[v]
             nextgen <- c(nextgen, u)
         }
         curgen <- unique(nextgen)
     }
     phylo$heights <- heights
     phylo$maxHeight <- max(heights)
     return(phylo)
}

.calculate.edgemap <- function(phylo){
     inEdgeMap <- rep(-1, length((phylo$Nnode + length(phylo$tip.label))))
     outEdgeMap <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), 2)
     parent <- 1:(phylo$Nnode + length(phylo$tip.label)) 
     daughters <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), 2)
     for (u in 1:(phylo$Nnode + length(phylo$tip.label))){
         #if (u!=length(phylo$tip.label)+1){ #if u not root
         tryCatch({ inEdgeMap[u] <- which( phylo$edge[,2]==u ) }, error = function(e) {inEdgeMap[u] <- u} )
         #} else{ 
         #     inEdgeMap[u] <- u
         #}
         if (u > length(phylo$tip.label)){
             outEdgeMap[u,] <- which( phylo$edge[,1]==u ) 
             daughters[u,] <- phylo$edge[outEdgeMap[u,],2]
         } else{ 
             outEdgeMap[u,] <- c(NA, NA)
             daughters[u,] <- c(NA, NA)
         }
         parent[u] <- phylo$edge[inEdgeMap[u],1]
     }
     phylo$inEdgeMap <- inEdgeMap
     phylo$outEdgeMap <- outEdgeMap
     phylo$parent = phylo$parents <- parent
     phylo$daughter = phylo$daughters <- daughters
     phylo$parentheight = phylo$parentheights <- phylo$heights[parent[1:(phylo$Nnode + length(phylo$tip.label))]]
     return(phylo)
}

.initialize.states <- function(phylo)
{
     phylo$m <- n.demes <- dim(phylo$sampleStates)[2]
     phylo$lstates <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), n.demes)
     phylo$mstates <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), n.demes)
     phylo$ustates <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), n.demes)
     phylo$lstates[1:nrow(phylo$sampleStates),] <- phylo$sampleStates
     phylo$mstates[1:nrow(phylo$sampleStates),] <- phylo$sampleStates
     phylo$ustates[1:nrow(phylo$sampleStates),] <- phylo$sampleStates
     return(phylo)
}

#' Create binary dated tree
#' binaryDatedTree class, includes heights for each node and other helper variables
#' @export
#~ binaryDatedTree <- function( x, ...) UseMethod("binaryDatedTree")
#~ binaryDatedTree.default <- function( phylo, sampleTimes, sampleStates){
binaryDatedTree <- function( phylo, sampleTimes, sampleStates=NULL, sampleStatesAnnotations=NULL){
     if (phylo$Nnode != length(phylo$tip.label) - 1 ) { stop('Object class phylo is not a binary tree.') }
     if (is.null(names(sampleTimes))) stop('sampleTimes vector must have names of tip labels')
     if (is.null(sampleStates) & !is.null(sampleStatesAnnotations) ) sampleStates <- .infer.sample.states.from.annotation(phylo, sampleStatesAnnotations)
     if (is.null(sampleStates) & is.null(sampleStatesAnnotations)) { sampleStates <- t(t( rep(1, length(phylo$tip.label)))) ; rownames( sampleStates) <- phylo$tip.label }
     if (is.null(rownames(sampleStates))) stop('sampleStates matrix must have row names of tip labels')
     
     phylo$sampleTimes <- sampleTimes[phylo$tip.label]
     phylo$sampleStates <- sampleStates[phylo$tip.label, ]
     if (is.vector(phylo$sampleStates)) phylo$sampleStates <- t(t( phylo$sampleStates))
     phylo <- .calculate.heights(phylo)
     phylo <- .calculate.edgemap(phylo)
     phylo <- .initialize.states(phylo)
     phylo$coalescentRates <- rep(0, (phylo$Nnode + length(phylo$tip.label)))
     phylo$coalescentSurvivalProbability <- rep(0, (phylo$Nnode + length(phylo$tip.label)))
     phylo$logCoalescentSurvivalProbability <- rep(-Inf, (phylo$Nnode + length(phylo$tip.label)))
     
     # correct any funniness due to very small or negative branch lengths
     inodes <- (length(phylo$tip.label)+1):length(phylo$heights)
     while ( any( phylo$heights[inodes] < phylo$heights[ phylo$daughters[inodes,1] ] ) |  any( phylo$heights[inodes] < phylo$heights[ phylo$daughters[inodes,2] ] )  )
     {
         phylo$heights[inodes] <- sapply( inodes, function(alpha) max(phylo$heights[alpha], phylo$heights[phylo$daughters[alpha,] ] ) ) 
     } 
     #  edge.lengths may have changed- repair them: 
     for (i in 1:nrow(phylo$edge))
     {
         phylo$edge.length[i] <- phylo$heights[phylo$edge[i,1]] - phylo$heights[phylo$edge[i,2]]
     }
     phylo$n <- length(phylo$tip.label)
     class(phylo) <- c("binaryDatedTree", "phylo")
     return(phylo)
}


make.fgy <- function(t0, t1, births, deaths, nonDemeDynamics,  x0,  migrations=NA,  parms=NA, fgyResolution = 2000, integrationMethod = 'rk4')
{
#~      generates a discrete numeric representation of the demographic process given ODE as type str
     demeNames <- rownames(births)
     n.demes <- nrow(births)
     nonDemeNames <- names(nonDemeDynamics)
     mm <- length(nonDemeNames)
     
     # reorder x0
     if (length(x0) != n.demes + mm) 
        stop('initial conditons incorrect dimension', x0, n.demes, mm) 
     if ( sum( !(c(demeNames, nonDemeNames) %in% names(x0)) )  > 0)
        stop('initial conditions vector incorrect names', names(x0), demeNames, nonDemeNames)
     y0 <- x0[c(demeNames, nonDemeNames)]
     
     # parse equations
     pbirths <- sapply(1:n.demes, function(k) 
                    sapply(1:n.demes, function(l)
             parse(text=births[k,l])
     ))

     if (any(is.na(migrations)))
     {
         migrations <- matrix('0', nrow=n.demes, ncol=n.demes)
         colnames(migrations)=rownames(migrations) <- demeNames
     }
     pmigrations <- sapply( 1:n.demes, function(k) 
            sapply(1:n.demes, function(l)
             parse(text=migrations[k,l])
         ))
     pdeaths <- sapply(1:n.demes, function(k) parse(text=deaths[k]) )
     if (mm > 0) {
         pndd <- sapply(1:mm, function(k) parse(text=nonDemeDynamics[k]) )
     } else {
         pndd <- NA
     }
     
     .birth.matrix <- function( x, t) 
     {
         with(as.list(x), 
          t(matrix( sapply( 1:n.demes^2, function(k) eval(pbirths[k]))
             , nrow=n.demes, ncol=n.demes
         )))
     }
     .migration.matrix <- function( x, t) 
     {
         with(as.list(x), 
          t(matrix( sapply( 1:n.demes^2, function(k) eval(pmigrations[k]))
             , nrow=n.demes, ncol=n.demes)))
     }
     tBirths <- function(x, t)
     {
         colSums( .birth.matrix(x,t) )
     }
     tMigrationsIn <- function(x,t)
     {
         colSums( .migration.matrix(x, t) )
     }
     tMigrationsOut <- function(x,t)
     {
         rowSums( .migration.matrix(x, t))
     }
     tDeaths <- function(x, t) 
     {
         with(as.list(x, t), 
           sapply(1:n.demes, function(k) eval(pdeaths[k]) )
         ) 
     }
     
     dNonDeme <- function(x, t) 
     {
         with(as.list(x, t), 
           sapply(1:mm, function(k) eval(pndd[k]) )  
         )
     }
     dx <- function(t, y, parms, ...) 
     {
         dxdeme <- setNames( tBirths(y, t) + tMigrationsIn(y, t) - tMigrationsOut(y,t) - tDeaths(y,t), demeNames)
         if (mm > 0)
         {
             dxnondeme <- setNames( dNonDeme(y, t), nonDemeNames )
         }else{
             dxnondeme <- NULL
         }
         
         list( c(dxdeme, dxnondeme) )
     }
     
     times1 <- seq(t0, t1, length.out=fgyResolution)
     dh <- times1[2] - times1[1]
     if (t0 - dh <= t0) { times0 <- c() }
     else{ times0 <-  seq(t0, t0 - dh, by=dh) }
     times <- unique( c( times0, times1) )
     ox <- ode(y=y0, times, func=dx, parms, method=integrationMethod)
     # note does not include first value, which is t0; 2nd value corresponds to root of tree
     Ys <- lapply( nrow(ox):(1+length(times0)), function(i) ox[i, demeNames] )
     Fs <- lapply( nrow(ox):(1+length(times0)), function(i) .birth.matrix(ox[i,], ox[i,1])  ) # dh *
     Gs <- lapply( nrow(ox):(1+length(times0)), function(i) .migration.matrix(ox[i,], ox[i,1])  ) #dh *
     
     list( rev(times), Fs, Gs, Ys , ox )
}

# are all the elements of a vector numerically equal?
# http://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector
vector.all.equal <- function (x) 
{
    diff(range(x)) < .Machine$double.eps ^ 0.5
}

# we have split the interval [min.h, max.h] into resolution equal pieces
# given a value h, return the index of the interval it falls into
get.index <- function (h, min.h, max.h, resolution) {
    pmin(1 + floor(resolution * h / (max.h - min.h)), resolution)    
}

simulate.binary.dated.tree.fgy <- function(times, births, migrations, demeSizes, sampleTimes, sampleStates, integrationMethod = 'rk4', n.reps=1, cluster=NULL)
{
    #NOTE assumes times in equal increments
    #~ TODO mstates, ustates not in returned tree 
    if (vector.all.equal(diff(sort(times))) != 1)
        warning('Tree simulator assumes times given in equal increments')

    # moved from binaryDatedTree
    if (!is.matrix(sampleStates)) 
        stop('sampleStates must be a matrix (not a data.frame)')

    # assign taxa names if they don't have them already
    if (length(names(sampleTimes)) == 0) {
        sampleNames <- as.character(1:length(sampleTimes))
        names(sampleTimes) <- rownames(sampleStates) <- sampleNames
    }

    # make sure sample times and sample states match up
    if (length(rownames(sampleStates)) == 0) 
        warning('simulate.binaryDatedTree.fgy: sampleStates should have row names')
    else if (!all(sort(rownames(sampleStates)) == sort(names(sampleTimes))))
        warning('names of sampleStates and sampleTimes are different')
    else
        sampleStates <- sampleStates[match(names(sampleTimes), rownames(sampleStates)),,drop=FALSE]

    n.taxa <- length(sampleTimes) 
    n.demes <- ncol(sampleStates)
    maxSampleTime <- max(sampleTimes)

    # order sample heights by time
    sampleHeights <- maxSampleTime - sampleTimes 
    ix <- order(sampleHeights)
    sampleHeights <- sampleHeights[ix]
    sortedSampleStates <- sampleStates[ix,,drop=FALSE] 
    
    maxHeight <- maxSampleTime -  min(times)
    times_ix <- order(-times)
    
    FGY_RESOLUTION <- length(times)
    # reverse order (present to past): 
    F_DISCRETE <- births[times_ix]
    G_DISCRETE <- migrations[times_ix]
    Y_DISCRETE <- demeSizes[times_ix]

    if (max(times) < maxSampleTime) 
        stop( 'Time axis does not cover the last sample time' )

    # construct forcing timeseries for ode's
    heights <- sort(maxSampleTime - times)
    heights <- heights[heights >= 0]

    fmat <- do.call(rbind, lapply(F_DISCRETE, c))
    gmat <- do.call(rbind, lapply(G_DISCRETE, c))
    ymat <- do.call(rbind, lapply(Y_DISCRETE, c))
    fgymat <- c(pmax(cbind(fmat, gmat, ymat), 0))
    
    Q0 <- c(diag(n.demes))
    .solve.Q.A.L <- function(h0, h1, A0, L0)
     { # uses C implementation
        parameters <- c(n.demes, maxHeight, length(heights), sum(A0), fgymat)
        y0 <- c(Q0, A0, L0)
        o <- ode(y=y0, c(h0,h1), func = "dQAL", parms = parameters, dllname =
                 "rcolgem", initfunc = "initfunc", method=integrationMethod)
        Q1 <- t(matrix(abs(o[nrow(o), 2:(1 + n.demes^2)]), nrow=n.demes)) #NOTE the transpose
        L1 <- o[nrow(o), ncol(o)]

        # clean output
        if (is.nan(L1)) L1 <- Inf
        if (sum(is.nan(Q1)) > 0) Q1 <- diag(n.demes)

        list(unname(Q1), unname(L1))
    }

    test.qal <- function (h0, h1, A0, L0, events, times)
    {
        parameters <- c(n.demes, maxHeight, length(heights), sum(A0), fgymat)
        y0 <- c(Q0, A0, L0, sum(A0))
        o <- ode(y = y0, times = times, func = "dQAL_RM", parms = parameters, 
                 dllname = "rcolgem", initfunc = "initfunc",
                 method = integrationMethod,
                 events = list(data=events))

        # extract Q matrices
        Q1 <- o[,2:(n.demes^2+1),drop=FALSE]
        Q1 <- lapply(1:nrow(Q1), function (i) matrix(Q1[i,], byrow=TRUE, nrow=n.demes))
        return(Q1)

        #Q1 <- lapply(o[, 3:(n.demes^2+2), drop=FALSE], matrix, nrow=n.demes)
        Q1 <- t(matrix(abs(o[nrow(o), 2:(1 + n.demes^2)]), nrow=n.demes)) #NOTE the transpose
        L1 <- o[nrow(o), ncol(o)-1]

        # clean output
        if (is.nan(L1)) {L1 <- Inf}
        if (sum(is.nan(Q1)) > 0) Q1 <- diag(n.demes)

        list(unname(Q1), unname(L1))
    }
    
    # cumSortedSampleStates[s, d]: at the time of collection of sample s (going
    # backwards in time), how many samples of deme d have been collected?
    cumSortedSampleStates <- apply(sortedSampleStates, 2, cumsum)

    # cumSortedNotSampleStates[s, d]: at the time of collection of sample s (going
    # backwards in time), how many samples of deme d remain to be collected?
    cumSortedNotSampledStates <- t(cumSortedSampleStates[n.taxa,] - t(cumSortedSampleStates))

    haxis <- seq(0, maxHeight, length.out=FGY_RESOLUTION)
    jitter.factor <- max(1e-6, maxHeight/1e6)

    run1 <- function(repl) {
        jitter.heights <- jitter(sampleHeights, factor=jitter.factor)
        # solve for A
        noisy.index <- approxfun(sort(jitter.heights), 1:n.taxa, method='constant',
                                 rule=2)
        dA <- function(h, A, parms, ...)
        {
            nsy <- cumSortedNotSampledStates[ noisy.index(h), ]
            cur.idx <- get.index(h, 0, maxHeight, FGY_RESOLUTION)
            cur.y <- Y_DISCRETE[[cur.idx]]
            cur.f <- F_DISCRETE[[cur.idx]]
            cur.g <- G_DISCRETE[[cur.idx]]

            A_Y <- (A - nsy) / cur.y
            A_Y[is.nan(A_Y)] <- 0
            csFpG <- colSums(cur.f + cur.g)
            list( setNames( c(
             cur.g %*% A_Y - csFpG * A_Y + (cur.f %*% A_Y) * pmax(1-A_Y, 0)
             ), names(A)
            ))
        }
        odA <-  ode(y = colSums(sortedSampleStates), times = haxis, func=dA,
                    parms = NA, method = integrationMethod)
        AplusNotSampled <- odA[,2:(n.demes+1),drop=FALSE]

        # Amono is the total number of lineages (of all demes) present at each time step
        Amono <- rowSums( AplusNotSampled )
        Amono[is.na(Amono)] <- min(Amono, na.rm=TRUE)

        # transform Amono to [0, 1], reversing the magnitudes
        Amono <- (max(Amono) - Amono) / (diff(range(Amono)))

        # choose node heights by interpolating between the discrete A values
        nodeHeights <- sort( approx( Amono, haxis, xout=runif(n.taxa-1, 0, 1) )$y ) # careful of impossible node heights

        # combine sample and coalescence times
        eventTimes <- c(unique(sampleHeights), nodeHeights)
        isSampleEvent <- c(rep(TRUE, length(unique(sampleHeights))), rep(FALSE, length(nodeHeights)))
        ix <- order(eventTimes)
        eventTimes <- eventTimes[ix]
        isSampleEvent <- isSampleEvent[ix]

        # initialize variables; tips in order of sampleHeights (which is sorted)
        Nnode <- n.taxa - 1
        edge.length <- rep(-1, Nnode + n.taxa - 1) # should not have root edge
        edge <- matrix(-1, (Nnode + n.taxa - 1), 2)

        if (length(names(sampleHeights))==0)
            tip.label <- as.character(1:n.taxa)
        else
            tip.label <- names(sampleHeights)

        heights       <- c(sampleHeights, rep(0, Nnode))
        parentheights <- rep(-1, (Nnode + n.taxa) )
        inEdgeMap     <- rep(-1, Nnode + n.taxa)
        outEdgeMap    <- matrix(-1, (Nnode + n.taxa), 2)
        parent        <- 1:(Nnode + n.taxa)
        daughters     <- matrix(-1, (Nnode + n.taxa), 2)
        lstates       <- matrix(-1, (Nnode + n.taxa), n.demes)
        mstates       <- matrix(-1, (Nnode + n.taxa), n.demes)
        ustates       <- matrix(-1, (Nnode + n.taxa), n.demes)
        lstates[1:n.taxa,] <- sortedSampleStates
        mstates[1:n.taxa,] <- lstates[1:n.taxa,]

        isExtant <- rep(FALSE, Nnode+n.taxa)
        extantLines <- which(sampleHeights == 0)
        isExtant[extantLines] <- TRUE

        if (length(extantLines) > 1){
            A0 <- colSums(sortedSampleStates[extantLines,,drop=FALSE])
        } else{
            A0 <- sortedSampleStates[extantLines,]
        }

        alpha <- n.taxa+1

        h.index <- get.index(eventTimes, 0, maxHeight, FGY_RESOLUTION)
        event.F <- F_DISCRETE[tail(h.index, -1)]
        event.Y <- Y_DISCRETE[tail(h.index, -1)]

        event.A0 <- AplusNotSampled[head(h.index, -1),,drop=FALSE] -
                    cumSortedNotSampledStates[ noisy.index(head(eventTimes, -1)),,drop=FALSE]

        n.events <- nrow(event.A0)
        event.data.Q0 <- data.frame(var = rep(1:n.demes^2, each=n.events),
                                    time = head(eventTimes, -1),
                                    value = rep(c(diag(n.demes)), each=n.events),
                                    method = "rep")
        event.data.A0 <- data.frame(var = rep((n.demes^2+1):(n.demes^2+n.demes), each=n.events),
                                    time = head(eventTimes, -1),
                                    value = c(event.A0),
                                    method = "rep")
        event.data.A0sum <- data.frame(var = rep(n.demes^2 + n.demes + 2, n.events),
                                       time = head(eventTimes, -1),
                                       value = rowSums(event.A0),
                                       method = "rep")
        event.data <- do.call(rbind, list(event.data.Q0, event.data.A0, event.data.A0sum))
        event.Q <- test.qal(head(eventTimes, 1), tail(eventTimes, 1),
                            event.A0[1,], 0, event.data, eventTimes)

        for (ih in 1:(length(eventTimes)-1)) {
            h0 <- eventTimes[ih]
            h1 <- eventTimes[ih+1]

            #get A0, process new samples, calculate state of new lines
            Q <- event.Q[[ih+1]]

            # update mstates
            nExtant <- sum(isExtant)
            mstates[isExtant,] <- mstates[isExtant,,drop=FALSE] %*% Q
            mstates[isExtant,] <- abs(mstates[isExtant,,drop=FALSE]) / rowSums(abs(mstates[isExtant,,drop=FALSE]))

            # recalculate A
            A <- colSums(mstates[isExtant,,drop=FALSE])

            if (isSampleEvent[ih+1]) {
                isExtant[which(sampleHeights == h1)] <- TRUE
            } else { #coalecent event
                F <- event.F[[ih]]
                Y <- event.Y[[ih]]

                if (nExtant > 1 && 0 %in% Y) 
                    # last two lineages have not coalesced by simulation time zero
                    # reject this tree and return NA to prevent divide-by-zero error
                    return(NA)

                # a is the proportion of extant lineages in each deme
                a <- A / Y 
                extantLines <- which(isExtant)

                .lambdamat <- (t(t(a)) %*% a) * F

                # determine transmission type from deme A to deme B
                # n.demes is number of demes
                # n.demes^2 is number of transmission types
                kl <- sample.int(n.demes^2, size=1, prob=c(.lambdamat) )
                k <- 1 + ((kl-1) %% n.demes)#row
                l <- 1 + floor( (kl-1) / n.demes ) #column

                # mstates stores probabilities of deme membership (columns)
                #  for all lineages (rows)
                probstates <- mstates[extantLines,,drop=FALSE]

                # which extant lineages are source and recipient?
                u_i <- sample.int(nExtant, size=1, prob = probstates[,k])
                probstates[u_i,] <- 0 # can't transmit to itself
                u <- extantLines[u_i]
                v <- sample(extantLines, size=1, prob = probstates[,l])

                uv <- c(u, v)
                ustates[uv,] <- mstates[uv,]

                a_u <- pmin(1, mstates[u,] / Y )
                a_v <- pmin(1, mstates[v,] / Y )
                lambda_uv <- (a_u %*% t(a_v) + a_v %*% t(a_u)) * F

                palpha <- rowSums(lambda_uv) / sum(lambda_uv)
                isExtant[alpha] <- TRUE
                isExtant[u] <- isExtant[v] <- FALSE
                lstates[alpha,] <- mstates[alpha,] <- palpha
                heights[alpha] <- h1

                inEdgeMap[uv] <- alpha
                outEdgeMap[alpha,] <- uv
                parent[uv] <- alpha
                parentheights[uv] <- h1
                daughters[alpha,] <- uv
                edge[u,] <- c(alpha, u)
                edge[v,] <- c(alpha, v)
                edge.length[u] <- h1 - heights[u]
                edge.length[v] <- h1 - heights[v]

                alpha <- alpha + 1
            } # end else
        } # end for loop

        self <- list(edge=edge, edge.length=edge.length, Nnode=Nnode,
                     tip.label=tip.label, heights=heights,
                     parentheights=parentheights, parent=parent,
                     daughters=daughters, lstates=lstates, mstates=mstates,
                     ustates=ustates, n.demes=n.demes, sampleTimes = sampleTimes,
                     sampleStates=sampleStates, maxSampleTime=maxSampleTime,
                     inEdgeMap=inEdgeMap, outEdgeMap=outEdgeMap)

        class(self) <- c("binaryDatedTree", "phylo")
        # reorder edges for compatibility with ape::phylo functions
        # (ideally ape would not care about the edge order, but actually
        # most functions assume a certain order)
        sampleTimes2 <- sampleTimes[names(sampleHeights)]
        sampleStates2 <- lstates[1:n.taxa,,drop=FALSE]
        rownames(sampleStates2) <- tip.label
        phylo <- read.tree(text=write.tree(self) )

        sampleTimes2 <- sampleTimes2[phylo$tip.label]
        sampleStates2 <- sampleStates2[phylo$tip.label,,drop=FALSE]
        binaryDatedTree(phylo, sampleTimes2, sampleStates = sampleStates2)
    } # end run1()

    if (any(is.null(cluster))) {
        # single-threaded mode
        result <- lapply(1:n.reps, run1)
    } else {
        result <- parLapply(cluster, 1:n.reps, run1)
    }

    # exclude NA values
    result[!is.na(result)]
}
