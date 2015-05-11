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
     phylo$m =m <- dim(phylo$sampleStates)[2]
     phylo$lstates <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), m)
     phylo$mstates <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), m)
     phylo$ustates <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), m)
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
     m <- nrow(births)
     nonDemeNames <- names(nonDemeDynamics)
     mm <- length(nonDemeNames)
     
     # reorder x0
     if (length(x0) != m + mm) 
        stop('initial conditons incorrect dimension', x0, m, mm) 
     if ( sum( !(c(demeNames, nonDemeNames) %in% names(x0)) )  > 0)
        stop('initial conditions vector incorrect names', names(x0), demeNames, nonDemeNames)
     y0 <- x0[c(demeNames, nonDemeNames)]
     
     # parse equations
     pbirths <- sapply(1:m, function(k) 
                    sapply(1:m, function(l)
             parse(text=births[k,l])
     ))

     if (any(is.na(migrations)))
     {
         migrations <- matrix('0', nrow=m, ncol=m)
         colnames(migrations)=rownames(migrations) <- demeNames
     }
     pmigrations <- sapply( 1:m, function(k) 
            sapply(1:m, function(l)
             parse(text=migrations[k,l])
         ))
     pdeaths <- sapply(1:m, function(k) parse(text=deaths[k]) )
     if (mm > 0) {
         pndd <- sapply(1:mm, function(k) parse(text=nonDemeDynamics[k]) )
     } else {
         pndd <- NA
     }
     
     .birth.matrix <- function( x, t) 
     {
         with(as.list(x), 
          t(matrix( sapply( 1:m^2, function(k) eval(pbirths[k]))
             , nrow=m, ncol=m
         )))
     }
     .migration.matrix <- function( x, t) 
     {
         with(as.list(x), 
          t(matrix( sapply( 1:m^2, function(k) eval(pmigrations[k]))
             , nrow=m, ncol=m
         )))
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
           sapply(1:m, function(k) eval(pdeaths[k]) )
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

# we have split the interval [min.h, max.h] into resolution equal pieces
# given a value h, return the index of the interval it falls into
get.index <- function (h, min.h, max.h, resolution) {
    pmin(1 + floor(resolution * h / (max.h - min.h)), resolution)    
}

simulate.binary.dated.tree.fgy <- function(times, births, migrations, demeSizes, sampleTimes, sampleStates, integrationMethod = 'rk4', n.reps=1, cluster=NULL)
{
    #NOTE assumes times in equal increments
    #~ TODO mstates, ustates not in returned tree 
    if (length(unique(diff(sort(times)))) != 1)
        warning('Tree simulator assumes times given in equal increments')

    # moved from binaryDatedTree
    if (!is.matrix(sampleStates)) 
        stop('sampleStates must be a matrix (not a data.frame)')

    n <- length(sampleTimes) 
    m <- ncol(sampleStates)
    maxSampleTime <- max(sampleTimes)
    if (length(names(sampleTimes)) == 0) {
        sampleNames <- as.character(1:length(sampleTimes))
        names(sampleTimes) <- sampleNames
        rownames(sampleStates) <- sampleNames
    }
    if (length(rownames(sampleStates)) == 0) 
       warning('simulate.binaryDatedTree.fgy: sampleStates should have row names')
    
    sampleHeights <- maxSampleTime - sampleTimes 
    ix <- order(sampleHeights)
    sortedSampleHeights <- sampleHeights[ix]
    sortedSampleStates <- sampleStates[ix,,drop=FALSE] 
    
    maxtime <- max(times)
    mintime <- min(times)
    maxHeight <- maxSampleTime -  mintime
    times_ix <- order(-times)
    
    FGY_RESOLUTION <- length(times)
    # reverse order (present to past): 
    F_DISCRETE <- births[times_ix]
    G_DISCRETE <- migrations[times_ix]
    Y_DISCRETE <- demeSizes[times_ix]

    # the following line accounts for any discrepancies between the maximum
    # time axis and the last sample time
    hoffset <- maxtime - maxSampleTime
    if (hoffset < 0) 
        stop( 'Time axis does not cover the last sample time' )

    # construct forcing timeseries for ode's
    heights <- sort(maxSampleTime - times)
    heights <- heights[heights <= maxHeight & heights >= 0]
    height.indices <- get.index(heights + hoffset, mintime, maxtime, FGY_RESOLUTION)

    fmat <- do.call(rbind, lapply(F_DISCRETE[height.indices], c))
    gmat <- do.call(rbind, lapply(G_DISCRETE[height.indices], c))
    ymat <- do.call(rbind, lapply(Y_DISCRETE[height.indices], c))
    fgymat <- c(pmax(cbind(fmat, gmat, ymat), 0))
    
    Q0 <- diag(m)
    .solve.Q.A.L <- function(h0, h1, A0, L0)
     { # uses C implementation
        parameters <- c(m, maxHeight, length(heights), sum(A0), fgymat)
        y0 <- c(as.vector(Q0), A0, L0)
        o <- ode(y=y0, c(h0,h1), func = "dQAL", parms = parameters, dllname = "rcolgem", initfunc = "initfunc", method=integrationMethod)
        Q1 <- t(matrix(abs(o[nrow(o),2:(1 + m^2)]), nrow=m)) #NOTE the transpose
        A1 <- o[nrow(o), (1 + m^2 + 1):(1 + m^2 +  m)]
        L1 <- o[nrow(o), ncol(o)]
        list(unname(Q1), unname(A1), unname(L1))
    }
    
    cumSortedSampleStates <- apply(sortedSampleStates, 2, cumsum)
    cumSortedNotSampledStates <- t(cumSortedSampleStates[n,] - t(cumSortedSampleStates) )

    jitter.factor <- max(1e-6, maxHeight/1e6)
    run1 <- function(repl) {
        jitter.heights <- jitter(sortedSampleHeights, factor=jitter.factor)
        # solve for A
        noisy.index <- approxfun(sort(jitter.heights), 1:n, method='constant',
                               rule=2)
        dA <- function(h, A, parms, ...)
        {
            nsy <- cumSortedNotSampledStates[ noisy.index(h), ]
            cur.y <- Y_DISCRETE[[get.index(h + hoffset, mintime, maxtime, FGY_RESOLUTION)]]
            cur.f <- F_DISCRETE[[get.index(h + hoffset, mintime, maxtime, FGY_RESOLUTION)]]
            cur.g <- G_DISCRETE[[get.index(h + hoffset, mintime, maxtime, FGY_RESOLUTION)]]

            A_Y <- (A - nsy) / cur.y
            A_Y[is.nan(A_Y)] <- 0
            csFpG <- colSums(cur.f + cur.g)
            list( setNames( as.vector(
             cur.g %*% A_Y - csFpG * A_Y + (cur.f %*% A_Y) * pmax(1-A_Y, 0)
             ), names(A)
            ))
        }
        h0 <- 0
        sampled.at.h <- function(h) which(sortedSampleHeights == h)

        haxis <- seq(0, maxHeight, length.out=FGY_RESOLUTION)
        odA <-  ode(y = colSums(sortedSampleStates), times = haxis, func=dA,
                   parms = NA, method = integrationMethod)

        haxis <- odA[,1]
        AplusNotSampled <- odA[,2:(m+1),drop=FALSE]
        Amono <- rowSums( AplusNotSampled )
        Amono[is.na(Amono)] <- min(Amono, na.rm=TRUE)
        Amono <- Amono - min(Amono)
        Amono <- (max(Amono)  - Amono) / max(Amono)
        nodeHeights <- sort( approx( Amono, haxis, xout=runif(n-1, 0, 1) )$y ) # careful of impossible node heights
        eventTimes <- c(unique(sampleHeights), nodeHeights)
        isSampleEvent <- c(rep(TRUE, length(unique(sampleHeights))), rep(FALSE, length(nodeHeights)))
        ix <- order(eventTimes)
        eventTimes <- eventTimes[ix]
        isSampleEvent <- isSampleEvent[ix]

        get.A <- function(h) {
            i <- get.index(h, 0, maxHeight, FGY_RESOLUTION)
            AplusNotSampled[i,] - cumSortedNotSampledStates[ noisy.index(h), ]
        }

        S <- 1
        L <- 0

        # initialize variables; tips in order of sortedSampleHeights
        Nnode <- n-1
        edge.length <- rep(-1, Nnode + n-1) # should not have root edge
        edge <- matrix(-1, (Nnode + n-1), 2)

        if (length(names(sortedSampleHeights))==0)
            tip.label <- as.character(1:n)
        else
            tip.label <- names(sortedSampleHeights)

        heights       <- rep(0, (Nnode + n) )
        parentheights <- rep(-1, (Nnode + n) )
        heights[1:n]  <- sortedSampleHeights
        inEdgeMap     <- rep(-1, Nnode + n)
        outEdgeMap    <- matrix(-1, (Nnode + n), 2)
        parent        <- 1:(Nnode + n)
        daughters     <- matrix(-1, (Nnode + n), 2)
        lstates       <- matrix(-1, (Nnode + n), m)
        mstates       <- matrix(-1, (Nnode + n), m)
        ustates       <- matrix(-1, (Nnode + n), m)
        ssm           <- matrix( 0, nrow=n, ncol=m)
        lstates[1:n,] <-  sortedSampleStates
        mstates[1:n,] <- lstates[1:n,]

        isExtant <- rep(FALSE, Nnode+n)
        isExtant[sampled.at.h(h0)] <- TRUE
        extantLines <- which(isExtant)

        if (length(extantLines) > 1){
            A0 <- colSums(as.matrix(sortedSampleStates[extantLines,], nrow=length(extantLines)) )
        } else{
            A0 <- sortedSampleStates[extantLines,]
        }
        lineageCounter <- n+1

        event.F <- F_DISCRETE[get.index(tail(eventTimes, -1) + hoffset, mintime, maxtime, FGY_RESOLUTION)]
        event.G <- G_DISCRETE[get.index(tail(eventTimes, -1) + hoffset, mintime, maxtime, FGY_RESOLUTION)]
        event.Y <- Y_DISCRETE[get.index(tail(eventTimes, -1) + hoffset, mintime, maxtime, FGY_RESOLUTION)]

        for (ih in 1:(length(eventTimes)-1)) {
            h0 <- eventTimes[ih]
            h1 <- eventTimes[ih+1]
            nExtant <- sum(isExtant)

            #get A0, process new samples, calculate state of new lines
            A0 <- get.A(h0) 
            out <- .solve.Q.A.L(h0, h1, A0,  L)
            Q <- out[[1]]
            A <- out[[2]]
            L <- out[[3]]

            # clean output
            if (is.nan(L)) {L <- Inf}
            if (sum(is.nan(Q)) > 0) Q <- diag(length(A))
            if (sum(is.nan(A)) > 0) A <- A0

            #update mstates
            if ( nExtant > 1)
            {
                mstates[isExtant,] <- t( t(Q) %*% t(mstates[isExtant,])  )
                mstates[isExtant,] <- abs(mstates[isExtant,]) / rowSums(abs(mstates[isExtant,,drop=FALSE]))
                #recalculate A
                A <- colSums(mstates[isExtant,,drop=FALSE])
            }
            else{
                mstates[isExtant,] <- t( t(Q) %*% mstates[isExtant,] )
                mstates[isExtant,] <- abs(mstates[isExtant,]) / sum(abs(mstates[isExtant,]))
                #recalculate A
                A <- mstates[isExtant,,drop=FALSE]
            }

            if (isSampleEvent[ih+1])
            {
                 sat_h1 <- sampled.at.h(h1)
                 isExtant[sat_h1] <- TRUE
                 heights[sat_h1] <- h1
            } else { #coalecent event
                .F <- event.F[[ih]]
                .G <- event.G[[ih]]
                .Y <- event.Y[[ih]]

                if (nExtant > 1 && 0 %in% .Y) 
                    # last two lineages have not coalesced by simulation time zero
                    # reject this tree and return NA to prevent divide-by-zero error
                    return(NA)

                a <- A / .Y
                extantLines <- which(isExtant)

                .lambdamat <- (t(t(a)) %*% a) * .F

                # determine transmission type from deme A to deme B
                # m is number of demes
                # m^2 is number of transmission types
                kl <- sample.int( m^2, size=1, prob=as.vector(.lambdamat) )
                k <- 1 + ((kl-1) %% m)#row
                l <- 1 + floor( (kl-1) / m ) #column

                # mstates stores probabilities of deme membership (columns)
                #  for all lineages (rows)
                probstates <- as.matrix(mstates[extantLines,], nrow=length(extantLines))

                # which extant lineages are source and recipient?
                u_i <- sample.int(  nExtant, size=1, prob = probstates[,k])
                probstates[u_i,] <- 0 # cant transmit to itself
                u <- extantLines[u_i]
                v <- sample(  extantLines, size=1, prob = probstates[,l])

                ustates[u,] <- mstates[u,]
                ustates[v,] <- mstates[v,]

                a_u <- pmin(1, mstates[u,] / .Y )
                a_v <- pmin(1, mstates[v,] / .Y )
                lambda_uv <- ((a_u) %*% t( a_v )) * .F + ((a_v) %*% t( a_u )) * .F

                palpha <- rowSums(lambda_uv) / sum(lambda_uv)
                alpha <- lineageCounter
                lineageCounter <- lineageCounter + 1
                isExtant[alpha] <- TRUE
                isExtant[u] <- isExtant[v] <- FALSE
                lstates[alpha,] <- mstates[alpha,] <- palpha
                heights[alpha] <- h1

                uv <- c(u, v)
                inEdgeMap[uv] <- alpha #lineageCounter
                outEdgeMap[alpha,] <- uv
                parent[uv] <- alpha
                parentheights[uv] <- h1
                daughters[alpha,] <- uv
                edge[u,] <- c(alpha, u)
                edge.length[u] <- h1 - heights[u]
                edge[v,] <- c(alpha,v)
                edge.length[v] <- h1 - heights[v]
            } # end else
        } # end for loop

        self <- list(edge=edge, edge.length=edge.length, Nnode=Nnode,
                     tip.label=tip.label, heights=heights,
                     parentheights=parentheights, parent=parent,
                     daughters=daughters, lstates=lstates, mstates=mstates,
                     ustates=ustates, m=m, sampleTimes = sampleTimes,
                     sampleStates= sampleStates, maxSampleTime=maxSampleTime,
                     inEdgeMap = inEdgeMap, outEdgeMap=outEdgeMap)

        class(self) <- c("binaryDatedTree", "phylo")
        # reorder edges for compatibility with ape::phylo functions
        # (ideally ape would not care about the edge order, but actually
        # most functions assume a certain order)
        sampleTimes2 <- sampleTimes[names(sortedSampleHeights)]
        sampleStates2 <- lstates[1:n,,drop=FALSE]
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
    na.omit(result)
}
