#' this file contains all R functions of the rcolgem package
#' @import ape
#' @import deSolve

#~ DONE  finite size correction for pair (i,j) of lineages at internal node (see written notes): 
#~	if i transmits and is type k: 
#~ 		pjk -> pjk * (Yk-1)/Yk
#~ 		pjl (l\neq k) -> pjl * Yk/(Yk-pjk)
#~ TODO option to correct for direct ancestor sampling if doing serial samples and there is a 0-branch length
#~ 	add these terms to the likelihood:
#~ 		if 0 bl at s_i: \sum_k pik Ak / Yk
#~ 		else: \sum_k pik (1-Ak/Yk)
#~ DONE validate input & raise warnings
#~ 	check sampleTimes compatible with edge.length; FGY functions defined over length of tree;
#~ DONE adapt coalescent simulator to heterochronous sample
#~ TODO add option to switch to semi-structured coalescent if difference in line states below threshold



#' @export
SIMULATIONTIMERESOLUTION<- 1e+04

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
		#	inEdgeMap[u] <- u
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

#' If the taxon labels end in _<label>, treat <label> as the discrete 
#' state of the taxon, eg location of sampling
.infer.sample.states.from.annotation <- function(phylo, sampleStatesAnnotations)
{
	annotations <- regmatches( phylo$tip.label, regexpr( '_[.,-]*[[:alnum:]]+$', phylo$tip.label) )
	annotations <- unname( sapply( annotations, function(a) substr(a, 2, nchar(a))) )
	sampleStates <- matrix( 0, nrow=length(phylo$tip.label), length(sampleStatesAnnotations ) )
	rownames(sampleStates) <- phylo$tip.label
	for (i in 1:length(phylo$tip.label))
	{
		k <- which( annotations[i] == sampleStatesAnnotations )
		sampleStates[i,k] <- 1
	}
	sampleStates
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


.extant.at.height <- function(h, tree)
{
	return( which( tree$heights <= h & tree$parentheights > h)  )
}

rescale.binaryDatedTree <- function(bdt, ef)
{
#~ bdt : a binaryDatedTree
#~ ef : numeric factory by which node heights are altered 
#~ returns a new binaryDatedTree instance
	ini  <- (bdt$n+1):(length(bdt$heights))
	bdt$heights[ini] <- ef * bdt$heights[ini]
	while( sum(bdt$heights[bdt$daughters[ini,2] ] > bdt$heights[ini] ) > 0 | sum(bdt$heights[bdt$daughters[ini,1] ] >   bdt$heights[ini] ) > 0){
		bdt$heights[ini] <- pmax(bdt$heights[bdt$daughters[ini,1] ],   bdt$heights[ini] )
		bdt$heights[ini] <- pmax(bdt$heights[bdt$daughters[ini,2] ],   bdt$heights[ini] )
	}
	
	for (i in 1:nrow(bdt$edge))
	{
		bdt$edge.length[i] <- bdt$heights[bdt$edge[i,1]] - bdt$heights[bdt$edge[i,2]]
	}
	bdt	
}

##########################
#CALCULATE INTERNAL STATES
.calculate.internal.states.unstructuredModel <- function(tree, globalParms , maxHeight=FALSE)
{
	
		F. <- globalParms$F.; G. <- globalParms$G.; Y. <- globalParms$Y. ; 
	   USE_DISCRETE_FGY <- TRUE
	   FGY_RESOLUTION <- globalParms$FGY_RESOLUTION
	   FGY_H_BOUNDARIES <- globalParms$FGY_H_BOUNDARIES
	   FINITESIZECORRECTIONS <- TRUE
	
	eventTimes <- unique( sort(tree$heights) )
	tree$maxHeight <- maxHeight
	if (maxHeight) { 
		eventTimes <- eventTimes[eventTimes<=maxHeight]}
	S <- 1
	L <- 0
	
	FGY_INDEX <- 1 
	
	get.fgy <- function(h) 
	{
		.Y		<- globalParms$Y.(h)
		.F		<- globalParms$F.(h)
		.G		<- globalParms$G.(h)
		list(.F=.F, .G=.G, .Y=.Y, FGY_INDEX=NA)
	}
#~ 	get.fgy <- function(h,t=NULL)
#~ 	{
#~ 		if (globalParms$USE_DISCRETE_FGY)
#~ 		{
#~ 			while (FGY_INDEX < globalParms$FGY_RESOLUTION &  h > globalParms$FGY_H_BOUNDARIES[FGY_INDEX])  
#~ 			{
#~ 				FGY_INDEX <- FGY_INDEX+1 ; 
#~ 			}
#~ 			if ( (FGY_INDEX > 1) ) 
#~ 			{
#~ 				while ( h< globalParms$FGY_H_BOUNDARIES[FGY_INDEX-1]) 
#~ 				{
#~ 					if (FGY_INDEX==1) break;
#~ 					FGY_INDEX <- FGY_INDEX-1 ; 
#~ 					#return(dA(h, A, parms))
#~ 				} 
#~ 			}
#~ 			.Y		<- globalParms$Y.(FGY_INDEX)
#~ 			.F		<- globalParms$F.(FGY_INDEX)
#~ 			.G		<- globalParms$G.(FGY_INDEX)
#~ 		} else
#~ 		{
#~ 			.Y		<- globalParms$Y.(t)
#~ 			.F		<- globalParms$F.(t)
#~ 			.G		<- globalParms$G.(t)
#~ 			FGY_INDEX <- NA
#~ 		}
#~ 		list(.F=.F, .G=.G, .Y=.Y, FGY_INDEX=FGY_INDEX)
#~ 	}
	
	#DEBUG
	#with( get.fgy(0), {print(.Y); print(.F)  })
	#browser()
	
	.dL.unstructuredModel <- function(h, L, parms, FGY_INDEX, ...){
		# conditional on no coalescent
		t <- parms$treeT - h
		A <- parms$A
		fgy <- get.fgy(h,t)
		with( fgy, 
		{
			.F * (A / .Y)^2# *  (A-1) / max((.Y-1), 0.01)
			#.F * (A / .Y) *  (A-1) / max((.Y-1), 0.01)
		}) -> dL 
		return( list(dL, FGY_INDEX=fgy$FGY_INDEX) )
	}
	.solve.L <- function(h0, h1, L,  A0 ) 
	{
		parameters 		<- list(treeT = tree$maxSampleTime, m = tree$m, A=A0)
			tryCatch({
						out0 <- ode(y=L, times=c(h0, h1), func=.dL.unstructuredModel, parms=parameters, FGY_INDEX=FGY_INDEX, method=globalParms$INTEGRATIONMETHOD ) 
					}, error = function(e) browser() )
			L1 <- unname( out0[nrow(out0),2] )
			FGY_INDEX <- unname( out0[nrow(out0),3] )
			return( list( L1, FGY_INDEX) )
	}
	
	
	
	for (ih in 1:(length(eventTimes)-1)){
		h0 <- eventTimes[ih]
		h1 <- eventTimes[ih+1]
		#get A0, process new samples, calculate lstate for new lines
		extantLines <- .extant.at.height(h0, tree)
		A0 <- length(extantLines)
		
		o <- .solve.L(h0, h1, L,  A0 ) 
		L <- o[[1]]
		FGY_INDEX <- o[[2]]
		
		newNodes <- which( tree$heights == h1)
		newNodes <- newNodes[newNodes > length(tree$sampleTimes)]
		if (length(newNodes) > 0)
		{
#~ 			fgy<- get.fgy(h1,  bdt$maxSampleTime - h1 )
			fgy<- get.fgy(h1 )
			.Y <- fgy$.Y
			.F <- fgy$.F
			FGY_INDEX <- fgy$FGY_INDEX
			
			S <- exp(-L)
			
			#TODO option to return -Inf in this situation: 
			#tryCatch( { if (sum(.Y) < length(extantLines) ) {S <- 0; } }, error = function(e) browser() )
			#print(paste(sum(.Y),length(extantLines)));
			
			for (alpha in newNodes){
				X1 <- A0 / .Y
				# option for finite size corrections:
				#if (FINITESIZECORRECTIONS){
				#	X1 <- X1 * .Y / (max(.Y,1.01)-1)
				#}
				
				tree$lstates[alpha,] <- 1
				tree$mstates[alpha,] <- 1
				ratekl <- .F * (A0 / .Y)^2# *  (A0-1) / max((.Y-1), 0.01)
				#ratekl <- .F * (A0 / .Y) *  (A0-1) / max((.Y-1), 0.01)
				tree$coalescentRates[alpha] <- ratekl
				tree$coalescentSurvivalProbability[alpha] <- S
			}
			
			L<-0
		}
	}
	return( tree )
}


.calculate.internal.states <- function(tree, fgyParms,  censorAtHeight=FALSE, forgiveAgtY = 0.20,   INTEGRATIONMETHOD = 'rk4'){
	F. <- fgyParms$F.; G. <- fgyParms$G.; Y. <- fgyParms$Y. ; 
	FGY_RESOLUTION <- fgyParms$FGY_RESOLUTION
	
#~ 	get.fgy <- function(h)
#~ 	{
#~ 			FGY_INDEX <- min( 1+floor( FGY_RESOLUTION * h / tree$maxHeight ), FGY_RESOLUTION)
#~ 			.Y		<- fgyParms$Y.(FGY_INDEX)
#~ 			.F		<- fgyParms$F.(FGY_INDEX)
#~ 			.G		<- fgyParms$G.(FGY_INDEX)
#~ 		list(.F=.F, .G=.G, .Y=.Y)
#~ 	}
	get.fgy <- function(h)
	{
			.Y		<- fgyParms$Y.(h)
			.F		<- fgyParms$F.(h)
			.G		<- fgyParms$G.(h)
		list(.F=.F, .G=.G, .Y=.Y)
	}
	
	
	# construct forcing timeseries for ode's
#~ 	times <- fgyParms$FGY_H_BOUNDARIES
	heights <- sort( fgyParms$FGY_H_BOUNDARIES )
	heights <- heights[heights <= tree$maxHeight]
	heights <- heights[heights >= 0]
	fgymat <- t( sapply( heights, function(h) 
	  with(get.fgy(h), {
		c( as.vector(.F), as.vector(.G), as.vector(.Y) )
	  })
	) )
#~ 	fgymat <- t( sapply( 1:FGY_RESOLUTION, function(i) 
#~ 	  c( as.vector(fgyParms$F_DISCRETE[[i]]) , 
#~ 	    as.vector( fgyParms$G_DISCRETE[[i]] ), 
#~ 	    fgyParms$Y_DISCRETE[[i]]) ) 
#~ 	)
#~ browser()
	fgymat <- pmax(fgymat, 0)

	.solve.Q.A.L <- function(h0, h1, A0, L0)
	{ # uses C implementation
		Q0 <- diag(tree$m)
		parameters 	<- c(tree$m, tree$maxHeight, length(heights), sum(A0), as.vector(fgymat))
		y0 <- c( as.vector( Q0), A0,  L0 ) #
		#o <- ode(y=y0, c(h0,h1), func = "dQAL", parms = parameters, dllname = "dQAL-6.3", initfunc = "initfunc", method=INTEGRATIONMETHOD )
		o <- ode(y=y0, c(h0,h1), func = "dQAL", parms = parameters, dllname = "rcolgem", initfunc = "initfunc", method=INTEGRATIONMETHOD )
		Q1 		<- t( matrix(  abs(o[nrow(o),2:(1 + tree$m^2)]) , nrow=tree$m) ) #NOTE the transpose
		A1 <- o[nrow(o), (1 + tree$m^2 + 1):(1 + tree$m^2 +  tree$m)]
		L1 <- o[nrow(o), ncol(o)]
		return ( list(  unname(Q1), unname(A1),  unname(L1) ) ) #
	}
	
	eventTimes <- unique( sort(tree$heights) )
	if (censorAtHeight) { 
		eventTimes <- eventTimes[eventTimes<=censorAtHeight]}
	S <- 1
	L <- 0
	
	#debugging symbols
	tree$A <- c() #h->A
	tree$lnS <- c()
	tree$lnr <- c()
	tree$ih <- c()
	
	
	# helper function for updating ancestral node and likelihood terms
	.update.alpha <- function(u,v, tree, fgy)
	{
		.F <- fgy$.F
		.G <- fgy$.G
		.Y <- fgy$.Y
		{
			.Y <- pmax(A, .Y)
			FklXpuk_Yk <- (.F * tree$ustates[u,]/.Y)
			FklXpvk_Yk <- (.F * tree$ustates[v,]/.Y)
			FklXpuk_Yk[is.nan(FklXpuk_Yk)] <- 0
			FklXpvk_Yk[is.nan(FklXpvk_Yk)] <- 0
			vk_Yk <- pmin(pmax(tree$ustates[v,]/.Y, 0),1); vk_Yk[is.nan(vk_Yk)] <- 0
			uk_Yk <- pmin(pmax(tree$ustates[u,]/.Y, 0),1); uk_Yk[is.nan(uk_Yk)] <- 0
			ratekl <- FklXpuk_Yk %*% vk_Yk + FklXpvk_Yk %*% uk_Yk
		}
		
		coalescentRate <- max( sum(ratekl) , 0)
		if (is.nan(L))
		{
			L <- (h1 - h0) * sum(ratekl) 
		}
		
		logCoalescentSurvivalProbability <- -L
		if (sum(ratekl)==0) {ratekl <- rep(1/tree$m, tree$m) * 1e-6}
		# definitions of alpha state
		p_a <- as.vector( ratekl / sum(ratekl) )
		
		list( p_a, logCoalescentSurvivalProbability, coalescentRate) 
	}
	
	.fsc.extantLines <- function(alpha, extantLines, A, tree )
	{ #	finite size corrections for lines not involved in coalescent
		p_a     <- tree$lstates[alpha,]
		p_i_mat <- tree$mstates[extantLines,]
		A_mat   <- t( matrix( A, nrow=tree$m, ncol=length(extantLines) )  )
		p_a_mat <- t( matrix(p_a, nrow=tree$m, ncol=length(extantLines)) )
		rho_mat <- A_mat /  (A_mat - p_i_mat)
		rterms  <- p_a_mat / (A_mat - p_i_mat)
		lterms  <- rho_mat %*% p_a
		lterms <- matrix(lterms, nrow=length(lterms), ncol=tree$m) #
		new_p_i <- p_i_mat * (lterms - rterms)
		new_p_i <- new_p_i / rowSums(new_p_i)
		corrupted <- is.nan(rowSums(new_p_i))
		new_p_i[corrupted,] <- p_i_mat[corrupted,]
		tree$mstates[extantLines,] <- new_p_i
		tree
	}
	
	
	for (ih in 1:(length(eventTimes)-1)){
		h0 <- eventTimes[ih]
		h1 <- eventTimes[ih+1]
		fgy <- get.fgy(h1)
		
		#get A0, process new samples, calculate state of new lines
		extantLines <- .extant.at.height(h0, tree)
		if (length(extantLines) > 1 ){
			A0 <- colSums(tree$mstates[extantLines,])
		} else if (length(extantLines)==1){ A0 <- tree$mstates[extantLines,] }
		out <- .solve.Q.A.L(h0, h1, A0,  L)
		Q <- out[[1]]
		A <- out[[2]]
#~ print(out)
#~ if (out[[3]] > 1000) browser()
		L <- out[[3]]

		# clean output
		if (is.nan(L)) {L <- Inf}
		if (sum(is.nan(Q)) > 0) Q <- diag(length(A))
		if (sum(is.nan(A)) > 0) A <- A0
		
		#update mstates
		if (length(extantLines) > 1)
		{
			tree$mstates[extantLines,] <- t( t(Q) %*% t(tree$mstates[extantLines,])  )
			tree$mstates[extantLines,] <- abs(tree$mstates[extantLines,]) / rowSums(abs(tree$mstates[extantLines,]))
			#recalculate A
			A <- colSums(tree$mstates[extantLines,])
		}
		else{
			tree$mstates[extantLines,] <- t( t(Q) %*% tree$mstates[extantLines,] )
			tree$mstates[extantLines,] <- abs(tree$mstates[extantLines,]) / sum(abs(tree$mstates[extantLines,]))
			#recalculate A
			A <- tree$mstates[extantLines,]
		}
		
		
		#if applicable: update ustate & calculate lstate of new line 
		newNodes <- which( tree$heights == h1)
		newNodes <- newNodes[newNodes > length(tree$sampleTimes)] # <- internal nodes
		
		with(fgy, 
		{
			# check inputs
			if (sum(is.na(.Y)) > 0) {warning('Model returned NA population values', .F, ' ', .G, ' ', .Y); 
				.Y[is.na(.Y)] <- 0}
			if (sum(is.na(.G)) > 0) {warning('Model returned NA population values', .F, ' ', .G, ' ', .Y); 
				.G[is.na(.G)] <- 0}
			if (sum(is.na(.F)) > 0) {warning('Model returned NA population values', .F, ' ', .G, ' ', .Y); 
				.F[is.na(.F)] <- 0}
				
			if (sum(is.nan(.Y)) > 0) {warning('Model returned NaN population values', .F, ' ', .G, ' ', .Y); 
				.Y[is.nan(.Y)] <- 0}
			if (sum(is.nan(.G)) > 0) {warning('Model returned NaN population values', .F, ' ', .G, ' ', .Y); 
				.G[is.nan(.G)] <- 0}
			if (sum(is.nan(.F)) > 0) {warning('Model returned NaN population values', .F, ' ', .G, ' ', .Y); 
				.F[is.nan(.F)] <- 0}
			
			if (sum(.Y < 0) > 0) {warning('Model returned negative population values', .F, ' ', .G, ' ', .Y); .Y <- pmax(.Y,0)}
			if (sum(.G < 0) > 0) {warning('Model returned negative population values', .F, ' ', .G, ' ', .Y); .G <- pmax(.G, 0)}
			if (sum(.F < 0) > 0) {warning('Model returned negative population values', .F, ' ', .G, ' ', .Y); .F <- pmax(.F, 0)}
			if ( (sum(.Y) < length(extantLines)) & (!forgiveAgtY) ) { L <- Inf }
			else if ( (sum(.Y) < length(extantLines)) & (length(extantLines)/length(tree$tip.label)) > forgiveAgtY) { L <- Inf }
			
			#for (alpha in newNodes){
			if (length(newNodes)==1)
			{
				#TODO should rewrite to use helper functions
				alpha <- newNodes
				u <- tree$daughters[alpha,1]
				v <- tree$daughters[alpha,2]
				{
					tree$ustates[u,] <- tree$mstate[u,]
					tree$ustates[v,] <- tree$mstate[v,]
					
					.Y <- pmax(A, .Y)
					FklXpuk_Yk <- (.F * tree$ustates[u,]/.Y)
					FklXpvk_Yk <- (.F * tree$ustates[v,]/.Y)
					FklXpuk_Yk[is.nan(FklXpuk_Yk)] <- 0
					FklXpvk_Yk[is.nan(FklXpvk_Yk)] <- 0
					vk_Yk <- pmin(pmax(tree$ustates[v,]/.Y, 0),1); vk_Yk[is.nan(vk_Yk)] <- 0
					uk_Yk <- pmin(pmax(tree$ustates[u,]/.Y, 0),1); uk_Yk[is.nan(uk_Yk)] <- 0
					ratekl <- FklXpuk_Yk %*% vk_Yk + FklXpvk_Yk %*% uk_Yk
				}
				
				
				tree$coalescentRates[alpha] <- max( sum(ratekl) , 0)
				if (is.nan(L))
				{
					L <- (h1 - h0) * sum(ratekl) 
				}
				tree$coalescentSurvivalProbability[alpha] <- exp(-L)
				tree$logCoalescentSurvivalProbability[alpha] <- -L
				if (tree$coalescentSurvivalProbability[alpha]==0) {warning('Zero likelihood. Try setting integrationMethod=adams')}
#~ if (tree$coalescentSurvivalProbability[alpha]==0) browser()
				if (sum(ratekl)==0) {ratekl <- rep(1/tree$m, tree$m) * 1e-6}
				# definitions of alpha state
				tree$lstates[alpha,] <- ratekl / sum(ratekl)
				tree$mstates[alpha,] <- ratekl / sum(ratekl)
				
				# debug
				tree$A <- rbind( tree$A, A)
				tree$lnS <- c( tree$lnS, -L )
				tree$lnr <- c( tree$lnr, log(sum(ratekl)) )
				tree$ih <- c( tree$ih, h1)
				# finite size corrections for lines not involved in coalescent
				{ #
					p_a     <- tree$lstates[alpha,]
					p_i_mat <- tree$mstates[extantLines,]
					A_mat   <- t( matrix( A, nrow=tree$m, ncol=length(extantLines) )  )
					p_a_mat <- t( matrix(p_a, nrow=tree$m, ncol=length(extantLines)) )
					rho_mat <- A_mat /  (A_mat - p_i_mat)
					rterms  <- p_a_mat / (A_mat - p_i_mat)
					lterms  <- rho_mat %*% p_a
					lterms <- matrix(lterms, nrow=length(lterms), ncol=tree$m) #
					new_p_i <- p_i_mat * (lterms - rterms)
					new_p_i <- new_p_i / rowSums(new_p_i)
					corrupted <- is.nan(rowSums(new_p_i))
					new_p_i[corrupted,] <- p_i_mat[corrupted,]
					tree$mstates[extantLines,] <- new_p_i
				}
			} else if (length(newNodes) > 1){
				#TODO should check that this is definite polytomy
				daughters <- setdiff( unique( tree$daughters[newNodes,] ), newNodes)
				uv_corates <- c()
				uv_survProbs <- c()
				uv_alphaStates <- c()
				tree$ustates[daughters,] <- tree$mstates[daughters,]
				for (u_i in 1:(length(daughters)-1)) for (v_i in (u_i+1):length(daughters))
				{
					u <- daughters[u_i]
					v <- daughters[v_i]
					o <- .update.alpha(u, v, tree, fgy) #list( p_a, logCoalescentSurvivalProbability, coalescentRate) 
					uv_corates <- c( uv_corates, o[[3]] )
					uv_survProbs <- c( uv_survProbs, exp( o[[2]] ) )
					uv_alphaStates <- rbind ( uv_alphaStates, o[[1]] )
				}
				corate <- mean(uv_corates)
				S <- mean(uv_survProbs)
				p_a <- colMeans( uv_alphaStates)
				p_a <- p_a / sum(p_a) 
				for (alpha in newNodes)
				{
					# debug
					tree$A <- rbind( tree$A, A)
					tree$lnS <- c( tree$lnS, -L )
					tree$lnr <- c( tree$lnr, log(corate) )
#~ if( is.infinite( log(corate) ) ) browser()
					tree$ih <- c( tree$ih, h1)
					# update alpha states
					tree$lstates[alpha,] <- p_a
					tree$mstates[alpha,] <- tree$lstates[alpha,]
					if ((tree$parentheights[alpha]-tree$heights[alpha]) == 0) tree$ustates[alpha,] <- tree$lstates[alpha,]
					# likelihood stuff
					tree$coalescentRates[alpha] <- corate
					tree$coalescentSurvivalProbability[alpha] <- S
					tree$logCoalescentSurvivalProbability[alpha] <- log(S)
#~ if (S==0) browser()
				}
				tree <- .fsc.extantLines(alpha, extantLines, A, tree )
#~ print(c(h0, h1, corate , S)) 
#~ print(daughters)
			}
			tree
		}) -> tree
		# if coalescent occurred, reset cumulative hazard function
		if (length(newNodes) > 0) {
		  L <- 0 }
	}
	 
	return(tree)
} 



.start.discrete.rates <- function(fgyResolution, maxSampleTime, globalParms) 
{ #version that does not generate global variables
	F. <- globalParms$F.; G. <- globalParms$G.; Y. <- globalParms$Y. ; 
	globalParms$FGY_RESOLUTION		<- fgyResolution
	#~ speed up calculation of FGY by discretizing & pre-caching
	globalParms$F.bak 				<- F.
	globalParms$G.bak 				<- G.
	globalParms$Y.bak 				<- Y.
	globalParms$USE_DISCRETE_FGY 	<- TRUE
	globalParms$maxHeight			<- maxSampleTime #NOTE: This does not allow FGY defined for t<0
	globalParms$FGY_H_BOUNDARIES 	<- seq(0, maxSampleTime, length.out = fgyResolution) 
#~ 	globalParms$FGY_H_BOUNDARIES 	<- globalParms$FGY_H_BOUNDARIES + globalParms$FGY_H_BOUNDARIES[2]/2 #TODO is there a better way to do this offset? 
	globalParms$FGY_INDEX 			<- 1 #update in desolve:ode
	globalParms$F_DISCRETE 			<- lapply( globalParms$FGY_H_BOUNDARIES, function(h) { F.(maxSampleTime-h) })
	globalParms$G_DISCRETE 			<- lapply( globalParms$FGY_H_BOUNDARIES, function(h) { G.(maxSampleTime-h) })
	globalParms$Y_DISCRETE 			<- lapply( globalParms$FGY_H_BOUNDARIES, function(h) { Y.(maxSampleTime-h) })
	globalParms$F. 					<- function(FGY_INDEX) { globalParms$F_DISCRETE[[FGY_INDEX]] } #note does not actually use arg t
	globalParms$G. 					<- function(FGY_INDEX) { globalParms$G_DISCRETE[[FGY_INDEX]] }
	globalParms$Y. 					<- function(FGY_INDEX) { globalParms$Y_DISCRETE[[FGY_INDEX]] }
	globalParms
}

#############################
.end.discrete.rates <- function(globalParms){
	#reset fgy functions
	# does not use global variables
	globalParms$F. <- globalParms$F.bak
	globalParms$G. <- globalParms$G.bak
	globalParms$Y. <- globalParms$Y.bak
	globalParms
}
#############################
#~ .gen.discrete.fgy.parameters <- function(FGY, fgyResolution, maxSampleTime, maxHeight ) 
#~ { #~ speed up calculation of FGY by discretizing & pre-caching
#~ 	parms <- list()
#~ 	F. <- FGY$F.; G. <- FGY$G.; Y. <- FGY$Y. ; 
#~ 	parms$FGY_RESOLUTION		<- fgyResolution
#~ 	parms$maxHeight				<- maxHeight #
#~ 	parms$get.index <- function(h) min(1 + parms$FGY_RESOLUTION * (h) / ( maxHeight ), parms$FGY_RESOLUTION )
#~ 	parms$FGY_H_BOUNDARIES 		<- seq(0, maxHeight, length.out = fgyResolution) #
#~ 	parms$F_DISCRETE 			<- lapply( parms$FGY_H_BOUNDARIES, function(h) { F.(maxSampleTime-h) })
#~ 	parms$G_DISCRETE 			<- lapply( parms$FGY_H_BOUNDARIES, function(h) { G.(maxSampleTime-h) })
#~ 	parms$Y_DISCRETE 			<- lapply( parms$FGY_H_BOUNDARIES, function(h) { pmax( 1e-16, Y.(maxSampleTime-h) ) })
#~ ##~ 	parms$F. 					<- function(FGY_INDEX) { parms$F_DISCRETE[[FGY_INDEX]] } #note does not actually use #arg t
#~ ##~ 	parms$G. 					<- function(FGY_INDEX) { parms$G_DISCRETE[[FGY_INDEX]] }
#~ ##~ 	parms$Y. 					<- function(FGY_INDEX) { parms$Y_DISCRETE[[FGY_INDEX]] }
#~ 	fgyParms$F. 					<- function(h) { fgyParms$F_DISCRETE[[parms$get.index(h)]] } #
#~ 	fgyParms$G. 					<- function(h) { fgyParms$G_DISCRETE[[parms$get.index(h)]] }
#~ 	fgyParms$Y. 					<- function(h) { fgyParms$Y_DISCRETE[[parms$get.index(h)]] }
#~ 	parms
#~ }
#############################
# CALCULATE LIKELIHOOD
#' @export
#~ coalescent.log.likelihood.fgy <- function(bdt, FGY, integrationMethod = 'rk4',  censorAtHeight=FALSE, forgiveAgtY=.2, fgyResolution = 2000)
#~ {
#~ # bdt : binaryDatedTree instance
#~ 	fgyParms <- .gen.discrete.fgy.parameters(FGY, fgyResolution, bdt$maxSampleTime, bdt$maxHeight ) 
#~ 	if (ncol( bdt$sampleStates ) ==1) 
#~ 	  tree <- .calculate.internal.states.unstructuredModel(bdt, maxHeight=censorAtHeight, globalParms = fgyParms)
#~ 	else
#~ 	  tree <- .calculate.internal.states(bdt, fgyParms,  censorAtHeight=censorAtHeight, forgiveAgtY = forgiveAgtY,   INTEGRATIONMETHOD = integrationMethod)
#~ 	
#~ 	i<- (length(tree$sampleTimes)+1):(tree$Nnode + length(tree$tip.label))
#~ 	if (censorAtHeight) { 
#~ 		internalHeights <- tree$heights[(length(tree$tip.label)+1):length(tree$heights)]
#~ 		i <- i[internalHeights <= censorAtHeight] 
#~ 	}
#~ 	ll <- sum( log(tree$coalescentRates[i]) ) + sum( tree$logCoalescentSurvivalProbability[i] )#sum( log(tree$coalescentSurvivalProbability[i]) ) 
#~ 	if (censorAtHeight) { ll<- tree$Nnode *  ll/length(i)}
#~ 	if (is.nan(ll) | is.na(ll) ) ll <- -Inf
#~ 	#return(list( ll, tree)) #for debugging
#~ 	ll
#~ }

coalescent.log.likelihood.fgy <- function(bdt, times, births, migrations, demeSizes, integrationMethod = 'rk4',  censorAtHeight=FALSE, forgiveAgtY=.2, returnTree=FALSE)
{
# bdt : binaryDatedTree instance
# note births & migrations should be rates in each time step
	maxtime <- times[length(times)]
	mintime <- times[1]
	if (mintime > (bdt$maxSampleTime - bdt$maxHeight)) {warning('Root of tree occurs before earliest time on time axis'); return(-Inf) }
	fgyParms <- list()
	fgyParms$FGY_RESOLUTION		<- length(times)
	fgyParms$maxHeight				<- bdt$maxHeight #
	# reverse order (present to past): 
	fgyParms$F_DISCRETE 			<- lapply(length(births):1, function(k)  births[[k]] )
	fgyParms$G_DISCRETE 			<- lapply(length(migrations):1, function(k)  migrations[[k]] )
	fgyParms$Y_DISCRETE 			<- lapply(length(demeSizes):1, function(k)  demeSizes[[k]] )
	# the following line accounts for any discrepancies between the maximum time axis and the last sample time
	fgyParms$hoffset = hoffset <- maxtime - bdt$maxSampleTime
	if (hoffset < 0) stop( 'Time axis does not cover the last sample time' )
	fgyParms$get.index <- function(h) min(1 + fgyParms$FGY_RESOLUTION * (h + hoffset) / ( maxtime - mintime ), fgyParms$FGY_RESOLUTION )
	fgyParms$F. 					<- function(h) { fgyParms$F_DISCRETE[[fgyParms$get.index(h)]] } #
	fgyParms$G. 					<- function(h) { fgyParms$G_DISCRETE[[fgyParms$get.index(h)]] }
	fgyParms$Y. 					<- function(h) { fgyParms$Y_DISCRETE[[fgyParms$get.index(h)]] }
	fgyParms$FGY_H_BOUNDARIES 		<- bdt$maxSampleTime - times[length(times):1] 
	
	if (ncol( bdt$sampleStates ) ==1) 
	  tree <- .calculate.internal.states.unstructuredModel(bdt, maxHeight=censorAtHeight, globalParms = fgyParms)
	else
	  tree <- .calculate.internal.states(bdt, fgyParms,  censorAtHeight=censorAtHeight, forgiveAgtY = forgiveAgtY,   INTEGRATIONMETHOD = integrationMethod)
	
	i<- (length(tree$sampleTimes)+1):(tree$Nnode + length(tree$tip.label))
	if (censorAtHeight) { 
		internalHeights <- tree$heights[(length(tree$tip.label)+1):length(tree$heights)]
		i <- i[internalHeights <= censorAtHeight] 
	}
	ll <- sum( log(tree$coalescentRates[i]) ) + sum( tree$logCoalescentSurvivalProbability[i] )#sum( log(tree$coalescentSurvivalProbability[i]) ) 
	if (censorAtHeight) { ll<- tree$Nnode *  ll/length(i)}
	if (is.nan(ll) | is.na(ll) ) ll <- -Inf
	#~ 	ifelse(returnTree, list( ll, tree), ll) #for debugging
	if(returnTree)
		return(list(ll, tree))
	else
		return(ll)
}


pseudoLogLikelihood0.fgy <-  function(bdt, times, births, migrations, demeSizes, forgiveAgtY=.2, integrationMethod = 'adams', rescaleTree=FALSE, returnRescaleFactor=FALSE, ks_objfun = FALSE)
{
# note births & migrations should be rates in each time step
#~ ks_obfjun : if TRUE, returns -KS statistic (NLFT and node heights)
	sampleStates <- bdt$sampleStates
	sampleTimes <- bdt$sampleTimes
	n <- length(sampleTimes)
	maxSampleTime <- max(sampleTimes)
	sampleHeights <- maxSampleTime - sampleTimes 
	ix <- sort(sampleHeights, index.return=TRUE)$ix
	sortedSampleHeights <- sampleHeights[ix]
	sortedSampleStates <- sampleStates[ix,]
	uniqueSortedSampleHeights <- unique(sortedSampleHeights)
	
	maxtime <- times[length(times)]
	mintime <- times[1]
	maxHeight <- maxSampleTime -  mintime
	
	fgyParms <- list()
	fgyParms$FGY_RESOLUTION		<- length(times)
	# reverse order (present to past): 
	fgyParms$F_DISCRETE 			<- lapply(length(births):1, function(k)  births[[k]] )
	fgyParms$G_DISCRETE 			<- lapply(length(migrations):1, function(k)  migrations[[k]] )
	fgyParms$Y_DISCRETE 			<- lapply(length(demeSizes):1, function(k)  demeSizes[[k]] )
	# the following line accounts for any discrepancies between the maximum time axis and the last sample time
	fgyParms$hoffset = hoffset <- maxtime - bdt$maxSampleTime
	if (hoffset < 0) stop( 'Time axis does not cover the last sample time' )
	fgyParms$get.index <- function(h) min(1 + fgyParms$FGY_RESOLUTION * (h + hoffset) / ( maxtime - mintime ), fgyParms$FGY_RESOLUTION )
	fgyParms$F. 					<- function(h) { fgyParms$F_DISCRETE[[fgyParms$get.index(h)]] } #
	fgyParms$G. 					<- function(h) { fgyParms$G_DISCRETE[[fgyParms$get.index(h)]] }
	fgyParms$Y. 					<- function(h) { fgyParms$Y_DISCRETE[[fgyParms$get.index(h)]] }
	fgyParms$FGY_H_BOUNDARIES 		<- bdt$maxSampleTime - times[length(times):1] 

	F. <- fgyParms$F.; G. <- fgyParms$G.; Y. <- fgyParms$Y. ; 
	FGY_RESOLUTION <- fgyParms$FGY_RESOLUTION
	
	get.fgy <- function(h)
	{
			.Y		<- fgyParms$Y.(h)
			.F		<- fgyParms$F.(h)
			.G		<- fgyParms$G.(h)
		list(.F=.F, .G=.G, .Y=.Y)
	}
	
	# construct forcing timeseries for ode's
	heights <- sort( fgyParms$FGY_H_BOUNDARIES )
	heights <- heights[heights <= maxHeight]
	heights <- heights[heights >= 0]
	fgyParms$heights <- heights
	fgymat <- t( sapply( fgyParms$heights, function(h) 
	  with(get.fgy(h), {
		c( as.vector(.F), as.vector(.G), as.vector(.Y) )
	  })
	) )
	fgymat <- pmax(fgymat, 0)
	
	inheights <- sort( bdt$heights[(n+1):length(bdt$heights)] )
	forgiveHeight <- inheights[ max(1, length(inheights) - floor(forgiveAgtY*length(inheights)))]
	
	# solve for A
	cumSortedSampleStates <-  sapply( 1:m, function(k) cumsum(sortedSampleStates[,k]) )
	cumSortedNotSampledStates <-  t(cumSortedSampleStates[n,] - t(cumSortedSampleStates) )
	#nsy.index <- approxfun( sortedSampleHeights , 1:n , method='constant', rule=2)
	nsy.index <-  approxfun(
	   sort( jitter( sortedSampleHeights
	     , factor=max(1e-6, sortedSampleHeights[length(sortedSampleHeights)]/1e6)) )
	   , 1:n , method='constant', rule=2)
	not.sampled.yet <- function(h)
	{
		cumSortedNotSampledStates[ nsy.index(h), ]
	}
	nsymat <- t( sapply( fgyParms$heights, function(h) 
		not.sampled.yet(h)
	) )
	dA <- function(h, A, parms, ...)
	{
		nsy <- not.sampled.yet(h) 
		with(get.fgy(h), 
		{ 
			A_Y <- (A-nsy) / .Y
			if (h < forgiveHeight){
				if ( sum(A_Y > 1)>0 ){
					return(list( rep(0, length(A)) ))
				}
			}
			csFpG 	<- colSums( .F + .G )
			list( .G %*% A_Y - csFpG * A_Y + (.F %*% A_Y) * pmax(1-A_Y, 0) )
		})
	}
	AIntervals <- c(sortedSampleHeights[2:length(uniqueSortedSampleHeights)], maxHeight)
	h0 	<- 0
	sampled.at.h <- function(h) which(sortedSampleHeights==h)
	extantLines <- sampled.at.h(h0)
	
	#~ 	AplusNotSampled <- ode( y = colSums(sortedSampleStates), times = haxis, func=dA, parms=NA, method = integrationMethod)[, 2:(m+1)]
	parameters 	<- c(m, maxHeight, length(heights), as.vector(fgymat), as.vector(nsymat))
	# could possibly speed this up by calling c func:
	if (length(rescaleTree)==2 & is.numeric(rescaleTree) ) {
		stretch_lb <- rescaleTree[1]
		stretch_ub <- rescaleTree[2]
		odeHeights <- seq(0, max(inheights), length.out=2000)
		dh <- odeHeights[2] - odeHeights[1]
		AplusNotSampled <- ode(y=colSums(sortedSampleStates), times = odeHeights, func = "dCA", parms = parameters, dllname = "rcolgem", initfunc = "initfunc_dCA", method=integrationMethod )[, 2:(m+1)]
		.stretchInternalNodeHeights <- function(inheights, ef ) 
		{
			# note : does not check for consistency with sample times
			#(inheights-min(inheights)) * ef + min(inheights)
			inheights * ef
		}
		objfun <- function(ef){
			new_inheights <- .stretchInternalNodeHeights(inheights, ef)
			m_dAs <- sapply(new_inheights, function(h)
			{
				ih <- min( length(odeHeights),  1 + floor( h / dh )  )
				-sum( dA(h, AplusNotSampled[ih,] , NA )[[1]] )  
			})
			return( sum(log(m_dAs)) )
		}
		rsApNS <- rowSums( AplusNotSampled )
		# CDF of node heights: 
		pCo <- approxfun( odeHeights, (rsApNS[1] - rsApNS) / (rsApNS[1] - rsApNS[length(rsApNS)] ),  rule =2 )
		objfun.ks <- function(ef=1){
			new_inheights <- .stretchInternalNodeHeights(inheights, ef)
			return( ks.test(x = new_inheights, y = pCo)$statistic )
		}
		# minimise KS statistic: 
		efoptimise <- optimise(f = objfun.ks, interval = rescaleTree
		  , maximum = FALSE
		  , tol = .Machine$double.eps^0.25)
		if (returnRescaleFactor){
			if (ks_objfun){
				return( list( -efoptimise$objective , efoptimise$minimum ) )
			} else{
				return( list( objfun(efoptimise$minimum) , efoptimise$minimum ) )
			}
			
		} else{
			if (ks_objfun){
				return( -efoptimise$objective  )
			} else{
				return(objfun(efoptimise$minimum) )
			}
		}
	} else if (!rescaleTree) {
		if (ks_objfun){
			odeHeights <- seq(0, max(inheights), length.out=2000)
			dh <- odeHeights[2] - odeHeights[1]
			AplusNotSampled <- ode(y=colSums(sortedSampleStates), times = odeHeights, func = "dCA", parms = parameters, dllname = "rcolgem", initfunc = "initfunc_dCA", method=integrationMethod )[, 2:(m+1)]
			rsApNS <- rowSums( AplusNotSampled )
			# CDF of node heights: 
			pCo <- approxfun( odeHeights, (rsApNS[1] - rsApNS) / (rsApNS[1] - rsApNS[length(rsApNS)] ),  rule =2 )
#~ plot( inheights, 1:length(inheights) / length(inheights) )
#~ lines( odeHeights, pCo(odeHeights), col='red')
#~ browser()
			return( -ks.test(x = inheights, y = pCo)$statistic )
		} else{
			AplusNotSampled <- ode(y=colSums(sortedSampleStates), times = c(0,inheights), func = "dCA", parms = parameters, dllname = "rcolgem", initfunc = "initfunc_dCA", method='adams' )[2:(1+length(inheights)), 2:(m+1)]
			m_dAs <- sapply(1:nrow(AplusNotSampled), function(ih) -sum(dA(inheights[ih], AplusNotSampled[ih,] , NA)[[1]])  ) 
			return( sum(log(m_dAs)) )
		}
	} else{
		stop('pseudoLogLikelihood0.fgy : rescaleTree must be FALSE or a 2-vector with scaling limits')
	}
}

solve.model.unstructured <- function(t0,t1, x0, births,  deaths, nonDemeDynamics, parms, timeResolution=1000, integrationMethod = 'rk4')
{
	m <- 1
	mm <- length(nonDemeDynamics)
	demeNames <- names(births)
	nonDemeNames <- names(nonDemeDynamics)
	if (length(x0)!=m + mm) stop('initial conditons incorrect dimension', x0, m, mm) 
	if ( sum( !(c(demeNames, nonDemeNames) %in% names(x0)) )  > 0)  stop('initial conditions vector incorrect names', names(x0), demeNames, nonDemeNames)
	y0 <- c( x0[c(demeNames, nonDemeNames)])
	pbirths <- parse(text=births)
	pmigrations <- NA
	pdeaths <- parse(text=deaths)
	pndd <- sapply(1:mm, function(k) parse(text=nonDemeDynamics[k]) )
	
	dNonDeme <- function(x, t) 
	{
		with(as.list(x, t), 
		  sapply(1:mm, function(k) eval(pndd[k]) )  
		)
	}
	
	times <- seq(t0,t1, length.out=timeResolution)
	dx <- function(t, y, parms, ...) 
	{
		with(as.list(y), 
		{
			f <- eval(pbirths) 
			dxdeme <- f - eval(pdeaths) 
			dxnondeme <- dNonDeme(y, t) 
			names(dxdeme) <- demeNames
			names(dxnondeme) <- nonDemeNames
			list( c(dxdeme, dxnondeme ) )
		})
	}
	ode(y=y0, times, func=dx, parms, method=integrationMethod)
}


#~ if (m==1) tree <- .calculate.internal.states.unstructuredModel(bdt, maxHeight=censorAtHeight, globalParms = fgyParms)
coalescent.log.likelihood.unstructuredModel <- function(bdt, births,  deaths, nonDemeDynamics,  t0, x0, parms=NA, fgyResolution = 2000, integrationMethod = 'euler',  censorAtHeight=FALSE, forgiveAgtY=.2 )
{
	m <- 1
	mm <- length(nonDemeDynamics)
	demeNames <- names(births)
	nonDemeNames <- names(nonDemeDynamics)
	if (length(x0)!=m + mm) stop('initial conditons incorrect dimension', x0, m, mm) 
	if ( sum( !(c(demeNames, nonDemeNames) %in% names(x0)) )  > 0)  stop('initial conditions vector incorrect names', names(x0), demeNames, nonDemeNames)
	y0 <- c( x0[c(demeNames, nonDemeNames)], Lambda=0)
	pbirths <- parse(text=births)
	pmigrations <- NA
	pdeaths <- parse(text=deaths)
	pndd <- sapply(1:mm, function(k) parse(text=nonDemeDynamics[k]) )
	
	dNonDeme <- function(x, t) 
	{
		with(as.list(x, t), 
		  sapply(1:mm, function(k) eval(pndd[k]) )  
		)
	}
	
	
	nodeHeights <- bdt$heights[(length(bdt$tip.label)+1):length(bdt$heights)]
	delta <- c(rep(1, length(bdt$sampleTimes)), rep(-1, length(nodeHeights)))
	delta_times <- c(bdt$maxSampleTime - bdt$sampleTimes, nodeHeights) 
	itimes <- sort(delta_times, index.return=TRUE)$ix
	ltt_times <- sort(delta_times)
	ltt <- cumsum( delta[itimes] )
	n.extant <- function(h)
	{
		#sum(sampleHeights < h) - sum(nodeHeights < h)
		#pmax( -(ltt_times - h), 0)
		i <- which( ltt_times - h < 0, arr.ind=TRUE) #TODO should return n at 0, but returns 1
		if (length(i) == 0) return(1)
		ltt[i[length(i)]]
	}
	
	co.rate <- function( f, Y, A)
	{
		if (A < 2) return(  0 )
		else return( (A * (A-1) /2) * (2*f) / max(2, Y*(Y-1)) )
	}
	
	times <- unique(sort( c( bdt$maxSampleTime - bdt$heights, seq(t0, bdt$maxSampleTime, length.out=fgyResolution) )) )
	dx <- function(t, y, parms, ...) 
	{
		with(as.list(y), 
		{
			f <- eval(pbirths) 
			dxdeme <- f - eval(pdeaths) 
			dxnondeme <- dNonDeme(y, t) 
			names(dxdeme) <- demeNames
			names(dxnondeme) <- nonDemeNames
			dL <- co.rate( f, y[demeNames], n.extant( bdt$maxSampleTime - t ) )
			list( c(dxdeme, dxnondeme, Lambda=unname(dL) ) )
		})
	}
	ox <- ode(y=y0, times, func=dx, parms, method=integrationMethod)
	
	xinterps <- lapply( 2:ncol(ox) , function(k) approxfun( rule=2, bdt$maxSampleTime - ox[,1], ox[,k] ) )
	xinterp <- function(h) {
		x <- sapply( 1:length(x0), function(k) xinterps[[k]](h)  )
		names(x) <- names(x0) 
		x
	}
	heights <- sort(unique( c(bdt$maxSampleTime - bdt$sampleTimes, nodeHeights) ) )
		dLambda <- function(t, y, parms, ...) 
		{
			h <- t
			t <- bdt$maxSampleTime - h
			x <- xinterp( h) 
			Y <- x[demeNames]
			f <- with(as.list(x), eval(pbirths) )
			A <- n.extant(h)
			pco <- (A * (A-1) /2) * (2 / max(2, Y*(Y-1)) )
			if (Y <= 2) pco <- (A * (A-1) /2)
			#if (Y <= 2) pco <- 1
			#if (pco > 1) pco <- 1 #maybe dont do this.. 
			if (A < 2) dL <- 0
			else dL <- pco * f 
			list( c(Lambda=unname( dL )))
		}
		oL  <- ode(y=c(Lambda=0), heights, func=dLambda, parms, method=integrationMethod )
		
		nodeIndices = nI <- unique( c( 1, sort( which(oL[,1] %in% nodeHeights) ) ) )
		thetas <- sapply( 1:(length(nI)-1), function(i) exp(-(oL[ nI[i+1],2]-oL[nI[i],2]) ) )
		corates <- sapply( nodeHeights, function(h) dLambda( h, NA, parms )[[1]] )
		ll <- unname( sum(log(corates)) + sum(log(thetas))  )
	return(ll)
}

make.fgy <- function(t0, t1, births, deaths, nonDemeDynamics,  x0,  migrations=NA,  parms=NA, fgyResolution = 2000, integrationMethod = 'rk4')
{
#~ 	generates a discrete numeric representation of the demographic process given ODE as type str
	demeNames <- rownames(births)
	m <- nrow(births)
	nonDemeNames <- names(nonDemeDynamics)
	mm <- length(nonDemeNames)
	
	#reorder x0
	if (length(x0)!=m + mm) stop('initial conditons incorrect dimension', x0, m, mm) 
	if ( sum( !(c(demeNames, nonDemeNames) %in% names(x0)) )  > 0)  stop('initial conditions vector incorrect names', names(x0), demeNames, nonDemeNames)
	y0 <- x0[c(demeNames, nonDemeNames)]
	
	#parse equations
	pbirths <- sapply( 1:m, function(k) 
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
	#pndd <- ifelse(mm>0, sapply(1:mm, function(k) parse(text=nonDemeDynamics[k]) ), NA )
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
	#print(system.time(
		#ox <- ode(y=y0, times, func=dx, parms, method=integrationMethod)
	#))
	ox <- ode(y=y0, times, func=dx, parms, method=integrationMethod, verbose=FALSE)
	# note does not include first value, which is t0; 2nd value corresponds to root of tree
	Ys <- lapply( nrow(ox):(1+length(times0)), function(i) ox[i, demeNames] )
	Fs <- lapply( nrow(ox):(1+length(times0)), function(i) .birth.matrix(ox[i,], ox[i,1])  ) # dh *
	Gs <- lapply( nrow(ox):(1+length(times0)), function(i) .migration.matrix(ox[i,], ox[i,1])  ) #dh *
	
	list( rev(times), Fs, Gs, Ys , ox )
}

coalescent.log.likelihood <- function( bdt, births, deaths, nonDemeDynamics,  t0, x0,  migrations=NA,  parms=NA, fgyResolution = 2000, integrationMethod = 'rk4',  censorAtHeight=FALSE, forgiveAgtY=.2, returnTree=FALSE)
{
	if (is.vector( births)) return(coalescent.log.likelihood.unstructuredModel(
	   bdt, births,  deaths, nonDemeDynamics,  t0, x0, parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod,  censorAtHeight=censorAtHeight, forgiveAgtY=forgiveAgtY) )
	if (is.na(t0)) t0 <- bdt$maxSampleTime - bdt$maxHeight
	if ( (bdt$maxSampleTime - bdt$maxHeight) < t0)  {warning('Root of tree occurs before earliest time on time axis'); return(-Inf) }
	demeNames <- rownames(births)
	m <- nrow(births)
	nonDemeNames <- names(nonDemeDynamics)
	mm <- length(nonDemeNames)
	
	
	t1 <- bdt$maxSampleTime
	tfgy <- make.fgy( t0, t1, births, deaths, nonDemeDynamics,   x0,  migrations=migrations,  parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod )
	
	# make fgy parms for discretized rate functions
	fgyParms <- list()
	fgyParms$FGY_RESOLUTION		<- fgyResolution
	fgyParms$maxHeight				<- bdt$maxHeight #
	fgyParms$FGY_H_BOUNDARIES 		<- seq(0, bdt$maxHeight, length.out = fgyResolution) 
	# reverse order (present to past): 
	fgyParms$F_DISCRETE 			<- tfgy[[2]] #Fs
	fgyParms$G_DISCRETE 			<- tfgy[[3]] #Gs
	fgyParms$Y_DISCRETE 			<- tfgy[[4]] #Ys
	fgyParms$hoffset = hoffset <- 0
	fgyParms$get.index <- function(h) min(1 + fgyParms$FGY_RESOLUTION * (h + hoffset) / bdt$maxHeight, fgyParms$FGY_RESOLUTION )
	fgyParms$F. 					<- function(h) { fgyParms$F_DISCRETE[[fgyParms$get.index(h)]] } #
	fgyParms$G. 					<- function(h) { fgyParms$G_DISCRETE[[fgyParms$get.index(h)]] }
	fgyParms$Y. 					<- function(h) { fgyParms$Y_DISCRETE[[fgyParms$get.index(h)]] }
	
#~ 	fgyparms <- list()
#~ 	fgyparms$FGY_RESOLUTION		<- fgyResolution
#~ 	fgyparms$maxHeight			<- bdt$maxHeight #
#~ 	fgyparms$FGY_H_BOUNDARIES 		<- seq(0, bdt$maxHeight, length.out = fgyResolution) 
#~ 	fgyparms$F_DISCRETE 			<- tfgy[[2]] #Fs
#~ 	fgyparms$G_DISCRETE 			<- tfgy[[3]] #Gs
#~ 	fgyparms$Y_DISCRETE 			<- tfgy[[4]] #Ys
#~ 	fgyparms$F. 					<- function(FGY_INDEX) { fgyparms$F_DISCRETE[[FGY_INDEX]] } #note does not actually use arg t
#~ 	fgyparms$G. 					<- function(FGY_INDEX) { fgyparms$G_DISCRETE[[FGY_INDEX]] }
#~ 	fgyparms$Y. 					<- function(FGY_INDEX) { fgyparms$Y_DISCRETE[[FGY_INDEX]] }
	
	tree <- .calculate.internal.states(bdt, fgyParms,  censorAtHeight=censorAtHeight, forgiveAgtY = forgiveAgtY,   INTEGRATIONMETHOD = integrationMethod)
	
	i<- (length(tree$sampleTimes)+1):(tree$Nnode + length(tree$tip.label))
	if (censorAtHeight) { 
		internalHeights <- tree$heights[(length(tree$tip.label)+1):length(tree$heights)]
		i <- i[internalHeights <= censorAtHeight] 
	}
	ll <- sum( log(tree$coalescentRates[i]) ) + sum( tree$logCoalescentSurvivalProbability[i] )
	if (censorAtHeight) { ll<- tree$Nnode *  ll/length(i)}
	if (is.nan(ll) | is.na(ll) ) ll <- -Inf
	#return(list( ll, tree)) #for debugging
	if (returnTree) return(list(ll, tree))
	ll
}

solve.model <- function(t0,t1, x0, births,  deaths, nonDemeDynamics, parms, migrations=NA, integrationMethod = 'rk4', timeResolution=1000)
{
	if (is.vector(births)) return(solve.model.unstructured(t0,t1, x0, births,  deaths, nonDemeDynamics, parms, timeResolution=timeResolution, integrationMethod = integrationMethod) )
	
	tfgy <- make.fgy(t0, t1, births, deaths, nonDemeDynamics,  x0,  migrations=migrations,  parms=parms, fgyResolution = timeResolution, integrationMethod = integrationMethod )
	
	tfgy[[5]]
}


#' Simulate a binary dated tree from given model
#' @export
simulate.binary.dated.tree <- function(births, deaths, nonDemeDynamics,  t0, x0, sampleTimes, sampleStates,  migrations=NA,  parms=NA, fgyResolution = 2000, integrationMethod = 'rk4')
{
	m <- nrow(births)
	if (is.null(m))  stop('Error: currently only models with at least two demes are supported')
	if (m < 2)  stop('Error: currently only models with at least two demes are supported')
	maxSampleTime <- max(sampleTimes)
	tfgy <- make.fgy( t0, maxSampleTime, births, deaths, nonDemeDynamics,  x0,  migrations=migrations,  parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod )
	#~ 	simulate.binary.dated.tree.fgy(tfgy, sampleTimes, sampleStates,  parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod)
	simulate.binary.dated.tree.fgy( tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, integrationMethod = integrationMethod )
}

simulate.binary.dated.tree.unstructured <- function(births, deaths, nonDemeDynamics,  t0, x0, sampleTimes,  parms=NA, fgyResolution = 2000, integrationMethod = 'rk4')
{
	if (!is.vector(births)) stop('Error: use simulate.binary.dated.tree if there is more than one deme')
	if (length(births)>1) stop('Error: use simulate.binary.dated.tree if there is more than one deme')
	m <- 1
	if (length(names(sampleTimes))==0){
		sampleNames <- as.character(1:length(sampleTimes))
		names(sampleTimes) <- sampleNames
	}
	births2 <- matrix('0', nrow=2, ncol=2)
	births2[1,1] <- births
	DEMENAMES <- c( names(births), paste(sep='_', names(births), '0') )
	colnames(births2) = rownames(births2) <- DEMENAMES 
	deaths2 <- c(deaths, '0')
	names(deaths2) <- DEMENAMES
	sampleStates <- cbind( rep(1, length(sampleTimes)), rep(0, length(sampleTimes)) )
	colnames(sampleStates) <- DEMENAMES
	rownames(sampleStates) <- names(sampleTimes)
	x02 <- c(x0, 1)
	names(x02) <- c(names(x0), DEMENAMES[2])
	maxSampleTime <- max(sampleTimes)
	tfgy <- make.fgy( t0, maxSampleTime, births2, deaths2, nonDemeDynamics,  x02,   parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod )
	simulate.binary.dated.tree.fgy( tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, integrationMethod = integrationMethod )
}

simulate.binary.dated.tree.fgy <- function( times, births, migrations, demeSizes, sampleTimes, sampleStates, integrationMethod = 'rk4', n.reps=1, n.cores=1)
{
	if (n.cores > 1) {
		require(parallel, quietly=TRUE)
		mc <- getOption('mc.cores', n.cores)
	}

#NOTE assumes times in equal increments
#~ TODO mstates, ustates not in returned tree 
s <- round(coef( lm(x ~ y,  data.frame(x = 1:length(times), y = sort(times))) )[1], digits=9)
if (s!=1) warning('Tree simulator assumes times given in equal increments')
	n <- length(sampleTimes)
	m <- ncol(sampleStates)
	maxSampleTime <- max(sampleTimes)
	if (length(names(sampleTimes))==0){
		sampleNames <- as.character(1:length(sampleTimes))
		names(sampleTimes) <- sampleNames
		rownames(sampleStates) <- sampleNames
	}
	if (length(rownames(sampleStates))==0) warning('simulate.binaryDatedTree.fgy: sampleStates should have row names')
	
	sampleHeights <- maxSampleTime - sampleTimes 
	ix <- sort(sampleHeights, index.return=TRUE)$ix
	sortedSampleHeights <- sampleHeights[ix]
	sortedSampleStates <- as.matrix(sampleStates[ix,]) # if only one deme, this would otherwise return a vector
	uniqueSortedSampleHeights <- unique(sortedSampleHeights)
	
	maxtime <- max(times)
	mintime <- min(times)
	maxHeight <- maxSampleTime -  mintime
	
	times_ix <- sort(times, index.return=TRUE, decreasing=TRUE)$ix
	
	fgyParms <- list()
	fgyParms$FGY_RESOLUTION		<- length(times)
	# reverse order (present to past): 
	fgyParms$F_DISCRETE 			<- lapply(times_ix, function(k)  births[[k]] )
	fgyParms$G_DISCRETE 			<- lapply(times_ix, function(k)  migrations[[k]] )
	fgyParms$Y_DISCRETE 			<- lapply(times_ix, function(k)  demeSizes[[k]] )
	# the following line accounts for any discrepancies between the maximum time axis and the last sample time
	fgyParms$hoffset = hoffset <- maxtime - maxSampleTime
	if (hoffset < 0) stop( 'Time axis does not cover the last sample time' )
	fgyParms$get.index <- function(h) min(1 + fgyParms$FGY_RESOLUTION * (h + hoffset) / ( maxtime - mintime ), fgyParms$FGY_RESOLUTION )
	fgyParms$F. 					<- function(h) { fgyParms$F_DISCRETE[[fgyParms$get.index(h)]] } #
	fgyParms$G. 					<- function(h) { fgyParms$G_DISCRETE[[fgyParms$get.index(h)]] }
	fgyParms$Y. 					<- function(h) { fgyParms$Y_DISCRETE[[fgyParms$get.index(h)]] }
	fgyParms$FGY_H_BOUNDARIES 		<- sort( maxSampleTime - times )

	F. <- fgyParms$F.; G. <- fgyParms$G.; Y. <- fgyParms$Y. ; 
	FGY_RESOLUTION <- fgyParms$FGY_RESOLUTION
	
	get.fgy <- function(h)
	{
			.Y		<- fgyParms$Y.(h)
			.F		<- fgyParms$F.(h)
			.G		<- fgyParms$G.(h)
		list(.F=.F, .G=.G, .Y=.Y)
	}
	
	# construct forcing timeseries for ode's
	heights <- sort( fgyParms$FGY_H_BOUNDARIES )
	heights <- heights[heights <= maxHeight]
	heights <- heights[heights >= 0]
	fgyParms$heights <- heights
	fgymat <- t( sapply( fgyParms$heights, function(h) 
	  with(get.fgy(h), {
		c( as.vector(.F), as.vector(.G), as.vector(.Y) )
	  })
	) )
	fgymat <- pmax(fgymat, 0)
	
	.solve.Q.A.L <- function(h0, h1, A0, L0)
	{ # uses C implementation
		Q0 <- diag(m)
		parameters 	<- c(m, maxHeight, length(fgyParms$heights), sum(A0), as.vector(fgymat))
		y0 <- c( as.vector( Q0), A0,  L0 ) #
		o <- ode(y=y0, c(h0,h1), func = "dQAL", parms = parameters, dllname = "rcolgem", initfunc = "initfunc", method=integrationMethod )
		Q1 		<- t( matrix(  abs(o[nrow(o),2:(1 + m^2)]) , nrow=m) ) #NOTE the transpose
		A1 <- o[nrow(o), (1 + m^2 + 1):(1 + m^2 +  m)]
		L1 <- o[nrow(o), ncol(o)]
		return ( list(  unname(Q1), unname(A1),  unname(L1) ) ) #
	}
	
	run1 <- function(repl) {
		# solve for A
		cumSortedSampleStates <-  sapply( 1:m, function(k) cumsum(sortedSampleStates[,k]) )
		cumSortedNotSampledStates <-  t(cumSortedSampleStates[n,] - t(cumSortedSampleStates) )
		#nsy.index <- approxfun( sortedSampleHeights , 1:n , method='constant', rule=2)
		nsy.index <-  approxfun(
		   sort( jitter( sortedSampleHeights
		     , factor=max(1e-6, sortedSampleHeights[length(sortedSampleHeights)]/1e6)) )
		   , 1:n , method='constant', rule=2)
		not.sampled.yet <- function(h)
		{
			cumSortedNotSampledStates[ nsy.index(h), ]
		}
		dA <- function(h, A, parms, ...)
		{
			nsy <- not.sampled.yet(h)
			with(get.fgy(h),
			{
				A_Y 	<- (A-nsy) / .Y
				A_Y[is.nan(A_Y)] <- 0
				csFpG 	<- colSums( .F + .G )
				list( setNames( as.vector(
				  .G %*% A_Y - csFpG * A_Y + (.F %*% A_Y) * pmax(1-A_Y, 0)
				  ), names(A)
				))
			})
		}
		AIntervals <- c(sortedSampleHeights[2:length(uniqueSortedSampleHeights)], maxHeight)
		h0 						<- 0
		sampled.at.h <- function(h) which(sortedSampleHeights==h)

		haxis <- seq(0, maxHeight, length.out=fgyParms$FGY_RESOLUTION)
		AplusNotSampled <- ode( y = colSums(sortedSampleStates), times = haxis, func=dA, parms=NA, method = integrationMethod)[, 2:(m+1)]
		AplusNotSampled <- as.matrix(AplusNotSampled, nrows=length(haxis))
		Amono <- rowSums( AplusNotSampled )
		Amono[is.na(Amono)] <- min(Amono[!is.na(Amono)] )
		Amono <- Amono - min(Amono)
		Amono <- (max(Amono)  - Amono) / max(Amono)
		nodeHeights <- sort( approx( Amono, haxis, xout=runif(n-1, 0, 1) )$y ) # careful of impossible node heights
		uniqueSampleHeights <- unique( sampleHeights )
		eventTimes <- c(uniqueSampleHeights, nodeHeights)
		isSampleEvent <- c(rep(TRUE, length(uniqueSampleHeights)), rep(FALSE, length(nodeHeights)))
		ix <- sort(eventTimes, index.return=TRUE)$ix
		eventTimes <- eventTimes[ix]
		isSampleEvent <- isSampleEvent[ix]

		get.A.index <- function(h){
			min(1 + floor(fgyParms$FGY_RESOLUTION * (h ) / ( maxHeight)), fgyParms$FGY_RESOLUTION )
		}
		get.A <- function(h){
			i <- get.A.index(h)
			AplusNotSampled[i,] - not.sampled.yet(h)
		}

		S <- 1
		L <- 0
		# initialize variables; tips in order of sortedSampleHeights
		Nnode <- n-1
		edge.length 	<- rep(-1, Nnode + n-1) # should not have root edge
		edge			<- matrix(-1, (Nnode + n-1), 2)
		if (length(names(sortedSampleHeights))==0)
		{
			tip.label 		<- as.character(1:n)
		} else{
			tip.label 		<- names(sortedSampleHeights)# as.character( 1:n )
		}
		maxSampleTime 	<- max(sampleTimes)
		heights 		<- rep(0, (Nnode + n) )
		parentheights 	<- rep(-1, (Nnode + n) )
		heights[1:n] 	<- sortedSampleHeights
		inEdgeMap 		<- rep(-1, Nnode + n)
		outEdgeMap 		<- matrix(-1, (Nnode + n), 2)
		parent 			<- 1:(Nnode + n)
		daughters 		<- matrix(-1, (Nnode + n), 2)
		lstates 		<- matrix(-1, (Nnode + n), m)
		mstates 		<- matrix(-1, (Nnode + n), m)
		ustates 		<- matrix(-1, (Nnode + n), m)
		ssm 			<- matrix( 0, nrow=n, ncol=m)
		lstates[1:n,]	<-  sortedSampleStates
		mstates[1:n,] 	<- lstates[1:n,]

		isExtant 		<- rep(FALSE, Nnode+n)
		isExtant[sampled.at.h(h0)] <- TRUE
		extantLines <- which( isExtant)
		nExtant <- sum(isExtant)

		#print(paste(date(), 'simulate tree'))

		if (length(extantLines) > 1){
			A0 <- colSums(as.matrix(sortedSampleStates[extantLines,], nrow=length(extantLines)) )
		} else{
			A0 <- sortedSampleStates[extantLines,]
		}
		lineageCounter <- n+1

		for (ih in 1:(length(eventTimes)-1)) {
			h0 <- eventTimes[ih]
			h1 <- eventTimes[ih+1]
			fgy <- get.fgy(h1)

			nExtant <- sum(isExtant)


			#get A0, process new samples, calculate state of new lines
			A0 <- get.A(h0) #A[get.A.index(h0),]
			out <- .solve.Q.A.L(h0, h1, A0,  L)
			Q <- out[[1]]
			A <- out[[2]]
			L <- out[[3]]

			# clean output
			if (is.nan(L)) {L <- Inf}
			if (sum(is.nan(Q)) > 0) { return(NA) }#browser() #Q <- diag(length(A))
			if (sum(is.nan(A)) > 0) A <- A0

			#update mstates
			if ( nExtant > 1)
			{
				mstates[isExtant,] <- t( t(Q) %*% t(mstates[isExtant,])  )
				mstates[isExtant,] <- abs(mstates[isExtant,]) / rowSums(as.matrix(abs(mstates[isExtant,]), nrow=length(isExtant)))
				#recalculate A
				A <- colSums(as.matrix(mstates[isExtant,], nrow=length(isExtant)))
			}
			else{
				mstates[isExtant,] <- t( t(Q) %*% mstates[isExtant,] )
				mstates[isExtant,] <- abs(mstates[isExtant,]) / sum(abs(mstates[isExtant,]))
				#recalculate A
				A <- mstates[isExtant,]
			}

			if (isSampleEvent[ih+1])
			{
				sat_h1 <- sampled.at.h(h1)
				isExtant[sat_h1] <- TRUE
				#extantLines <- c(extantLines, sat_h1)
				heights[sat_h1] <- h1
			} else{ #coalecent event
				#with(fgy, {
				#attach(fgy)
				.F <- fgy$.F
				.G <- fgy$.G
				.Y <- fgy$.Y

				if (nExtant>1 && is.element(0, .Y)) {
				    # last two lineages have not coalesced by simulation time zero
				    # reject this tree and return NA to prevent divide-by-zero error
				    return(NA)
				}

					a <- A / .Y

		#nextant <- length(extantLines)
		extantLines <- which(isExtant)


	#~ 			if (F){
	#~ 				astates <- matrix(  pmin(1, t(t(mstates[isExtant,])/.Y ) ),  nrow=nExtant  )
	#~ 				tryCatch(
	#~ 					{ .lambdamat <- astates %*% .F %*% t(astates) }
	#~ 				 , error = function(e) browser())
	#~ 				diag(.lambdamat) <- 0
	#~ 				uivi <- sample.int( nExtant^2, size=1, prob=as.vector(.lambdamat) )
	#~ 				u_i <- 1 + ((uivi-1) %% nExtant)#row
	#~ 				v_i <- 1 + floor( (uivi-1) / nExtant ) #column
	#~ 				u <-  extantLines[ u_i ]
	#~ 				v <- extantLines[ v_i ]
	#~ 			}
	tryCatch({
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
	}, error = function(e) browser() )
	#~ tryCatch( {
					ustates[u,] <- mstates[u,]
					ustates[v,] <- mstates[v,]
	#~ }, error = function(e) browser() )
					#lambda_uv  <- as.vector( astates[u_i,] %*% .F %*% astates[v_i,] ) + as.vector( astates[v_i,] %*% .F %*% astates[u_i,] )
					#lambda_uv <- ((astates[u_i,]) %*% t( astates[v_i,] )) * .F + ((astates[v_i,]) %*% t( astates[u_i,] )) * .F

					a_u <- pmin(1, mstates[u,] / .Y )
					a_v <- pmin(1, mstates[v,] / .Y )
					lambda_uv <- ((a_u) %*% t( a_v )) * .F + ((a_v) %*% t( a_u )) * .F

					palpha <- rowSums(lambda_uv) / sum(lambda_uv)
	#~ browser()
					alpha <- lineageCounter
					lineageCounter <- lineageCounter + 1
					isExtant[alpha] <- TRUE
					isExtant[u] <- FALSE
					isExtant[v] <- FALSE
		#extantLines <- c(extantLines[-c(v_i, u_i)], alpha)
					lstates[alpha,] = mstates[alpha,] <- palpha
					heights[alpha] <- h1

					uv <- c(u, v)
					inEdgeMap[uv] <- alpha#lineageCounter
					outEdgeMap[alpha,] <- uv
					parent[uv] <- alpha
					parentheights[uv] <- h1
					daughters[alpha,] <- uv
					edge[u,] <- c(alpha, u) #
					edge.length[u] <- h1 - heights[u]
					edge[v,] <- c(alpha,v) #
					edge.length[v] <- h1 - heights[v]
			} # end else
		} # end for loop

        tryCatch ({
    		self<- list(edge=edge, edge.length=edge.length, Nnode=Nnode, tip.label=tip.label, heights=heights, parentheights=parentheights, parent=parent, daughters=daughters, lstates=lstates, mstates=mstates, ustates=ustates, m = m, sampleTimes = sampleTimes, sampleStates= sampleStates, maxSampleTime=maxSampleTime, inEdgeMap = inEdgeMap, outEdgeMap=outEdgeMap)
	    	class(self) <- c("binaryDatedTree", "phylo")
	    	#~ <reorder edges for compatibility with ape::phylo functions>
	    	#~ (ideally ape would not care about the edge order, but actually most functions assume a certain order)
	    	sampleTimes2 <- sampleTimes[names(sortedSampleHeights)]
	    	sampleStates2 <- as.matrix(lstates[1:n,], nrow=n); rownames(sampleStates2) <- tip.label
	        #~ browser()
		    phylo <- read.tree(text=write.tree(self) )

    		sampleTimes2 <- sampleTimes2[phylo$tip.label];
    		sampleStates2 <- as.matrix(sampleStates2[phylo$tip.label,], nrow=length(phylo$tip.label))
    		bdt <- binaryDatedTree(phylo, sampleTimes2, sampleStates = sampleStates2)
		}, error = function(e) browser())
		return(bdt)
	} # end run1()

	if (n.cores > 1) {
		result <- mclapply(1:n.reps, run1)
	} else {
		result <- lapply(1:n.reps, run1)
	}

    # exclude NA values
	return (result[!is.na(result)])
}



#' Simulate a binary dated tree from given model
#' @export
#~ TODO : Error in kids[[parent[i]]] : attempt to select more than one element
simulate.binary.dated.tree.0 <- function(births, deaths, nonDemeDynamics,  t0, x0, sampleTimes, sampleStates,  migrations=NA,  parms=NA, fgyResolution = 2000, integrationMethod = 'rk4')
{
	require(ape)
#~ same attributes are defined as binaryDatedTree
#~ <preliminaries>
	n 			<- length(sampleTimes)
	Nnode 		<-  n-1
	
	FGY_INDEX <- 2
	
	demeNames <- rownames(births)
	m <- nrow(births)
	if (m < 2)  stop('Error: currently only models with at least two demes are supported')
	nonDemeNames <- names(nonDemeDynamics)
	mm <- length(nonDemeNames)
	sortedSampleTimes <- sort(sampleTimes)
	maxSampleTime <- max(sampleTimes)
	sampleHeights <- maxSampleTime- sampleTimes 
	sortedSampleHeights <- sort(sampleHeights)
	uniqueSortedSampleHeights <- unique(sortedSampleHeights)
	ussh_index <- 2
	
	tfgy <- make.fgy( t0, maxSampleTime, births, deaths, nonDemeDynamics,  x0,  migrations=migrations,  parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod )
	
	# make fgy parms for discretized rate functions
	globalParms <- list()
	globalParms$FGY_RESOLUTION			<- fgyResolution
	globalParms$maxHeight				<- maxSampleTime-t0 #
	globalParms$FGY_H_BOUNDARIES 		<- tfgy[[1]] #seq(0, maxSampleTime-t0, length.out = fgyResolution) 
	globalParms$F_DISCRETE 			<- tfgy[[2]] #Fs
	globalParms$G_DISCRETE 			<- tfgy[[3]] #Gs
	globalParms$Y_DISCRETE 			<- tfgy[[4]] #Ys
	globalParms$F. 					<- function(FGY_INDEX) { globalParms$F_DISCRETE[[FGY_INDEX]] } #note does not actually use arg t
	globalParms$G. 					<- function(FGY_INDEX) { globalParms$G_DISCRETE[[FGY_INDEX]] }
	globalParms$Y. 					<- function(FGY_INDEX) { globalParms$Y_DISCRETE[[FGY_INDEX]] }
	
	#edge <- c()
	#edge.length <- c()
	edge.length 	<- rep(-1, Nnode + n-1) # should not have root edge
	edge			<- matrix(-1, (Nnode + n-1), 2)
	tip.label 		<- as.character( 1:n )
	maxSampleTime  	= T <- max(sampleTimes)
	heights 		<- rep(0, (Nnode + n) )
	parentheights 	<- rep(-1, (Nnode + n) )
	heights[1:n] 	<- maxSampleTime - sampleTimes
	inEdgeMap 		<- rep(-1, Nnode + n)
	outEdgeMap 		<- matrix(-1, (Nnode + n), 2)
	parent 			<- 1:(Nnode + n) 
	daughters 		<- matrix(-1, (Nnode + n), 2)
	lstates 		<- matrix(-1, (Nnode + n), m)
	mstates 		<- matrix(-1, (Nnode + n), m)
	ustates 		<- matrix(-1, (Nnode + n), m)
	ssm <- matrix( 0, nrow=n, ncol=m)
	lstates[1:n,] <- t( sapply( 1:n, function(i) diag(m)[sampleStates[i],] ) )
	mstates[1:n,] 	<- lstates[1:n,]

#~ </preliminaries>
	
#~ <survival time to next event>
	calculate.rates <- function(h, parms, globalParms)
	{
		# eliminate diag elements for migration
		t 					<- parms$maxSampleTime - h
			.G 					<- globalParms$G.(FGY_INDEX) 
			.F 					<- globalParms$F.(FGY_INDEX)
			.Y 					<- globalParms$Y.(FGY_INDEX) 
		X1 <- parms$A / .Y
		X2 					<-  (.Y - parms$A ) / .Y
		X1[is.nan(X1)] 		<- 0
		X2[is.nan(X2)] 		<- 0
		X1[is.infinite(X1)] <- 0
		X2[is.infinite(X2)] <- 0
		X1mat 				<- matrix(X1,nrow=parms$m,ncol=parms$m, byrow=TRUE)
		tX1mat 				<- matrix(X1,nrow=parms$m,ncol=parms$m, byrow=FALSE)
		tX2mat 				<- matrix(X2, nrow=parms$m,ncol=parms$m, byrow=FALSE)
		lambdaCoalescent 	<-  .F * X1mat * tX1mat
		lambdaMigration 	<-  t(.G * X1mat  )
		lambdaInvisibleTransmission	<-  t(.F * tX2mat * X1mat)
		diag(lambdaMigration) 		<- 0
		diag(lambdaInvisibleTransmission) <- 0
		#browser()
		return( list( X1=X1, X2=X2, lambdaCoalescent=lambdaCoalescent, lambdaMigration=lambdaMigration, lambdaInvisibleTransmission=lambdaInvisibleTransmission) )
	}
	dtheta <- function(h, theta, parms, globalParms, ...){
		if (is.nan(theta)) return(list(NaN))
		if (theta <= parms$r) return(list(NaN)) #terminate early if survival time found
		rates <- calculate.rates(h, parms, globalParms)
		return( list( -theta * (sum(rates$lambdaCoalescent) + sum(rates$lambdaMigration) + sum(rates$lambdaInvisibleTransmission)) ) )
	}
#~ </survival time to next event>

.calculate.A <- function(mstates, extant)
{
	if (length(extant) > 1) return( colSums(mstates[extant,]) )
	mstates[extant,]
}
	
#~ <simulate tree>
#TODO these trees are still biased? check rate matrices
#~ 	extant <- 1:n
	extant <- (1:n)[sampleHeights==0]
	lineageCounter <- n+1 # next lineage will have this index
	A <- .calculate.A(mstates, extant)
	h <- 0
	notdone <- TRUE
	lastExtant <- lineageCounter-1 #DEBUG
	nextSampleHeight <-  ifelse( ussh_index > length(uniqueSortedSampleHeights), Inf, uniqueSortedSampleHeights[ussh_index])
	while(notdone){
		r <- runif(1)
		theta0 <- 1
		parms <- list(A=A, m=m, maxSampleTime=maxSampleTime, r=r)
		
		# find appropriate time axis
		rates <- calculate.rates( h, parms, globalParms) 
		parms$rates <- rates
		lambda <- (sum(rates$lambdaCoalescent) + sum(rates$lambdaMigration) + sum(rates$lambdaInvisibleTransmission))
		
		# solve survival time
		# NOTE more accurate to interpolate approxfun( thetaTimes, o[,2] ), & invert to find o[o[,2]==r,1]
		
			eventTime <- rexp(1, rate=lambda) + h
			nextBoundaryTime <- min( nextSampleHeight, globalParms$FGY_H_BOUNDARIES[FGY_INDEX] )
			while (FGY_INDEX < globalParms$FGY_RESOLUTION &  eventTime > nextBoundaryTime) {
					eventTime <- nextBoundaryTime
					if (eventTime == nextSampleHeight)
					{
						extant <- c(extant, (1:n)[sampleHeights==eventTime])
						A <- .calculate.A(mstates, extant)
						parms$A <- A
						tryCatch({
								rates <- calculate.rates( eventTime, parms, globalParms) ; parms$rates <- rates
							}, error = function(e) browser() )
						lambda <- (sum(rates$lambdaCoalescent) + sum(rates$lambdaMigration) + sum(rates$lambdaInvisibleTransmission))
						eventTime <- eventTime + rexp(1, rate=lambda)
						ussh_index <- ussh_index + 1
						nextSampleHeight <-  ifelse( ussh_index > length(uniqueSortedSampleHeights), Inf, uniqueSortedSampleHeights[ussh_index])
						nextBoundaryTime <- min( nextSampleHeight, globalParms$FGY_H_BOUNDARIES[FGY_INDEX] )
					}
					else{
						FGY_INDEX <- FGY_INDEX+1 ; 
						tryCatch({
									rates <- calculate.rates( eventTime, parms, globalParms) ; parms$rates <- rates
								}, error = function(e) browser() )
						lambda <- (sum(rates$lambdaCoalescent) + sum(rates$lambdaMigration) + sum(rates$lambdaInvisibleTransmission))
						eventTime <- eventTime + rexp(1, rate=lambda)
						nextBoundaryTime <- min( nextSampleHeight, globalParms$FGY_H_BOUNDARIES[FGY_INDEX] )
					}
			}
		 
		# which event? 2m2 possibilities
		cumulativeRatesVector <- cumsum( c( as.vector(rates$lambdaCoalescent), as.vector(rates$lambdaMigration+rates$lambdaInvisibleTransmission)) )
		# as.vector unravels matrices by column
		cumulativeRatesVector <- cumulativeRatesVector / cumulativeRatesVector[length(cumulativeRatesVector) ]
		rwhichevent <- runif(1)
		indexEvent <- match(TRUE, cumulativeRatesVector > rwhichevent)
		if (is.na(indexEvent)){ 
			warning('Rates at ', eventTime, 'are NA or all rates are zero')
			indexEvent <- round( runif(1) * length(cumulativeRatesVector) )
		}
		#cumulativeRatesVector
		#indexEvent
		
		if (indexEvent <= m*m){ #coalescent event
			l <- 1+floor( (indexEvent-1) / m )
			k <-  1+ ( (indexEvent-1) %% m) #1+ ((indexEvent+1) %% m)
			kextant <- extant[ mstates[extant, k]==1 ]
			lextant <- extant[ mstates[extant, l]==1 ]
			if (k==l){ # guard against A[k]<2
				if (length(kextant) >=2) {ij <- sample(kextant, 2) }
				else { ij <- sample(extant, 2) }
			} else{
				# horribly ugly code thanks to stupid default behavior of 'sample':
				if ((length(kextant)>1) & (length(lextant)>1)) { ij <- c( sample(kextant, 1), sample(lextant, 1) ) }
				else if (length(kextant)==0 | length(lextant)==0) { warning('kextant length zero'); ij <- sample(extant, 2) }
				else if (length(kextant)==1 & length(lextant)==1) { ij <- c(kextant, lextant) }
				else if (length(kextant)==1) { ij <- c(kextant, sample(lextant, 1)) }
				else if (length(lextant)==1) { ij <- c(sample(kextant, 1), lextant) }
				else { browser() }
			}
			# set new lineage
			extant <- c(extant, lineageCounter)
			extant <- extant[extant!=ij[1] & extant!=ij[2]]
			heights[lineageCounter] <- eventTime
			lstates[lineageCounter,] <- diag(m)[k,]
			mstates[lineageCounter,] <- diag(m)[k,]
			inEdgeMap[ij] <- lineageCounter
			outEdgeMap[lineageCounter,] <- ij
			parent[ij] <- lineageCounter
			parentheights[ij] <- eventTime
			daughters[lineageCounter,] <- ij
			#edge <- rbind(edge, c(lineageCounter, ij[1]))
			#edge.length <- c(edge.length, eventTime - heights[ij[1]] )
			#edge <- rbind(edge, c(lineageCounter, ij[2]))
			#edge.length <- c(edge.length, eventTime - heights[ij[2]] )
			edge[ij[1],] <- c(lineageCounter, ij[1]) # 
			edge.length[ij[1]] <- eventTime - heights[ij[1]] 
			edge[ij[2],] <- c(lineageCounter, ij[2]) # edge <- rbind(edge, c(lineageCounter, ij[2]))
			edge.length[ij[2]] <- eventTime - heights[ij[2]] 
			lineageCounter <- lineageCounter+1
			#print(paste(date(), 'coalescent', h, eventTime,  k, l, length(extant)))
			#if (length(extant)>= lastExtant) browser() #DEBUG
			#lastExtant <- length(extant) #DEBUG
			#if (rates$lambdaCoalescent[k,l]==0) {print('coalescent'); browser()}
		} else{ #migration event
			indexEvent <- indexEvent- m*m
			l <- 1+floor( (indexEvent-1) / m )
			k <- 1+ ( (indexEvent-1) %% m) #1+ ( (indexEvent+1) %% m)
			kextant <- extant[ mstates[extant, k]==1 ]
			if (length(kextant) > 1) { i <- sample(kextant, 1) }
			else { i <- sample(extant, 1) }
			mstates[i,] <- diag(m)[l,]
			#print(paste( date(),'migration', h, eventTime, k, l) )
			#print(paste(  k, l) )
			#migrations <- rates$lambdaMigration+rates$lambdaInvisibleTransmission
			#if (migrations[k,l]==0) {print('migrations'); browser()}
		}
		if(length(extant) <= 1 ) {notdone <- FALSE}
		else{
			A <- .calculate.A(mstates, extant)
			h <-  eventTime
		}
	}
#~ </simulate tree>
	#~ assemble class
	self<- list(edge=edge, edge.length=edge.length, Nnode=Nnode, tip.label=tip.label, heights=heights, parentheights=parentheights, parent=parent, daughters=daughters, lstates=lstates, mstates=mstates, ustates=ustates, m = m, sampleTimes = sampleTimes, sampleStates= sampleStates, maxSampleTime=maxSampleTime, inEdgeMap = inEdgeMap, outEdgeMap=outEdgeMap)
	class(self) <- c("binaryDatedTree", "phylo")
#~ </>
	
	
#~ <reorder edges for compatibility with ape::phylo functions> 
#~ (ideally ape would not care about the edge order, but actually most functions assume a certain order)
	sampleTimes2 <- sampleTimes; names(sampleTimes2) <- tip.label
	sampleStates2 <- lstates[1:n,]; rownames(sampleStates2) <- tip.label
	phylo <- read.tree(text=write.tree(self) )
	sampleTimes2 <- sampleTimes2[phylo$tip.label]; sampleTimes2 <- unname(sampleTimes2)
	sampleStates2 <- sampleStates2[phylo$tip.label,]; sampleStates2 <- unname(sampleStates2)
	names(sampleTimes2) <- phylo$tip.label
	rownames(sampleStates2) <- phylo$tip.label
	binaryDatedTree(phylo, sampleTimes2, sampleStates = sampleStates2)
#~ </>
}

simulate.binary.dated.tree.2 <- function(births, deaths, nonDemeDynamics,  t0, x0, sampleTimes, sampleStates,  migrations=NA,  parms=NA, fgyResolution = 2000, integrationMethod = 'rk4')
{
	m <- nrow(births)
	if (m < 2)  stop('Error: currently only models with at least two demes are supported')
	maxSampleTime <- max(sampleTimes)
	tfgy <- make.fgy( t0, maxSampleTime, births, deaths, nonDemeDynamics,  x0,  migrations=migrations,  parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod )
	#~ 	simulate.binary.dated.tree.fgy(tfgy, sampleTimes, sampleStates,  parms=parms, fgyResolution = fgyResolution, integrationMethod = integrationMethod)
	simulate.binary.dated.tree.fgy.2( tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]], sampleTimes, sampleStates, integrationMethod = integrationMethod )
}

simulate.binary.dated.tree.fgy.2 <- function( times, births, migrations, demeSizes, sampleTimes, sampleStates, integrationMethod = 'rk4')
{
#NOTE assumes times in equal increments
#~ TODO mstates, ustates not in returned tree 
s <- round(coef( lm(x ~ y,  data.frame(x = 1:length(times), y = sort(times))) )[1], digits=9)
if (s!=1) warning('Tree simulator assumes times given in equal increments')
	n <- length(sampleTimes)
	m <- ncol(sampleStates)
	maxSampleTime <- max(sampleTimes)
	if (length(names(sampleTimes))==0){
		sampleNames <- as.character(1:length(sampleTimes))
		names(sampleTimes) <- sampleNames
		rownames(sampleStates) <- sampleNames
	}
	if (length(rownames(sampleStates))==0) warning('simulate.binaryDatedTree.fgy: sampleStates should have row names')
	
	sampleHeights <- maxSampleTime - sampleTimes 
	ix <- sort(sampleHeights, index.return=TRUE)$ix
	sortedSampleHeights <- sampleHeights[ix]
	sortedSampleStates <- sampleStates[ix,]
	uniqueSortedSampleHeights <- unique(sortedSampleHeights)
	
	maxtime <- max(times)
	mintime <- min(times)
	maxHeight <- maxSampleTime -  mintime
	
	times_ix <- sort(times, index.return=TRUE, decreasing=TRUE)$ix
	
	fgyParms <- list()
	fgyParms$FGY_RESOLUTION		<- length(times)
	# reverse order (present to past): 
	fgyParms$F_DISCRETE 			<- lapply(times_ix, function(k)  births[[k]] )
	fgyParms$G_DISCRETE 			<- lapply(times_ix, function(k)  migrations[[k]] )
	fgyParms$Y_DISCRETE 			<- lapply(times_ix, function(k)  demeSizes[[k]] )
	# the following line accounts for any discrepancies between the maximum time axis and the last sample time
	fgyParms$hoffset = hoffset <- maxtime - maxSampleTime
	if (hoffset < 0) stop( 'Time axis does not cover the last sample time' )
	fgyParms$get.index <- function(h) min(1 + fgyParms$FGY_RESOLUTION * (h + hoffset) / ( maxtime - mintime ), fgyParms$FGY_RESOLUTION )
	fgyParms$F. 					<- function(h) { fgyParms$F_DISCRETE[[fgyParms$get.index(h)]] } #
	fgyParms$G. 					<- function(h) { fgyParms$G_DISCRETE[[fgyParms$get.index(h)]] }
	fgyParms$Y. 					<- function(h) { fgyParms$Y_DISCRETE[[fgyParms$get.index(h)]] }
	fgyParms$FGY_H_BOUNDARIES 		<- sort( maxSampleTime - times )

	F. <- fgyParms$F.; G. <- fgyParms$G.; Y. <- fgyParms$Y. ; 
	FGY_RESOLUTION <- fgyParms$FGY_RESOLUTION
	
	get.fgy <- function(h)
	{
			.Y		<- fgyParms$Y.(h)
			.F		<- fgyParms$F.(h)
			.G		<- fgyParms$G.(h)
		list(.F=.F, .G=.G, .Y=.Y)
	}
	
	# construct forcing timeseries for ode's
	heights <- sort( fgyParms$FGY_H_BOUNDARIES )
	heights <- heights[heights <= maxHeight]
	heights <- heights[heights >= 0]
	fgyParms$heights <- heights
	fgymat <- t( sapply( fgyParms$heights, function(h) 
	  with(get.fgy(h), {
		c( as.vector(.F), as.vector(.G), as.vector(.Y) )
	  })
	) )
	fgymat <- pmax(fgymat, 0)
	
	.solve.Q.A.L <- function(h0, h1, A0, L0)
	{ # uses C implementation
		Q0 <- diag(m)
		parameters 	<- c(m, maxHeight, length(fgyParms$heights), sum(A0), as.vector(fgymat))
		y0 <- c( as.vector( Q0), A0,  L0 ) #
		o <- ode(y=y0, c(h0,h1), func = "dQAL", parms = parameters, dllname = "rcolgem", initfunc = "initfunc", method=integrationMethod )
		Q1 		<- t( matrix(  abs(o[nrow(o),2:(1 + m^2)]) , nrow=m) ) #NOTE the transpose
		A1 <- o[nrow(o), (1 + m^2 + 1):(1 + m^2 +  m)]
		L1 <- o[nrow(o), ncol(o)]
		return ( list(  unname(Q1), unname(A1),  unname(L1) ) ) #
	}
	
	
	# solve for A
	cumSortedSampleStates <-  sapply( 1:m, function(k) cumsum(sortedSampleStates[,k]) )
	cumSortedNotSampledStates <-  t(cumSortedSampleStates[n,] - t(cumSortedSampleStates) )
	#nsy.index <- approxfun( sortedSampleHeights , 1:n , method='constant', rule=2)
	nsy.index <-  approxfun(
	   sort( jitter( sortedSampleHeights
	     , factor=max(1e-6, sortedSampleHeights[length(sortedSampleHeights)]/1e6)) )
	   , 1:n , method='constant', rule=2)
	not.sampled.yet <- function(h)
	{
		cumSortedNotSampledStates[ nsy.index(h), ]
	}
	dA <- function(h, A, parms, ...)
	{
		nsy <- not.sampled.yet(h) 
		with(get.fgy(h), 
		{ 
			A_Y 	<- (A-nsy) / .Y
			csFpG 	<- colSums( .F + .G )
			list( setNames( as.vector(
			  .G %*% A_Y - csFpG * A_Y + (.F %*% A_Y) * pmax(1-A_Y, 0) 
			  ), names(A)
			))
		})
	}
	AIntervals <- c(sortedSampleHeights[2:length(uniqueSortedSampleHeights)], maxHeight)
	h0 						<- 0
	sampled.at.h <- function(h) which(sortedSampleHeights==h)
	extantLines <- sampled.at.h(h0)
	haxis <- seq(0, maxHeight, length.out=fgyParms$FGY_RESOLUTION)
	AplusNotSampled <- ode( y = colSums(sortedSampleStates), times = haxis, func=dA, parms=NA, method = integrationMethod)[, 2:(m+1)]
	Amono <- rowSums( AplusNotSampled)
	Amono[is.na(Amono)] <- min(Amono[!is.na(Amono)] )
	Amono <- Amono - min(Amono)
	Amono <- (max(Amono)  - Amono) / max(Amono)
	nodeHeights <- sort( approx( Amono, haxis, xout=runif(n-1, 0, 1) )$y ) # careful of impossible node heights
	uniqueSampleHeights <- unique( sampleHeights )
	eventTimes <- c(uniqueSampleHeights, nodeHeights) 
	isSampleEvent <- c(rep(TRUE, length(uniqueSampleHeights)), rep(FALSE, length(nodeHeights)))
	ix <- sort(eventTimes, index.return=TRUE)$ix
	eventTimes <- eventTimes[ix]
	isSampleEvent <- isSampleEvent[ix]
	
	get.A.index <- function(h){
		min(1 + floor(fgyParms$FGY_RESOLUTION * (h ) / ( maxHeight)), fgyParms$FGY_RESOLUTION )
	}
	get.A <- function(h){
		i <- get.A.index(h)
		AplusNotSampled[i,] - not.sampled.yet(h) 
	}
	
	S <- 1
	L <- 0
	# initialize variables; tips in order of sortedSampleHeights
	Nnode <- n-1
	edge.length 	<- rep(-1, Nnode + n-1) # should not have root edge
	edge			<- matrix(-1, (Nnode + n-1), 2)
	if (length(names(sortedSampleHeights))==0)
	{
		tip.label 		<- as.character(1:n)
	} else{
		tip.label 		<- names(sortedSampleHeights)# as.character( 1:n )
	}
	maxSampleTime 	<- max(sampleTimes)
	heights 		<- rep(0, (Nnode + n) )
	parentheights 	<- rep(-1, (Nnode + n) )
	heights[1:n] 	<- sortedSampleHeights
	inEdgeMap 		<- rep(-1, Nnode + n)
	outEdgeMap 		<- matrix(-1, (Nnode + n), 2)
	parent 			<- 1:(Nnode + n) 
	daughters 		<- matrix(-1, (Nnode + n), 2)
	lstates 		<- matrix(-1, (Nnode + n), m)
	mstates 		<- matrix(-1, (Nnode + n), m)
	ustates 		<- matrix(-1, (Nnode + n), m)
	ssm 			<- matrix( 0, nrow=n, ncol=m)
	lstates[1:n,]	<-  sortedSampleStates
	mstates[1:n,] 	<- lstates[1:n,]
	
	print(paste(date(), 'simulate tree'))
	# simulate
	if (length(extantLines) > 1){
		A0 <- colSums(sortedSampleStates[extantLines,] ) 
	} else{
		A0 <- sortedSampleStates[extantLines,]
	}
	lineageCounter <- n+1
	for (ih in 1:(length(eventTimes)-1)){
		h0 <- eventTimes[ih]
		h1 <- eventTimes[ih+1]
		fgy <- get.fgy(h1)
		
		#get A0, process new samples, calculate state of new lines
		A0 <- get.A(h0) #A[get.A.index(h0),]
		out <- .solve.Q.A.L(h0, h1, A0,  L)
		Q <- out[[1]]
		A <- out[[2]]
		L <- out[[3]]
		
		# clean output
		if (is.nan(L)) {L <- Inf}
		if (sum(is.nan(Q)) > 0) browser() #Q <- diag(length(A))
		if (sum(is.nan(A)) > 0) A <- A0
		
		#update mstates
		if (length(extantLines) > 1)
		{
			mstates[extantLines,] <- t( t(Q) %*% t(mstates[extantLines,])  )
			mstates[extantLines,] <- abs(mstates[extantLines,]) / rowSums(abs(mstates[extantLines,]))
			#recalculate A
			A <- colSums(mstates[extantLines,])
		}
		else{
			mstates[extantLines,] <- t( t(Q) %*% mstates[extantLines,] )
			mstates[extantLines,] <- abs(mstates[extantLines,]) / sum(abs(mstates[extantLines,]))
			#recalculate A
			A <- mstates[extantLines,] 
		}
		
		if (isSampleEvent[ih+1])
		{
			sat_h1 <- sampled.at.h(h1)
			extantLines <- c(extantLines, sat_h1)
			heights[sat_h1] <- h1
		} else{ #coalecent event
			#with(fgy, {
			#attach(fgy)
			.F <- fgy$.F
			.G <- fgy$.G
			.Y <- fgy$.Y
				a <- A / .Y
				
				nextant <- length(extantLines)
				astates <- matrix(  pmin(1, t(t(mstates[extantLines,])/.Y ) ),  nrow=nextant  )
				tryCatch( 
					{ .lambdamat <- astates %*% .F %*% t(astates) }
				 , error = function(e) browser())
				diag(.lambdamat) <- 0
				uivi <- sample.int( nextant^2, size=1, prob=as.vector(.lambdamat) )
				u_i <- 1 + ((uivi-1) %% nextant)#row
				v_i <- 1 + floor( (uivi-1) / nextant ) #column
				u <-  extantLines[ u_i ] 
				v <- extantLines[ v_i ] 
				ustates[u,] <- mstates[u,]
				ustates[v,] <- mstates[v,]
				#lambda_uv  <- as.vector( astates[u_i,] %*% .F %*% astates[v_i,] ) + as.vector( astates[v_i,] %*% .F %*% astates[u_i,] )
				lambda_uv <- ((astates[u_i,]) %*% t( astates[v_i,] )) * .F + ((astates[v_i,]) %*% t( astates[u_i,] )) * .F 
				palpha <- rowSums(lambda_uv) / sum(lambda_uv)
#~ browser()
				alpha <- lineageCounter 
				lineageCounter <- lineageCounter + 1
				extantLines <- c(extantLines[-c(v_i, u_i)], alpha)
				lstates[alpha,] = mstates[alpha,] <- palpha
				heights[alpha] <- h1
				
				uv <- c(u, v) 
				inEdgeMap[uv] <- lineageCounter
				outEdgeMap[alpha,] <- uv
				parent[uv] <- alpha
				parentheights[uv] <- h1
				daughters[alpha,] <- uv
				edge[u,] <- c(alpha, u) # 
				edge.length[u] <- h1 - heights[u] 
				edge[v,] <- c(alpha,v) # 
				edge.length[v] <- h1 - heights[v] 
		}
	}
	 
	self<- list(edge=edge, edge.length=edge.length, Nnode=Nnode, tip.label=tip.label, heights=heights, parentheights=parentheights, parent=parent, daughters=daughters, lstates=lstates, mstates=mstates, ustates=ustates, m = m, sampleTimes = sampleTimes, sampleStates= sampleStates, maxSampleTime=maxSampleTime, inEdgeMap = inEdgeMap, outEdgeMap=outEdgeMap)
	class(self) <- c("binaryDatedTree", "phylo")
	
	#~ <reorder edges for compatibility with ape::phylo functions> 
	#~ (ideally ape would not care about the edge order, but actually most functions assume a certain order)
	sampleTimes2 <- sampleTimes[names(sortedSampleHeights)]
	sampleStates2 <- lstates[1:n,]; rownames(sampleStates2) <- tip.label
#~ browser()
	phylo <- read.tree(text=write.tree(self) )
	sampleTimes2 <- sampleTimes2[phylo$tip.label];
	sampleStates2 <- sampleStates2[phylo$tip.label,]; 
	bdt <- binaryDatedTree(phylo, sampleTimes2, sampleStates = sampleStates2)
	bdt
}




#~ setOldClass("binaryDatedTree")
#~ plot.binaryDatedTree <- function(x,...) { 
#~ 	plot.phylo( read.tree(text=write.tree(x))  , ...)
#~ }
#~ setMethod("plot.phylo", signature(x="binaryDatedTree"), plot.binaryDatedTree)



#######################################
# SUMMARY STATISTIC GENERATOR
# TODO adapt for heterochronous samples
#~ m + m2+ m3
#' Calculate cluster size moments from model
#' @param sampleTime 
#' @param sampleStates	
#' @param maxTime
#' @param minTime		
#' @param timeResolution
#' @param discretizeRates
#' @param fgyResolution
#' @param integrationMethod
#' @export
#'
calculate.cluster.size.moments.from.model <- function(sampleTime, sampleStates, FGY = NULL , maxTime=NA, minTime = NA, timeResolution = 50, discretizeRates=FALSE, fgyResolution = 100 , integrationMethod = 'adams')
{
	require(deSolve)
	globalParms <- list( USE_DISCRETE_FGY = discretizeRates ,  INTEGRATIONMETHOD=integrationMethod )
	if (!is.null(FGY)) { globalParms <- modifyList( globalParms, FGY) }
	else {globalParms <- modifyList( globalParms, list( F. = F., G. = G., Y. = Y. ))} 
	
	n 		<- nrow(sampleStates) 
	m 		<- ncol(sampleStates)
	sampleTimes <- rep(sampleTime, n)
	if (is.na(maxTime)) maxTime <- sampleTime 
	if (is.na(minTime)) minTime <- 0
	heights <- seq( 0, maxTime-minTime, length.out = timeResolution)
	

	if (discretizeRates) 
	{ 
		globalParms <- .start.discrete.rates(fgyResolution, maxSampleTime = max(sampleTimes), globalParms = globalParms) 
	}
	
	
	FGY_INDEX <- 1 
	
	get.fgy <- function(h,t=NULL)
	{
		if (globalParms$USE_DISCRETE_FGY)
		{
			FGY_INDEX <- min( 1+floor( globalParms$FGY_RESOLUTION * h/ globalParms$maxHeight ), globalParms$FGY_RESOLUTION)
			.Y		<- globalParms$Y.(FGY_INDEX)
			.F		<- globalParms$F.(FGY_INDEX)
			.G		<- globalParms$G.(FGY_INDEX)
		} else
		{
			.Y		<- globalParms$Y.(t)
			.F		<- globalParms$F.(t)
			.G		<- globalParms$G.(t)
		}
		list(.F=.F, .G=.G, .Y=.Y,  FGY_INDEX=FGY_INDEX)
	}
	
	dXiXjXij <- function(h, XiXjXij, parms, ...)
	{
		#X are vectors indexed by ancestor type
		Xi	<- XiXjXij[1:m]
		Xj 	<- XiXjXij[(m+1):(2*m)]
		Xij	<- XiXjXij[(2*m+1):(length(XiXjXij))]
		t 	<- parms$maxTime - h
		with(get.fgy(h, t), 
		{
			Xi_Y 	<- Xi / .Y
			Xj_Y 	<- Xj / .Y
			Xij_Y 	<- Xij / .Y
			csFpG 	<- colSums( .F + .G )
			dXi 	<- .F %*% Xi_Y + .G %*% Xi_Y - csFpG * Xi_Y
			dXj 	<- .F %*% Xj_Y + .G %*% Xj_Y - csFpG * Xj_Y
			dXij 	<- .F %*% Xij_Y + .G %*% Xij_Y - csFpG * Xij_Y  + (.F %*% Xj_Y )*Xi_Y  + (.F %*% Xi_Y ) * Xj_Y
			list( c( dXi, dXj, dXij ) )
		})
	}
	
	dA <- function(h, A, parms, ...)
	{
		with(get.fgy(h), 
		{ 
			A_Y 	<- A / .Y
			csFpG 	<- colSums( .F + .G )
			list( .G %*% A_Y - csFpG * A_Y + (.F %*% A_Y) * pmax(1-A_Y, 0) )
		})
	}
	
	
	# solve for A
	A0						<- colSums(sampleStates) 
	parameters 				<- list(maxTime = maxTime)
	o 						<- ode(y=A0, times=heights, func=dA, parms=parameters,  method=globalParms$INTEGRATIONMETHOD) 
	FGY_INDEX 	<- 1
	A 						<- o[,2:(ncol(o))]
	
	# solve Xi Xj and Xij at same time, 3m variables, avoids approxfun nastiness
	# solve second aggregated moment X2 & moments
	Mij_h 		<- array(0, dim = c(m, m, timeResolution), dimnames=list(paste('state',1:m,sep='.'),paste('state',1:m,sep='.'),c()) )
	Mi_h		<- matrix(0, m, timeResolution,dimnames=list(paste('state',1:m,sep='.'),c()))
	for (i in 1:m)
	{ # first tip type
		for (j in i:m)
		{ #second tip type
			#if (i!=j){
			parameters 		<- list(maxTime = maxTime)#, globalParms = globalParms) 
			Xij_h0 			<- rep(0, m)
			if(i==j) 
				Xij_h0[i]	<- A0[i]
			XiXjXij_h0 		<- c( diag(m)[i,] * A0, diag(m)[j,] * A0, Xij_h0) # each X is vector m(ancestor type) X 1
			o 				<- ode(y=XiXjXij_h0, times=heights, func=dXiXjXij,parms=parameters,  method=globalParms$INTEGRATIONMETHOD)
			# moments
			FGY_INDEX 	<- 1
			Xij_h 			<- o[, (1 + 2*m+1) : ncol(o) ] 
			if(i==j)
			{
				Xi_h		<- o[, 1+1:m]	
				tryCatch( 
					Mi_h[i,]<- rowSums( Xi_h   ) / rowSums(A)	# why avg over cols? / total lineages
							, error = function(e) browser() )
			}		
			# covariances
			tryCatch( 					
				Mij_h[i,j,]	<- rowSums( Xij_h   ) / rowSums(A) 	# timeResolution X 1
				, error = function(e) browser() )
			Mij_h[j,i,] 	<- Mij_h[i,j,]			
		}
	}
	
	if (discretizeRates) { globalParms <- .end.discrete.rates(globalParms)}
	list( heights=heights, Mi= Mi_h, Mij= Mij_h, A = A)
}

calculate.cluster.size.distr.from.tree<- function(bdt, heights)
{
	require(data.table)
	m 		<- ncol( bdt$sampleStates )	
	is.tip 	<- function(bdt, u) 
	{ 
		ifelse(is.na(bdt$daughters[u,1]), TRUE, FALSE)
	}
	identify.tips <- function(bdt, u) 
	{
		if (is.tip(bdt, u)) return(u)
		uv	<- bdt$daughters[u,]
		c(identify.tips(bdt, uv[1]), identify.tips(bdt, uv[2]) ) 
	}
	calculate.sizes.u <- function(bdt, u) 
	{
		tips	<- identify.tips(bdt, u)
		colSums( bdt$sampleStates[tips,,drop=FALSE] )									
	}
	calculate.sizes <- function(bdt, extant)
	{
		if(length(extant)==0)	return( matrix(NA,ncol( bdt$sampleStates ),0) )
		sapply( extant, function(u) calculate.sizes.u(bdt, u) ) 
	}		
	distr	<- lapply( heights, function(h)
			{
				extant 				<- .extant.at.height(h, bdt) 
				sizes_h				<- calculate.sizes(bdt, extant)
				rownames(sizes_h)	<- paste('state',1:m,sep='.')				
				as.data.table( t( rbind(sizes_h, height=ifelse(ncol(sizes_h),h,numeric(0)), ntip=ifelse(ncol(sizes_h),colSums(sizes_h),numeric(0))) ) )
			})
	#print(sapply(distr, function(x) ncol(x) ))
	#print(tail(distr,1))
	distr	<- do.call("rbind",distr)
	distr	
}

#' Calculate cluster size moments from tree
#' @param bdt		binary dated tree
#' @param heights	vector numeric, heights at which to calculate cluster sizes	
#' @export
#'
calculate.cluster.size.moments.from.tree <- function(bdt, heights)
{
# bdt : binaryDatedTree
# heights : vector numeric, heights at which to calculate cluster sizes
	m 		<- ncol( bdt$sampleStates )
	Mi		<- matrix(0, m, length(heights))
	Mij 	<- array(0, dim=c(m, m, length(heights)) )
	is.tip 	<- function(u) 
	{ 
		ifelse(is.na(bdt$daughters[u,1]), TRUE, FALSE)
	}
	identify.tips <- function(u) 
	{
		if (is.tip(u)) return(u)
		uv	<- bdt$daughters[u,]
		c(identify.tips(uv[1]), identify.tips(uv[2]) ) 
	}
	calculate.moment1 <- function(u, i) 
	{
		tips	<- identify.tips(u)
		sum( bdt$sampleStates[tips,i] )				
	}	
	calculate.mean.moment1 <- function(extant, i)
	{
		if(length(extant)==0)	return(NA)
		mean( sapply( extant, function(u) calculate.moment1(u, i) ) ) 
	}
	calculate.moment2 <- function(u, i, j) 
	{
		tips	<- identify.tips(u)
		itips 	<- sum( bdt$sampleStates[tips,i] )
		jtips 	<- sum( bdt$sampleStates[tips,j] )
		itips * jtips
	}
	calculate.mean.moment2 <- function(extant, i, j)
	{
		if (length(extant)==0) return(NA)
		mean( sapply( extant, function(u) calculate.moment2(u, i, j) ) ) 
	}
	for(ih in 1:length(heights))
	{
		h 		<- heights[ih]
		extant 	<- .extant.at.height(h, bdt) 
		if(length(extant) <=1)
		{
			for (i in 1:m)
			{
				Mi[i,ih]	<- Mi[i,ih-1]
				for (j in i:m)
				{
					Mij[i,j,ih] <- Mij[i,j,ih-1]
					Mij[j,i,ih] <- Mij[i,j,ih]
				}
			}
		} 
		else
		{
			for (i in 1:m)
			{
				Mi[i,ih]	<- calculate.mean.moment1(extant,i)
				for (j in i:m)
				{
					Mij[i,j,ih] <- calculate.mean.moment2(extant,i,j)
					Mij[j,i,ih] <- Mij[i,j,ih]
				}
			}
		}
	}
	
	list( Mi=Mi, Mij=Mij )
}
