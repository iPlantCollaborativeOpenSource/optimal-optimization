#The function ace was modified to use the optimx package and allow diagnostic tracking

## ace.R (2013-03-18)

##   Ancestral Character Estimation

## Copyright 2005-2013 Emmanuel Paradis and Ben Bolker

## This file is part of the R-package ape.
## See the file ../COPYING for licensing issues.
ace <-
function(x, phy, type = "continuous",
method = if (type == "continuous") "REML" else "ML",
CI = TRUE, model = if (type == "continuous") "BM" else "ER",
scaled = TRUE, kappa = 1, corStruct = NULL, ip = 0.1,
use.expm = FALSE, optimx.alg = "L-BFGS-B")
{
    if (!inherits(phy, "phylo"))
	stop('object "phy" is not of class "phylo"')
    if (is.null(phy$edge.length))
	stop("tree has no branch lengths")
    type <- match.arg(type, c("continuous", "discrete"))
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    if (nb.node != nb.tip - 1)
	stop('"phy" is not rooted AND fully dichotomous.')
    if (length(x) != nb.tip)
	stop("length of phenotypic and of phylogenetic data do not match.")
    if (!is.null(names(x))) {
        if(all(names(x) %in% phy$tip.label))
		x <- x[phy$tip.label]
        else warning("the names of 'x' and the tip labels of the tree do not match: the former were ignored in the analysis.")
    }
    obj <- list()
    if (kappa != 1) phy$edge.length <- phy$edge.length^kappa
    if (type == "continuous") {
        switch(method, "REML" = {
			   minusLogLik <- function(sig2) {
			   if (sig2 < 0) return(1e100)
			   V <- sig2 * vcv(phy)
## next three lines borrowed from dmvnorm() in 'mvtnorm'
			   distval <- stats::mahalanobis(x, center = mu, cov = V)
			   logdet <- sum(log(eigen(V, symmetric = TRUE, only.values = TRUE)$values))
			   (nb.tip * log(2 * pi) + logdet + distval)/2
			   }
			   mu <- rep(ace(x, phy, method="pic")$ace[1], nb.tip)
			   out <- nlm(minusLogLik, 1, hessian = TRUE)
			   sigma2 <- out$estimate
			   se_sgi2 <- sqrt(1/out$hessian)
			   tip <- phy$edge[, 2] <= nb.tip
			   minus.REML.BM <- function(p) {
			   x1 <- p[phy$edge[, 1] - nb.tip]
			   x2 <- numeric(length(x1))
			   x2[tip] <- x[phy$edge[tip, 2]]
			   x2[!tip] <- p[phy$edge[!tip, 2] - nb.tip]
			   -(-sum((x1 - x2)^2/phy$edge.length)/(2 * sigma2) -
				 nb.node * log(sigma2))
			   }
			   out <- nlm(function(p) minus.REML.BM(p),
						  p = rep(mu[1], nb.node), hessian = TRUE)
			   obj$resloglik <- -out$minimum
			   obj$ace <- out$estimate
			   names(obj$ace) <- nb.tip + 1:nb.node
			   obj$sigma2 <- c(sigma2, se_sgi2)
			   if (CI) {
			   se <- .getSEs(out)
			   tmp <- se * qt(0.025, nb.node)
			   obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
			   }
			   }, "pic" = {
			   if (model != "BM")
			   stop('the "pic" method can be used only with model = "BM".')
## See pic.R for some annotations.
			   phy <- reorder(phy, "pruningwise")
			   phenotype <- numeric(nb.tip + nb.node)
			   phenotype[1:nb.tip] <- if (is.null(names(x))) x else x[phy$tip.label]
			   contr <- var.con <- numeric(nb.node)
			   ans <- .C("pic", as.integer(nb.tip), as.integer(nb.node),
						 as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]),
						 as.double(phy$edge.length), as.double(phenotype),
						 as.double(contr), as.double(var.con),
						 as.integer(CI), as.integer(scaled),
						 PACKAGE = "ape")
			   obj$ace <- ans[[6]][-(1:nb.tip)]
			   names(obj$ace) <- nb.tip + 1:nb.node
			   if (CI) {
			   se <- sqrt(ans[[8]])
			   tmp <- se * qnorm(0.025)
			   obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
			   }
			   }, "ML" = {
			   if (model == "BM") {
			   tip <- phy$edge[, 2] <= nb.tip
			   dev.BM <- function(p) {
			   if (p[1] < 0) return(1e100) # in case sigma is negative
			   x1 <- p[-1][phy$edge[, 1] - nb.tip]
			   x2 <- numeric(length(x1))
			   x2[tip] <- x[phy$edge[tip, 2]]
			   x2[!tip] <- p[-1][phy$edge[!tip, 2] - nb.tip]
			   -2 * (-sum((x1 - x2)^2/phy$edge.length)/(2*p[1]) -
					 nb.node * log(p[1]))
			   }
			   out <- nlm(function(p) dev.BM(p),
						  p = c(1, rep(mean(x), nb.node)), hessian = TRUE)
			   obj$loglik <- -out$minimum / 2
			   obj$ace <- out$estimate[-1]
			   names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
			   se <- .getSEs(out)
			   obj$sigma2 <- c(out$estimate[1], se[1])
			   if (CI) {
			   tmp <- se[-1] * qt(0.025, nb.node)
			   obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
			   }
			   }
			   }, "GLS" = {
			   if (is.null(corStruct))
			   stop('you must give a correlation structure if method = "GLS".')
			   if (class(corStruct)[1] == "corMartins")
			   M <- corStruct[1] * dist.nodes(phy)
			   if (class(corStruct)[1] == "corGrafen")
			   phy <- compute.brlen(attr(corStruct, "tree"),
									method = "Grafen",
									power = exp(corStruct[1]))
			   if (class(corStruct)[1] %in% c("corBrownian", "corGrafen")) {
			   dis <- dist.nodes(attr(corStruct, "tree"))
			   MRCA <- mrca(attr(corStruct, "tree"), full = TRUE)
			   M <- dis[as.character(nb.tip + 1), MRCA]
			   dim(M) <- rep(sqrt(length(M)), 2)
			   }
			   varAY <- M[-(1:nb.tip), 1:nb.tip]
			   varA <- M[-(1:nb.tip), -(1:nb.tip)]
			   V <- corMatrix(Initialize(corStruct, data.frame(x)),
							  corr = FALSE)
			   invV <- solve(V)
			   o <- gls(x ~ 1, data.frame(x), correlation = corStruct)
			   GM <- o$coefficients
			   obj$ace <- drop(varAY %*% invV %*% (x - GM) + GM)
			   names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
			   if (CI) {
			   se <- sqrt((varA - varAY %*% invV %*% t(varAY))[cbind(1:nb.node, 1:nb.node)])
			   tmp <- se * qnorm(0.025)
			   obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
			   }
			   })
    } else { # type == "discrete"
        if (method != "ML")
		stop("only ML estimation is possible for discrete characters.")
        if (!is.factor(x)) x <- factor(x)
        nl <- nlevels(x)
        lvls <- levels(x)
        x <- as.integer(x)
		if (is.character(model)) {
            rate <- matrix(NA, nl, nl)
            if (model == "ER") np <- rate[] <- 1
            if (model == "ARD") {
                np <- nl*(nl - 1)
                rate[col(rate) != row(rate)] <- 1:np
            }
            if (model == "SYM") {
                np <- nl * (nl - 1)/2
                sel <- col(rate) < row(rate)
                rate[sel] <- 1:np
                rate <- t(rate)
                rate[sel] <- 1:np
            }
        } else {
            if (ncol(model) != nrow(model))
			stop("the matrix given as 'model' is not square")
            if (ncol(model) != nl)
			stop("the matrix 'model' must have as many rows as the number of categories in 'x'")
            rate <- model
            np <- max(rate)
        }
        index.matrix <- rate
        tmp <- cbind(1:nl, 1:nl)
        index.matrix[tmp] <- NA
        rate[tmp] <- 0
        rate[rate == 0] <- np + 1 # to avoid 0s since we will use this as numeric indexing
		
        liks <- matrix(0, nb.tip + nb.node, nl)
        TIPS <- 1:nb.tip
        liks[cbind(TIPS, x)] <- 1
        phy <- reorder(phy, "pruningwise")
		
		## E <- if (use.expm) expm::expm else ape::matexpo
        E <- if (use.expm) {
            library(expm)
            get("expm", "package:expm")
        } else ape::matexpo
		
        Q <- matrix(0, nl, nl)
        dev <- function(p, output.liks = FALSE) {
            if (any(is.nan(p)) || any(is.infinite(p))) return(1e50)
			## from Rich FitzJohn:
            comp <- numeric(nb.tip + nb.node) # Storage...
            Q[] <- c(p, 0)[rate]
            diag(Q) <- -rowSums(Q)
            for (i  in seq(from = 1, by = 2, length.out = nb.node)) {
                j <- i + 1L
                anc <- phy$edge[i, 1]
                des1 <- phy$edge[i, 2]
                des2 <- phy$edge[j, 2]
                v.l <- E(Q * phy$edge.length[i]) %*% liks[des1, ]
                v.r <- E(Q * phy$edge.length[j]) %*% liks[des2, ]
                v <- v.l * v.r
                comp[anc] <- sum(v)
                liks[anc, ] <- v/comp[anc]
            }
            if (output.liks) return(liks[-TIPS, ])
            dev <- -2 * sum(log(comp[-TIPS]))
            if (is.na(dev)) Inf else dev
        }
        
        #Optimization 
#       out <- nlminb(rep(ip, length.out = np), function(p) dev(p),lower = rep(0, np), upper = rep(1e50, np)) #original
    #Bounded    
		if(optimx.alg == "bobyqa" | optimx.alg == "L-BFGS-B" | optimx.alg == "nlminb" | optimx.alg == "spg" | optimx.alg =="Rcgmin" | optimx.alg == "Rvmmin"){ 
			begin.time <- proc.time()
			out <- optimx(rep(ip, length.out = np), dev, lower = rep(0, np), upper = rep(1e50, np), method=optimx.alg)
			time <- as.numeric(proc.time()[3]-begin.time[3])/(60)
		}else{
		#Unbounded  
			begin.time <- proc.time()
			out <- optimx(rep(ip, length.out = np), dev, method=optimx.alg)
			time <- as.numeric(proc.time()[3]-begin.time[3])/(60)			
		}
		obj$loglik <- -out$fvalues$fvalues/2
		obj$rates <- out$par$par
		obj$time <- time
		obj$itns <- out$itns$itns
		obj$method <- out$method$method
#       oldwarn <- options("warn")
#       options(warn = -1)
#       out.nlm <- nlm(function(p) dev(p), p = obj$rates, iterlim = 1,
#                       stepmax = 0, hessian = TRUE)
#       options(oldwarn)
#       obj$se <- .getSEs(out.nlm)
#       obj$index.matrix <- index.matrix
#       if (CI) {
#           obj$lik.anc <- dev(obj$rates, TRUE)
#           colnames(obj$lik.anc) <- lvls
#       }
		return(obj)
		
	}
}
#   obj$call <- match.call()
#   class(obj) <- "ace"
