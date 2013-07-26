#The function fitDiscrete was modified to use the optimx package and allow diagnostic tracking

fitDiscreteMod <- function(phy, data, model=c("ER", "SYM", "ARD"), treeTransform=c("none", "lambda", "kappa", "delta", "linearChange", "exponentialChange", "twoRate"), data.names=NULL, plotlnl=FALSE, qLimits=c(0.0001, 1000), pLimits=c(0.00001, 10), ip, optimx.alg)
{
	
	model<-match.arg(model)
	
	treeTransform<-match.arg(treeTransform)
	
	res<-list()
	for(i in 1:1)
	{
		if(treeTransform=="none") {
			f<-function(x) {
				likelihoodDiscrete(phy, data[,i], exp(x), model)
			}
			nep=0; pLow=-10; pHigh=log(1); pStart=NULL;		
		}	
		if(treeTransform=="lambda") {
			f<-function(x) {
				likelihoodDiscrete(phy, data[,i], exp(x[-1]), lambda=exp(x[1]), model=model)
			}	
			nep=1; pLow=-10; pHigh=log(1); pStart=0.1;
		}
		if(treeTransform=="delta") {
			f<-function(x) {
				if(x[1]<= -10.000000 | x[1]>= 2.302585 | x[2] <= -5.831357 | x[1] >= 5.681568){
					return(1000000)
				}
				likelihoodDiscrete(phy, data[,2], exp(x[-1]), delta=exp(x[1]), model=model)
			}
			nep=1; pLow=-10; pHigh=log(10);pStart=ip;
		}
		if(treeTransform=="kappa") {
			f<-function(x) {
				likelihoodDiscrete(phy, data[,i], exp(x[-1]), kappa=exp(x[1]), model=model)
			}
			nep=1; pLow=-10; pHigh=log(1);pStart=0.1;
		}
		if(treeTransform=="linearChange") {
			f<-function(x) {
				likelihoodDiscrete(phy, data[,i], exp(x[-1]), endRate=exp(x[1]), linear=T, model=model)
			}
			nep=1; pLow=-10; pHigh=log(10);pStart=0.1;
		}
		if(treeTransform=="exponentialChange") {
			f<-function(x) {
				likelihoodDiscrete(phy, data[,i], exp(x[-1]), endRate=exp(x[1]), model=model)
			}
			nep=1; pLow=-10; pHigh=log(10);pStart=0.1;
		}
		if(treeTransform=="twoRate") {
			f<-function(x) {
				likelihoodDiscrete(phy, data[,i], exp(x[-(1:2)]), breakPoint=x[1], endRate=exp(x[2]), model=model)
			}
			mv<-max(branching.times(phy))	
			nep=2; pLow=c(mv/1000, -10); pHigh=c(mv, 10);pStart=c(0.1, 0.1);
		}
		nRateCats<-getRateCats(data[,i], model)
		
		outTries<-list()
		totalbl<-sum(phy$edge.length)
		minQ=log(0.01/totalbl)
		maxQ=log(1000/totalbl)
		ntries<-20
		ltry<-numeric(ntries)
		ltry[]<-NA
		lsol<-matrix(nrow= ntries, ncol=nRateCats+nep)
		sp<-numeric(nRateCats)
		qTries<-exp(-7:2)
		
		if(nep==0) {
			lower=rep(minQ, nRateCats)
			upper=rep(maxQ, nRateCats)
		} else {
			lower=c(pLow, rep(minQ, nRateCats))
			upper=c(pHigh, rep(maxQ, nRateCats))
		}
		
#			cat("Finding the maximum likelihood solution\n")
#			cat("[0        50      100]\n")
#			cat("[")
		
		if (optimx.alg == "bobyqa" | optimx.alg == "L-BFGS-B" | optimx.alg == "nlminb" | optimx.alg == "spg" | optimx.alg =="Rcgmin" | optimx.alg == "Rvmmin"){
			
			begin.time <- proc.time()
			for(j in 1:10) {
				sp<-c(pStart, log(rep(qTries[j], nRateCats)))
#				te<-try(outTries[[j]]<-optim(f, par=sp, method="L",  lower=lower, upper=upper), silent=T)
				te<-try(outTries[[j]]<-optimx(f, par=sp, method=optimx.alg, lower=lower, upper=upper), silent=TRUE)
				if(class(te)!="try-error") {
					ltry[j]<-outTries[[j]]$fvalues$fvalues
					lsol[j,]<-exp(outTries[[j]]$par$par)
				}
#				cat(".")
				
			}
			for(j in 1:10) {
				sp<-c(pStart, runif(nRateCats, minQ, maxQ))
				te<-try(outTries[[10+j]]<-optimx(f, par=sp, method=optimx.alg, lower=lower, upper=upper), silent=TRUE)
				if(class(te)!="try-error") {
					ltry[10+j]<-outTries[[10+j]]$fvalues$fvalues
					lsol[10+j,]<-exp(outTries[[10+j]]$par$par)
				}
#				cat(".")
				
			}
			time <- as.numeric(proc.time()[3]-begin.time[3])/(60)
		}else{
			begin.time <- proc.time()
			for(j in 1:10) {
				sp<-c(pStart, log(rep(qTries[j], nRateCats)))
#				te<-try(outTries[[j]]<-optim(f, par=sp, method="L",  lower=lower, upper=upper), silent=T)
				te<-try(outTries[[j]]<-optimx(f, par=sp, method=optimx.alg), silent=TRUE)
				
				if(class(te)!="try-error") {
					ltry[j]<-outTries[[j]]$fvalues$fvalues
					lsol[j,]<-exp(outTries[[j]]$par$par)
				}
#				cat(".")
				
			}
			for(j in 1:10) {
				sp<-c(pStart, runif(nRateCats, minQ, maxQ))
				te<-try(outTries[[10+j]]<-optimx(f, par=sp, method=optimx.alg), silent=TRUE)
				if(class(te)!="try-error") {
					ltry[10+j]<-outTries[[10+j]]$fvalues$fvalues
					lsol[10+j,]<-exp(outTries[[10+j]]$par$par)
				}
#				cat(".")
				
			}
			time <- as.numeric(proc.time()[3]-begin.time[3])/(60)
		}
		
		
#		cat("]\n")
		ok<-!is.na(ltry)
		b<-min(which(ltry==min(ltry[ok])))
		out<-outTries[[b[1]]]	
		
#		if(treeTransform=="none" & model=="ER") {
#			res[[i]]<-list(lnl=-out$fvalues$fvalues, q=-getQ(exp(out$par$par), nRateCats, model)[1,1],time)
#		} else if(treeTransform=="none" & model!="ER") {
#			res[[i]]<-list(lnl=-out$fvalues$fvalues, q=getQ(exp(out$par$par), nRateCats, model),time)
#			
#		} else if(treeTransform=="twoRate") {
#			res[[i]]<-list(lnl=-out$value, q=getQ(exp(out$par[-(1:2)]), nRateCats, model), breakpoint=out$par[1], endRate=exp(out$par[2]), message=out$message)
#		} else 	res[[i]]<-list(lnl=-out$fvalues$fvalues, q=getQ(exp(out$par$par[-1]), nRateCats, model), treeParam=exp(out$par$par[1]), time=time)
		
		res<-NULL
		res$loglik <- -out$fvalues$fvalues
		res$rates <- exp(out$par$par[1])
		res$time <- time
		res$method<- out$method$method
		
#		if(!is.null(colnames(data))) names(res)[i]<-colnames(data)[i] else names(res)[i]<-paste("Trait", i, sep="")
		
		if((model=="SYM" | model=="ARD") & nRateCats>2) {
			cat("Plotting surfaces currently not supported for SYM and ARD models unless your character has only two states.\n")
			plotLnlSurf=F
		}
		
		if(plotlnl) {
			cat("Calculating surface\n")
			if(qLimits[1]<=0) {
				cat("Q must be positive, resetting lower plotting limit to 0.00000001")
				qLimits[1]=0.00000001
			}
			if(treeTransform=="none") {
				if(model=="ER") {
					qx<-exp(seq(log(qLimits[1]), log(qLimits[2]), length=50))
					lnl<-numeric(50)
					for(j in 1:50)
					lnl[j]<- -f(log(qx[j]))
					
					lnlDiff<- -out$value-lnl
					plot(qx, lnl, log="x", type="l", xlab="Rate (q)", ylab="lnL")
				} else {
					qx<-exp(seq(log(qLimits[1]), log(qLimits[2]), length=20))
					qy<-exp(seq(log(qLimits[1]), log(qLimits[2]), length=20))
					lnl<-matrix(nrow=20, ncol=20)
					for(j in 1:20)
					for(k in 1:20)
					lnl[j,k]<- -f(log(c(qx[j], qy[k])))
					
					lnlDiff<- -out$value-lnl
					contour(qx, qy, lnl, xlab="Forward Rate", ylab="Backward Rate")
				}	
			} else {
				px<-exp(seq(log(pLimits[1]), log(pLimits[2]), length=20))
				qy<-exp(seq(log(qLimits[1]), log(qLimits[2]), length=20))
				lnl<-matrix(nrow=20, ncol=20)
				for(j in 1:20)
				for(k in 1:20)
				try(lnl[j,k]<- -f(log(c(px[j], qy[k]))))
				
				lnlDiff<- -out$value-lnl
				contour(px, qy, lnlDiff, levels=c(1, 2, 3, 4, 8, 10, 20, 30, 40, 50, 100, 500, 1000), xlab="Tree Transform parameter estimate", ylab="Rate (q)", main="lnL Surface")		
			}
			
		}
		
	}
	return(res)
}


###Felsensteins pruning algorithm
### Modified June 13 2008 to work on log-likelihoods instead of likelihoods
likelihoodDiscrete<-function(phy, tip.data, q, model="ER", delta=1, lambda=1,  kappa=1, endRate=1, linear=F, breakPoint=0, f=1, rtt.rescale=0, total.rescale=F, returnFull=F)
{
	if(!is.factor(tip.data)) tip.data<-factor(tip.data)
	Q<-getQ(q, n=nlevels(tip.data), model=model)
	if (class(phy) != "phylo")
	stop("object \"phy\" is not of class \"phylo\"");
	#new2old.phylo(phy)->phy ##converts back to old ape tree format with negative values denoting internal nodes
	tmp <- as.numeric(phy$edge)
	nb.tip <- max(tmp) #number of tips
	nb.node <- -min(tmp) #number of internal nodes
	nb.states <- nlevels(tip.data) #numbers of states in character
	l <- matrix(-Inf, nrow=nrow(phy$edge), ncol=nb.states) #makes matrix that will store likelihoods of each character state at each node
	root <- numeric(nb.states) #makes vector that will store root likelihoods
	m <- 1:length(phy$tip.label) ##identifies elements of tip.data matrix that corresponds with the tip.label
	if (delta != 1)
	deltaTree(phy, delta) -> phy;
	if (lambda != 1)
	lambdaTree(phy, lambda) -> phy;
	if (kappa != 1)
	kappaTree(phy, kappa) -> phy;
	if (endRate != 1) {
		if(breakPoint!=0) {
			tworateTree(phy, breakPoint, endRate) -> phy;
		} else if(linear==T) {
			linearchangeTree(phy, endRate) -> phy;
		} else exponentialchangeTree(phy, endRate)->phy;
	}
	
	
	
#When comparing deltas across different qs, it might be useful to rescale the total tree length to one	
	if(rtt.rescale!=0)	
	rescaleTree(phy, rtt.rescale) -> phy
	
	new2old.phylo(phy)->phy
	
	for(i in 1:nrow(phy$edge)) #for each edge
	if(as.numeric(phy$edge[i,2])>0) {
		
		l[i,tip.data[m[as.numeric(phy$edge[i,2])]]] <- 0 
	}
#if the edge is connected to a terminal taxon, you set the likelihood of the tip value equal to 1 and all others equal to zero.
	times <- branching.times(old2new.phylo(phy)) #get node to tip distances
	-1*(1:(max(as.numeric(names(times)))-min(as.numeric(names(times)))+1))->names(times)
	times = max(times) - times #convert into root to node tips
	if(total.rescale) {
		sum(phy$edge.length) -> total.tree
		phy$edge.length <- phy$edge.length/total.tree
	}
	while(1) {
		
		if(sum(as.numeric(phy$edge[,2])>0)==2) break 
		
		#obtain ancestors of current tips
		x <- match(phy$edge[,1], phy$edge[,2])[as.numeric(phy$edge[,2])>0] #finds nodes connected to terminal taxa
		#find last node with two tip descendent
		a <- max(x[duplicated(x)])
		t <- which(phy$edge[,1]==phy$edge[a,2])
		bl <- phy$edge.length[t]
		age = times[which(names(times)==phy$edge[a,2])]
		l[a,] <- frag.like(l[t,], bl, Q)
		#next line effectively prunes out the tips just used
		phy$edge[a,2]<-1
		phy$edge[t,2]<-0
		
		
	}
	t <- which(as.numeric(phy$edge[,2])>0)
	bl <- phy$edge.length[t]
	root <- frag.like(l[t,], bl, Q)
	neglnl=-logspace_sum(root)+log(nb.states)
	if(returnFull==F) {
		return(neglnl)
	} else return(list(neglnl=neglnl, root=root, l=l))
}

getQ<-function(q, n, model)
{
	if(model=="ER") Q=evenQ(n)*q
	if(model=="SYM") {
		if(length(q)!=n*(n-1)/2) stop("You must supply the correct number of rate categories.")
		Q<-diag(n)
		xx=1
		for(i in 2:n) {
			for(j in 1:(i-1)) {
				Q[i,j]<-Q[j,i]<-q[xx]
				xx<-xx+1
			}
		}	
		for(i in 1:n) diag(Q)[i]<- -sum(Q[i,-i])
	}
	
	if(model=="ARD") {
		if(length(q)!=n*(n-1)) stop("You must supply the correct number of rate categories.")
		Q<-diag(n)
		xx=1
		for(i in 1:n) {
			for(j in (1:n)[-i]) {
				Q[i,j]<-q[xx]
				xx<-xx+1
			}
		}	
		for(i in 1:n) diag(Q)[i]<- -sum(Q[i,-i])
	}
	
	return(Q)
}

getRateCats<-function(data, model)
{
	if(model=="ER") return(1)
	n<-nlevels(factor(data))
	if(model=="SYM") return(n*(n-1)/2)
	if(model=="ARD") return(n*(n-1))
	
}

##evenQ is an internal function of GEIGER
##This function makes a template for the calculate of a rate matrix that will set all transitions to the same value (based on the total number of states).  One can then multiply the resulting matrix by the overall rate to get the rate matrix for a particular analysis.


`evenQ` <-
function(n)

{
	
	q<--diag(n)
	
	q[q==0]<-1/(n-1)
	
	return(q)
	
}

#This function calculates the likelihood on one branch of a tree

`frag.like` <-
function(tip.like, bl, q)

{
	
	nb.states<-ncol(tip.like)
	
	r<-rep(0, nb.states)
	
	d<-length(bl)
	
	p<-list(d)
	
	for(i in 1:d)
	
	p[[i]]<-MatrixExp.eig(q*bl[i])
	
	for(i in 1:nb.states)
	
	for(j in 1:d) 
	
	r[i]<-r[i]+logspace_sum(log(p[[j]][i,])+tip.like[j,])
	
	return(r)
	
}

#This function is required so that the matrix exponentiation never blows up during the likelihood calculation - but it only works for  symmetric matrices (ie evenQ matrices)


`MatrixExp.simple` <-
function(Q)
{
	n<-nrow(Q)
	res<-matrix(0, nrow=n, ncol=n)
	q<-Q[1,2]
	for(i in 1:n)
	res[i, i]<-1/n+(n-1)/n*exp(-n*q)
	res[res==0]<-1/n-1/n*exp(-n*q)
	return(res)
}

## Code lifted from ape function "ace"
MatrixExp.eig<-
function(Q)
{
	tmp <- eigen(Q, symmetric = FALSE)
	P1 <- tmp$vectors %*% diag(exp(tmp$values)) %*% solve(tmp$vectors)
	return(P1)
	
}

hasZeroLengthTips<-function(phy)
{
	ntips<-length(phy$tip.label)
	tips<-phy$edge[,2]<=ntips
	nn<-phy$edge.length[tips]==0
	if(sum(nn)>0) return(T)
	return(F)
}


logspace_add<-function(logx, logy) {
	if(logx==-Inf) return(logy) else max(logx, logy) + log1p(exp (-abs (logx - logy)));
}

logspace_sum<-function(logx) {
	r<-logx[1]
	if(length(logx)>1)
	for(i in 2:length(logx))
	r<-logspace_add(r, logx[i])
	r	
}