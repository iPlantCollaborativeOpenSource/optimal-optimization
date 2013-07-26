#Wrapper for easier use of ouch 

#written by Jeremy M. Beaulieu

ouch.wrap<-function(tree, trait, model=c("brown", "ou1", "ousm"), ip, optimx.alg){

	ot <- ape2ouch(tree)
	otd <- as(ot,"data.frame")
	trait <- data.frame(regime=trait[,2], phenoTrait=trait[,3], labels=trait[,1])
	otd <- merge(otd,trait, by="labels", all=TRUE)
	rownames(otd) <- otd$nodes
	ot <- with(otd,ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))
	
	if (is.character(model)) {
		
		if (model == "brown"){
			#evaluates a brownian motion model
			obj <- brown(tree=ot,data=otd["phenoTrait"])
		}
		#Note I only modified this option:
		if (model == "ou1"){
			#evaluates an OU model with a single, global selective regime
			otd$regimes <- as.factor("global")
			obj <- hansen(tree=ot,data=otd["phenoTrait"],regimes=otd["regimes"],sqrt.alpha=ip[1],sigma=ip[2], method=optimx.alg)
		}
		
		if (model == "ousm"){
			#evaluates an OU model with multiple selective regimes
			nb.nodes <- tree$Nnode
			otd$regime[1:nb.nodes] <- otd$labels[1:nb.nodes]
			otd$regime <- as.factor(otd$regime)
			obj <- hansen(tree=ot,data=otd["phenoTrait"],regimes=otd["regime"],sqrt.alpha=1,sigma=1, maxit=10000)
		}
		
	}
	summary(obj)
}
