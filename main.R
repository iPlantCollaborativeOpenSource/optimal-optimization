##Master controller for Testing Optimizers##

#written by Jeremy M. Beaulieu


#----------------------------------#
#-----     Dependencies       -----#
#----------------------------------#	

library(ape)
library(optimx)
library(data.table)
library(mvtnorm)

#Installs our modified version of ouch:
detach("package:ouch", unload=TRUE)
install.packages("dependencies/ouch_0.00.tar.gz",repos=NULL,type="source",lib="/tmp")
library(ouch, lib.loc="/tmp")

source("dependencies/ouch.wrap.R")
source("dependencies/ace.modified.R")
source("dependencies/fitContinuous.modified.R")
source("dependencies/fitDiscrete.modified.R")

#----------------------------------#
#-----     Test functions     -----#
#----------------------------------#	

###Function for testing ace discrete###
DiscreteAce <- function(phy,data){
	data<-data.frame(data[,2], data[,2], row.names=data[,1])
	data<-data[phy$tip.label,]
	optimx.alg<-c("spg","Rcgmin","Rvmmin","bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")
	res<-c()
	#20 points evenly spaced between the two ends of 1st and 99th quantile of an exponential centered on 0.1:
	ip.vect<-round(seq(0.001,0.50,by=.005),2)
	ip.vect[1]<-0.001
	subSet<-seq(5,100,by=5)
	ip.vect<-ip.vect[subSet]
	for(j in 1:12){
		for(i in 1:length(ip.vect)){
			tmp<-Inf
			try(tmp<-ace(data[,1],phy,type='discrete',ip=ip.vect[i],optimx.alg=optimx.alg[j]),silent=TRUE)
			if(is.list(tmp)){
				res<-rbind(res,c(ip.vect[i], tmp$loglik,tmp$rates,tmp$time, tmp$method))
			}
		}
	}
	res<-as.data.table(res)
	return(res)
}

###Function for testing fitDiscrete###
FitDiscTest <- function(phy,data){
	data<-data.frame(data[,2], data[,2], row.names=data[,1])
	data<-data[phy$tip.label,]
	optimx.alg<-c("spg","Rcgmin","Rvmmin","bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")
	res<-c()
	#20 points evenly spaced between the two ends of 1st and 99th quantile of an exponential centered on 0.1:
	ip.vect<-round(seq(0.001,0.50,by=.005),2)
	ip.vect[1]<-0.001
	subSet<-seq(5,100,by=5)
	ip.vect<-ip.vect[subSet]
	for(j in 1:12){
		for(i in 1:length(ip.vect)){
			tmp<-Inf
			try(tmp<-fitDiscreteMod(phy,data, model="ER", treeTransform="delta", ip=ip.vect[i], optimx.alg=optimx.alg[j]),silent=TRUE)
			if(is.list(tmp)){
				res<-rbind(res,c(ip.vect[i], tmp$loglik, tmp$rates, tmp$time, tmp$method))
			}
		}
	}
	res<-as.data.table(res)
	return(res)
}

###Function for testing fitContinuous###
FitContTest <- function(phy,data){
	data<-data.frame(data[,2], data[,2], row.names=data[,1])
	data<-data[phy$tip.label,]
	optimx.alg<-c("spg","Rcgmin","Rvmmin","bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")
	res<-c()
	#20 points evenly spaced between the two ends of 1st and 99th quantile of an exponential centered on 0.5:
	ip.vect<-round(seq(0.005,2.50,by=.025),2)
	ip.vect[1]<-0.01
	subSet<-seq(1,100,by=5)
	ip.vect<-ip.vect[subSet]
	for(j in 1:12){
		for(i in 1:length(ip.vect)){
			tmp<-Inf
			try(tmp<-fitContinuousMod(phy,data[,1],ip=ip.vect[i],optimx.alg=optimx.alg[j]),silent=TRUE)
			if(is.list(tmp)){
				res<-rbind(res,c(ip.vect[i], tmp$loglik, tmp$beta, tmp$delta, tmp$time, tmp$method))
			}
		}
	}
	res<-as.data.table(res)
	return(res)
}

#Function for testing ouch:
OuchTest <- function(phy,data){
	optimx.alg<-c("spg","Rcgmin","Rvmmin","bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")
	res<-c()
	#20 points evenly spaced 0.01 and 9.5:
	ip.vect<-exp(seq(log(0.01),log(20),by=1))
	ip.vect[1]<-0.01
	ip.vect.grid<-expand.grid(ip.vect,ip.vect)
	for(j in 1:12){
		for(i in 1:length(ip.vect.grid[,1])){
			tmp<-Inf
			try(tmp<-ouch.wrap(phy,data,model="ou1", ip=c(ip.vect.grid[i,1],ip.vect.grid[i,2]), optimx.alg=optimx.alg[j]),silent=TRUE)
			if(is.list(tmp)){
				res<-rbind(res,c(ip.vect[i], tmp$loglik, tmp$alpha, tmp$sigma.squared, tmp$optimizer.message, optimx.alg[j]))
			}
		}
	}
	res<-as.data.table(res)
	return(res)
}


#----------------------------------#
#-----           Run          -----#
#----------------------------------#

#####Geospiza####
tree<-read.tree("examples/geospiza.tre")
trait.tot<-read.delim("examples/geospiza.traits")
trait<-trait.tot[,c(1,2)]
pp1<-DiscreteAce(tree,trait)
save(pp1,file="aceDISCRETE_Geo.Rsave")
pp2<-FitDiscTest(tree,trait)
save(pp2,file="fitDISC_Geo.Rsave")
trait<-trait.tot[,c(1,3)]
pp3<-FitContTest(tree,trait)
save(pp3,file="fitCONT_Geo.Rsave")
trait<-trait.tot[,c(1:3)]
pp4<-OuchTest(tree,trait)
save(pp4,file="ouch_Geo.Rsave")

#####Aquilegia######
tree<-read.tree("examples/Aquilegia.new.tre")
trait.tot<-read.delim("examples/Aquilegia.traits")
trait<-trait.tot[,c(1,2)]
pp1<-DiscreteAce(tree,trait)
save(pp1,file="aceDISCRETE_Aqui.Rsave")
pp2<-FitDiscTest(tree,trait)
save(pp2,file="fitDISC_Aqui.Rsave")
trait<-trait.tot[,c(1,3)]
pp3<-FitContTest(tree,trait)
save(pp3,file="fitCONT_Aqui.Rsave")
pp4<-OuchTest(tree,trait.tot)
save(pp4,file="ouch_Aqui.Rsave")

#####Monocots#######
tree<-read.tree("examples/monocots.tre")
trait.tot<-read.delim("examples/monocots.gs")
trait<-trait.tot[,c(1,4)]
pp1<-DiscreteAce(tree,trait)
save(pp1,file="aceDISCRETE_Mono.Rsave")
pp2<-FitDiscTest(tree,trait)
save(pp2,file="fitDISC_Mono.Rsave")
trait<-trait.tot[,c(1,3)]
pp3<-FitContTest(tree,trait)
save(pp3,file="fitCONT_Mono.Rsave")
trait<-trait.tot[,c(1:3)]
pp4<-OuchTest(tree,trait)
save(pp4,file="ouch_Mono.Rsave")

#####Mammals#######
tree<-read.tree("examples/mammal.litter.tre")
trait<-read.delim("examples/mammal.litter")
pp1<-DiscreteAce(tree,trait)
save(pp1,file="aceDISCRETE_Mamm.Rsave")
pp2<-FitDiscTest(tree,trait)
save(pp2,file="fitDISC_Mamm.Rsave")
tree<-read.tree("examples/mammal.bs.tre")
trait.tot<-read.delim("examples/mammal.bs")
trait<-trait.tot[,c(1,3)]
pp3<-FitContTest(tree,trait)
save(pp3,file="fitCONT_Mamm.Rsave")
pp4<-OuchTest(tree,trait.tot)
save(pp4,file="ouch_Mamm.Rsave")

#Detaches the custom ouch
detach("package:ouch", unload=TRUE)
