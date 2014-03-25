# Inferences on population growth with the coalescent
# A skyline approach
# by Anders Gonçalves da Silva (C) 2014
# email: andersgs@gmail.com
# 20 March 2014

#load the ape library
# it has a number of relevant functions we will need
library(ape)

#population history

sim_ncoal_growth = function(N,n_genes,alpha,gen=0,max){
  #create some storage space of the times (measured in generations)
  # notice we only need (n_genes-1) slots 
  tau = numeric((n_genes-1))
  
  #for the sake of uniformity, let us re-name our n_genes variable k
  k = n_genes
  
  #now, we loop backwards until there is only a single ancestor left
  # we start off with n_genes, and we sample the time it will take for the 
  # n_genes to become (n_genes - 1). We sample that time from a geometric 
  # distribution with parameter p = (k 2)*(1/(2*N))
  while(k > 1){
  	if(sum(tau)>gen){
  		N = N*exp(-alpha*(sum(tau)-gen))
  		if((alpha < 0) && (N > max)){
  			N=max
  		} else if((alpha > 0) && (N < max)){
  				N=max
  			}
  		print(N)
  	}
    tau[k-1] = rgeom(n=1,prob=(choose(k,2)*(1/(2*N))))
    k = k - 1
  }
  #return the vector of times
  return(tau)
}

#simulate a genealogy
sim_times = sim_ncoal_growth(n_diploid_ind,n_genes,alpha=-5*10^-5,gen=1000,max=10^5)

#transform it into a coalescentIntervals object
x=list(lineages=seq(25,2),interval.length=rev(sim_times),interval.count=length(sim_times),total.depth=sum(sim_times))
attr(x,"class")<-"coalescentIntervals"

#plot the classic skyline plot
plot(skyline(x))

#plot the generalized skyline plot
plot(skyline(x,-1))

#let us try with 1000 loci
locs = matrix(0,nrow=1000,ncol=24)
for(i in 1:1000){
locs[i,]=sim_ncoal_growth(n_diploid_ind,n_genes,alpha=-6*10^-5,gen=1000,max=10^5)
}

#get mean interval lengths across genealogies
sim_times = apply(locs,2,mean)

#transform to a coalescentIntervals object
x=list(lineages=seq(25,2),interval.length=rev(sim_times),interval.count=length(sim_times),total.depth=sum(sim_times))
attr(x,"class")<-"coalescentIntervals"

#calculate the generalized skyline plot
plot(skyline(x,-1))

#use reversible-jump MCMC
mcmc.out <- mcmc.popsize(x,nstep=100000)
popsize <- extract.popsize(mcmc.out)

plot(popsize)
abline(h=c(10^4,10^5))
