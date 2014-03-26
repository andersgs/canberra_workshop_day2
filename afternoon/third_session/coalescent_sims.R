# Inferences on effective population size with the coalescent
# A maximum likelihood approach
# by Anders Gonçalves da Silva (C) 2014
# email: andersgs@gmail.com
# 20 March 2014

# simplest model: discrete-time coalescent (Kingman's n-coalescent)
# also known as the coalescent of n genes
# Hein et al. (2010) pg 21.

#Core idea:
# The time between coalescent events of a set of k genes is modeled as a geometric 
#  distribution:
#
# Tk ~ Geom(p)
# where p = (k 2)*1/(2*N) - where (k 2) is the number of combinations of 2 that 
#  are possible given k genes, and N is the number of diploid individuals in the
#  population, thus 2*N is the number of haploid genomes in the population
#
# So, the distribution changes with the number of genes that are still remain to
#  Coalesce

# Let us first see how the probability distribution changes with each k
n_genes = 5
n_diploid_ind = 10000

# a set of potential coalescent times, in number of generations
coal_times=seq(1,40000,100)
coal_times

# a function to calculate the probability of observing a certain coalescent time
#  given number of genes (k) and effective population size (N)
prob_coal_times = function(k,time, N){
  return(dgeom(time,prob=(choose(k,2)*(1/(2*N)))))
}

ps = sapply(2:n_genes,prob_coal_times,time=coal_times,N=n_diploid_ind)

#now let us plot the distributions
# notice how flat k=2 is...
matplot(x=coal_times,y=ps,type='l',lwd=2,xlab='Coalescent Time (T)',ylab="Pr(T|k,N)")
legend(10000,4e-4,c(2,3,4,5),col=c(1,2,3,4),lty=c(1,2,3,4),lwd=2,title='k genes to coalesce')

################################################################################
# writing a simple simulator of the n-coalescent
#
# the simulator will take N, the diploid populatons size, and n_genes, the
#  number of genes sampled. The n_genes are from a single locus.

#simulation function
sim_ncoal = function(N,n_genes){
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
    tau[k-1] = rgeom(n=1,prob=(choose(k,2)*(1/(2*N))))
    k = k - 1
  }
  #return the vector of times
  return(tau)
}

# to estimate N we need a likelihood function
coal_like = function(N,obs_times,n_genes){
  log_like = 0
  k = n_genes
  while(k > 1){
    log_like = log_like + dgeom(x=obs_times[k-1],prob=(choose(k,2)*(1/(2*N))),log=T)
    k = k - 1
  }
  return(log_like)
}

# simulate our coalescent times
n_genes = 25
n_diploid_ind = 10000

sim_times = sim_ncoal(n_diploid_ind,n_genes)

#distribution of time to MRCA
dist_total_height=replicate(n=10000,sum(sim_ncoal(n_diploid_ind,n_genes)))
hist(dist_total_height,breaks=100,xlab=c("Total tree hight (generations)"),ylab="Density",freq=F)

#add a line to the histogram that shows where our current simulation sits
abline(v=sum(sim_times),col='red',lwd=2)

#plotting the likelihood surface
#we are doing this numerically
#so, first we choose a number of potential population sizes
poss_N = seq(1000,50000,100)

#calculate the log-likelihood for each population size
loglike = sapply(poss_N,coal_like,obs_times=sim_times, n_genes=n_genes)

#plot the surface
#on the log-scale
plot(poss_N,loglike,type='l',xlab='Population size (N)',ylab="logLike(N|obsTimes,k)")

#on the natural scale
plot(poss_N,exp(loglike-max(loglike))/sum(exp(loglike-max(loglike))),type='l',xlab='Population size (N)',ylab="Like(N|obsTimes,k)")

#estimate N using maximum likelihood
ml_N = optim(par=100,coal_like,obs_times=sim_times,n_genes=n_genes,control=list(fnscale=-1),hessian=T,method='Brent',lower=100,upper=100000)
ml_N$par
#standard deviation
sqrt(solve(-1*ml_N$hessian))

###############################################################################
## What effect do extra loci have?

#this function wraps around the sim_ncoal function, and repeats it n_loci times
geom_coal_multi = function(n_diploid_ind,n_genes,n_loci){
  sims = replicate(n=n_loci,sim_ncoal(n_diploid_ind,n_genes))
  return(sims)
}

#to estimate N, we need a likelihood function
coal_like_multi = function(N,obs_times,n_genes,n_loci){
  log_like = 0
  k = n_genes
  
  #because we assume that each tree is an independent realisation of the coalescent process
  # we have to sum the log-likelihood of each coalescent event across trees
  while(k>1){
    log_like = log_like + sum(dgeom(x=obs_times[k-1,],prob=(choose(k,2)*(1/(2*N))),log=T))
    k = k -1
  }
  return(log_like)
}


# simulate some new data
n_loc = 100
n_genes = 50
n_diploid_ind = 10000

sim_times = geom_coal_multi(n_diploid_ind,n_genes,n_loc)
sim_times

#where do our 10 trees sit in the distribution of possible tree heights
dist_total_height=replicate(n=10000,sum(sim_ncoal(n_diploid_ind,n_genes)))
hist(dist_total_height,breaks=100,xlab=c("Total tree hight (generations)"),ylab="Density",freq=F)

#add a lines to the histogram that shows where our current simulations sit
abline(v=apply(sim_times,2,sum),col='red',lwd=2)

# ploting the likelihood surface
#choose potential N
poss_N = seq(1000,50000,100)

#calculate the log-likelihood for these values
loglike = sapply(poss_N,coal_like_multi,obs_times=sim_times, n_genes=n_genes)

#plot the surface
#on the log-scale
plot(poss_N,loglike,type='l',xlab='Population size (N)',ylab="logLike(N|obsTimes,k)")

#on the natural scale
plot(poss_N,exp(loglike-max(loglike))/sum(exp(loglike-max(loglike))),type='l',xlab='Population size (N)',ylab="Like(N|obsTimes,k)")

#estimate N using maximum likelihood
mlN = optim(par=100,coal_like_multi,obs_times=sim_times,n_genes=n_genes,control=list(fnscale=-1),hessian=T,method='Brent',lower=100,upper=100000)
mlN$par
#standard deviation
sqrt(solve(-1*mlN$hessian))
