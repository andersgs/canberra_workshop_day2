#simulate drift - demonstrate loss of heterozygosity but overall maintenance of allele frequencies
# by Anders Gonçalves da Silva (C) 2014
# email: andersgs@gmail.com
# 20 March 2014

#a function to simulate drift
sim_drift = function(initial_geno_count=c(4,8,4), generations=200,n_pop=50){  
  # 	use: sim_drift(initial_geno_count=c(4,8,4),generations=200,n_pop=50)
  #
  #		This will run 200 generations with 16 individuals are allowed to breed in each
  #		generation
  #		
  #		The experiment is replicated across 50 populations
  #		
  #		Initially, each population gets 4 'AA', 8 'Aa', and 4 'aa' individuals
  #
  #		The function outputs a list with allele frequency of 'A' for each population at
  # 	each timestep (sim_ps); the frequency of 'Aa' for each population at each step
  #		(sim_hs); mean frequency of 'A' across populations for each step (p_hat); and the
  # 	mean heterozygosity across populations for each time step (hs_hat).
  
  #auxiliary function to calculate allele frequency
  new_p = function(genos,n_ind){
    fAA = sum(genos=='AA')
    fAa = sum(genos=='Aa')
    return((2*fAA+fAa)/(2*n_ind))
  }
  
  #set some variables
  #number of individual to sample each generation
  n_ind = sum(initial_geno_count)
  
  #vector of initial population allele frequencies of 'A'
  cur_p = (2*initial_geno_count[1]+initial_geno_count[2])/(2*n_ind)
  cur_p = rep(cur_p,times=n_pop)
  
  #create some storage space and output matrices and vectors
  
  #store allele frequency in each generation
  ps = matrix(0,ncol=n_pop,nrow=generations)
  ps[1,] = cur_p
  
  #store the mean allele frequency across populations for each generation
  mean_p = numeric(generations)
  mean_p[1] = mean(cur_p)
  
  #store heterozygosity for across populations for each generation
  hs = matrix(0,ncol=n_pop,nrow=generations)
  hs[1,] = rep((initial_geno_count[2]/n_ind),times=n_pop)
  #store the mean heterozygosity across population for each generation
  mean_h = numeric(generations)
  mean_h[1] = mean(hs[1,])
  
  #store the genotypes for each individual in each population in a particular generation
  #this matrix gets re-written in each generation
  pop_genotypes = matrix(0,ncol=n_ind,nrow=n_pop)
  
  #initiate experiment
  #determine genotypic composition of first generation of all populations
  init_samp = c(rep('AA',initial_geno_count[1]),rep('Aa',initial_geno_count[2]),rep('aa',initial_geno_count[3]))
  #populate our experiment
  for(i in 1:n_pop){
    pop_genotypes[i,] = init_samp
  }
  
  #iterate over generations
  for(gen in 2:generations){
    # set the frequency of 'a'
    cur_q = 1-cur_p
    #create a data.frame with HW expected genotypic frequencies for each population
    # if one wanted to simulate inbreeding, it would be possible to add this here
    gf = data.frame(AA=cur_p^2,Aa=(2*cur_p*cur_q),aa=cur_q^2)
    
    #based on HW genotypic frequencies, generate a new set of n individuals for the next
    # generation
    for(p in 1:n_pop){
      pop_genotypes[p,] <- sample(c('AA','Aa','aa'),size=n_ind,replace=T,prob=gf[p,])
    }
    
    #calculate the new generations allele frequency for 'A'
    cur_p = ps[gen,] = apply(pop_genotypes,1,new_p,n_ind)
    
    #calculate mean allele frequency and observed heterzogysity 
    mean_p[gen] = mean(cur_p)
    hs[gen,] = apply(pop_genotypes,1,function(inds) sum(inds=='Aa')/n_ind)
    mean_h[gen] = mean(hs[gen,])
    
    #print progress to screen every 10 generations
    if(gen%%10 == 0){
      	plot(ps[1:gen,1],type='l',lwd=2,ylim=c(0,1),xlab='Generations',xlim=c(1,generations),ylab='Allele frequency')
      matlines(ps[1:gen,2:n_pop],lwd=2)
      title(paste("p =",mean_p[gen],"   h=",mean_h[gen]))
      Sys.sleep(0.02)
    }
  }
  # return list of results  
  return(list(sim_p=ps,sim_hs=hs,p_hat=mean_p,hs_hat=mean_h))
}

#run a simulation where we sample 16 individuals per population, for 200 populations and 200 generations. We start with 16 heterozygote individuals

drift1 = sim_drift(initial_geno_count=c(0,16,0),n_pop=200)
