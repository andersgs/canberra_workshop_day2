#Calculating FST with from SNP data using R
# and testing sampling strategies.
# by Anders Gonçalves da Silva (C) 2014
# email: andersgs@gmail.com
# 20 March 2014

#libraries needed
library(reshape2)
library(adegenet)
library(hierfstat)



possible_alleles = c('A','T','C','G')


sim_data = function(npop,nind,nloc,nal,mig,mut,eff_size){
  # use: sim_data(npop=2,nind=10,nloc=10,nal=2,mig=0.0001,mut=10^-9,eff_size=10^4)
  #
  #   The example above will generate a dataset where 2 populations were sampled
  #   (npop=2), 10 diploid individuals were taken (nind=10), and 10 loci
  #   were genotyped (nloc=10), each with a maximum of 2 alleles (nal=2) ---
  #   it is possible that monomorphic loci will be generated ---. The samples
  #   were taken assuming no selection, effective population size is 10000
  #   (that is 10000 diploid individuals, leadign to 20000 haploid genomes), that
  #   the two populations were exchaning migrants at a rate of 0.0001 individuals
  #   per generation, and were experiencing a mutation rate of 10^9/locus/generation
  #
  #   The function outputs a data.frame in long-format, with 5 columns:
  #     1. Population of origin (from 1 to npop)
  #     2. Individual (from 1 to nind*npop)
  #     3. Locus (from 1 to nloc)
  #     4. Allele 1 (A,T,C,or G)
  #     5. Allele 2 (A,T,C,or G)
  
  #create population, individual, and locus labels
  pop = rep(1:npop,each=nind*nloc)
  ind = rep(1:(npop*nind),each=nloc)
  loc = rep(seq(1:nloc),nind*npop)
  
  #figure out alleles for each locus. assumes each SNP type is equally likely.
  alleles = sapply(1:nloc,function(x) sample(possible_alleles,nal,replace=F))
  
  #sample allele frequencies for each locus from a Dirichlet distribution
  allele_frequencies = rdirichlet(n=nloc,alpha=rep(1,nal))
  
  #weigh the frequencies by the expected migration/mutation
  alpha = allele_frequencies*4*eff_size*(mig+mut)
  
  #create some storage for the alleles at each locus for each individual
  genotype1 <- character(npop*nind*nloc)
  genotype2 <- character(npop*nind*nloc)
  
  #sample alleles
  for(p in 1:npop){
    for(l in 1:nloc){
      tmp_alles = alleles[,l]
      al_freq = rdirichlet(n=1, alpha=alpha[l,])
      for(i in 1:nind){
        pos = (p-1)*nind*nloc+(i-1)*nloc+l
        gen = sample(tmp_alles,2,al_freq[1,],replace=T)
        genotype1[pos] = gen[1]
        genotype2[pos] = gen[2]
      }
    }
  }
  
  #return the simulated data
  return(data.frame(pop=pop,ind=ind,loc=loc,geno1=genotype1,geno2=genotype2))
}

test_df=sim_data(npop=10,nind=20,nloc=1000,nal=2,mig=0.0001,mu=10^-9,eff_size=10^4)

#reading the data into R
test_df = read.table("morning/first_session/test_data.txt",header=T)

#transforming it into a Genind object

test_genotype = paste(test_df[,4],test_df[,5],sep='.')

test_df$genotype<-test_genotype

test_df_wide=dcast(test_df,pop+ind~loc,value.var='genotype')

test_genind=df2genind(X=test_df_wide[,3:12],sep='\\.',ind.names=test_df_wide[,2],pop=test_df_wide[,1],loc.names=colnames(test_df_wide[,3:12]))

fstat(test_genind)

#simulate replicate data to see how things change depending on choice of sampling design
sim_reps = function(reps=100,npop=10,nind=20,nloc=100,nal=2,mig=0.0001,mut=10^-9,eff_size=10^4){
  #similar to sim_data, but repeats the process reps number of times
  # then calculates FST (Weir and Cockerham's 1984 version)
  
  #a results vector to hold the simulated FST values
  res = numeric(reps)
  
  #calculate the theoretically expected FST value
  exp = (1/(1+(4*eff_size*(mig+mut))))
  #store it in a nicely formatted number for printing
  exp_p = formatC(round(exp,2),2,format='f')
  
  #generate the datasets and then calculate FST
  for(i in 1:reps){
    test_df=sim_data(npop=npop,nind=nind,nloc=nloc,nal=nal,mig=mig,mu=mut,eff_size=eff_size)
    
    test_genotype = paste(test_df[,4],test_df[,5],sep='.')
    
    test_df$genotype<-test_genotype
    
    test_df_wide=dcast(test_df,pop+ind~loc,value.var='genotype')
    
    test_genind=df2genind(X=test_df_wide[,3:(nloc+2)],sep='\\.',ind.names=test_df_wide[,2],pop=test_df_wide[,1],loc.names=colnames(test_df_wide[,3:(nloc+2)]))
    
    #print some summary information
    obs = fstat(test_genind)[1,1]
    res[i] = obs
    bias = formatC(round(obs-exp,2),2,format='f')
    obs = formatC(round(obs,2),2,format='f')
    print(paste("Expected:",exp_p,"Observed:",obs, "Bias:",bias))
  }  
  return(res)
}

#low_pops, low_inds, many loci - FST = 0.1
test1 = sim_reps(reps=100,npop=2,nind=10,nloc=50,nal=2,mig=0.000225,mut=10^-9,eff_size=10^4)

#low_pops, high_inds, not that many loci - FST = 0.1
test2 = sim_reps(reps=100,npop=2,nind=50,nloc=10,nal=2,mig=0.000225,mu=10^-9,eff_size=10^4)

#high_pops, low_inds, few loci - FST = 0.1
test3 = sim_reps(reps=100,npop=10,nind=10,nloc=10,nal=2,mig=0.000225,mu=10^-9,eff_size=10^4)


#plot the distributions:
boxplot(list(test1,test2,test3))
abline(h=0.1,col='red')


####
#test_fstat = replicate(n=100,wc(sim.genot(size=20,nbal=2,nbloc=100,nbpop=10,N=20000,mig=0.0001,mut=10^-9,f=0))$FST)
