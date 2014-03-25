#creating structure input files using R.
# by Anders Gonçalves da Silva (C) 2014
# email: andersgs@gmail.com
# 20 March 2014

#load necessary libraries
library(reshape2)
library(gtools)

#some functions and definitions
possible_alleles = c('A','T','C','G')

#translate nucleotides to integers --- structure codes alleles as integers
nucl_switch=function(nucl) 
    {switch(nucl,
       "A" = 1,
       "T" = 2,
       "C" = 3,
       "G" = 4,
       "N" = 0
       )
}

# a function to simulate the data
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

#### generate some data
test_df=sim_data(npop=2,nind=20,nloc=500,nal=2,mig=0.0001,mu=10^-9,eff_size=10^4)

#melt then cast into structure mould
test_wide=melt(test_df,id.vars=c("pop","ind","loc"))
test_struc=dcast(test_wide,pop+ind+variable~loc)

#what the data looks like
final_struct = data.frame(test_struc[,c(2,1),],apply(test_struc[,c(4:503)],c(1,2),nucl_switch))

#write it to a file you can use
write.table(x=final_struct,file="infile",quote=F,row.names=F,col.names=F)


