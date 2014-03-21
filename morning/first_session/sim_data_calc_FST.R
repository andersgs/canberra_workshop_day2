#libraries needed
library(reshape2)
library(adegenet)
library(hierfstat)

possible_alleles = c('A','T','C','G')


sim_data = function(npop,nind,nloc,nal,mig,mut,eff_size){
  #print((1/(1+(4*eff_size*(mig+mut)))))
  #nfun <- function(x,y) paste(x,y,sep='.')
  pop = rep(1:npop,each=nind*nloc)
  ind = rep(1:(npop*nind),each=nloc)
  loc = rep(seq(1:nloc),nind*npop)
  alleles = sapply(1:nloc,function(x) sample(possible_alleles,nal,replace=F))
  allele_frequencies = rdirichlet(n=nloc,alpha=rep(1,nal))
  alpha = allele_frequencies*4*eff_size*(mig+mut)
  genotype1 <- character(npop*nind*nloc)
  genotype2 <- character(npop*nind*nloc)
  for(p in 1:npop){
    for(l in 1:nloc){
      tmp_alles = alleles[,l]
      al_freq = rdirichlet(n=1, alpha=alpha[l,])
      #gen_freq = as.numeric(outer(al_freq,al_freq))
      #gens = outer(tmp_alles,tmp_alles,nfun)
      for(i in 1:nind){
        pos = (p-1)*nind*nloc+(i-1)*nloc+l
        gen = sample(tmp_alles,2,al_freq[1,],replace=T)
        genotype1[pos] = gen[1]
        genotype2[pos] = gen[2]
      }
    }
  }
  return(data.frame(pop=pop,ind=ind,loc=loc,geno1=genotype1,geno2=genotype2))
}

test_df=sim_data(npop=10,nind=20,nloc=1000,nal=2,mig=0.0001,mu=10^-9,eff_size=10^4)

#transforming it into a Genind object

test_genotype = paste(test_df[,4],test_df[,5],sep='.')

test_df$genotype<-test_genotype

test_df_wide=dcast(test_df,pop+ind~loc,value.var='geno1')

test_genind=df2genind(X=test_df_wide[,3:7],sep='\\.',ind.names=test_df_wide[,2],pop=test_df_wide[,1],loc.names=colnames(test_df_wide[,3:7]))

fstat(test_genind)

sim_reps = function(reps=100,npop=10,nind=20,nloc=100,nal=2,mig=0.0001,mut=10^-9,eff_size=10^4){
  res = numeric(reps)
  exp = formatC(round((1/(1+(4*eff_size*(mig+mut)))),2),2,format='f')
  for(i in 1:reps){
    test_df=sim_data(npop=npop,nind=nind,nloc=nloc,nal=nal,mig=mig,mu=mut,eff_size=eff_size)
    
    test_genotype = paste(test_df[,4],test_df[,5],sep='.')
    
    test_df$genotype<-test_genotype
    
    test_df_wide=dcast(test_df,pop+ind~loc,value.var='genotype')
    
    test_genind=df2genind(X=test_df_wide[,3:(nloc+2)],sep='\\.',ind.names=test_df_wide[,2],pop=test_df_wide[,1],loc.names=colnames(test_df_wide[,3:(nloc+2)]))
    
    obs = fstat(test_genind)[1,1]
    res[i] = obs
    obs = formatC(round(obs,2),2,format='f')
    bias = formatC(round(obs-exp,2),2,format='f')
    print(paste("Expected:",exp,"Observed:",obs, "Bias:",bias))
  }  
  return(res)
}

#low_pops, low_inds, many loci - FST = 0.1
test1 = sim_reps(reps=100,npop=2,nind=10,nloc=10000,nal=2,mig=0.000225,mut=10^-9,eff_size=10^4)

#low_pops, high_inds, not that many loci - FST = 0.1
test2 = sim_reps(reps=100,npop=2,nind=50,nloc=1000,nal=2,mig=0.000225,mu=10^-9,eff_size=10^4)

#high_pops, med_inds, few loci - FST = 0.1
test3 = sim_reps(reps=100,npop=10,nind=20,nloc=100,nal=2,mig=0.000225,mu=10^-9,eff_size=10^4)

test_fstat = replicate(n=100,wc(sim.genot(size=20,nbal=2,nbloc=100,nbpop=10,N=20000,mig=0.0001,mut=10^-9,f=0))$FST)
