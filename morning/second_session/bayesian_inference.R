#Bayesian inference - demonstrate the effect of changing prior and amount of data on 
#	the posterior distribution
# by Anders Gonçalves da Silva (C) 2014
# email: andersgs@gmail.com
# 20 March 2014

# The exercise consists of simulating flipping a biased coin, and trying to estimate
# the posterior distribution of theta, the probability of flipping a 'head'

#functions
#the log_likelihood function: calculates the log-likelihood of observing
# al=number of heads, be=number of tails in al+be trials
# given a value of theta.
#
# it assumes a Beta(alpha,beta) distribution for the counts of 'heads' and 'tails'

llike = function(theta,al,be){
  return(dbeta(theta,shape1=al,shape2=be,log=T))
}

#the prior function: calculate the probability of theta given prior knowledge on the
# expectation of number of 'heads' and 'tails' that one might observe.
#
# the probability of theta is assumed to derive from a Beta distribution with parameters
# alpha and beta, where alpha is the prior number of expected 'heads' and beta is the prior
# number of expected 'tails'. 
#
# a flat, or non-informative, prior is modelled as Beta(1,1). If there is strong belief that
# the coin is fair (i.e., theta = 0.5) then the prior could be Beta(10,10) or Beta(25,25).
# the higher the values of alpha and beta, the more information the prior will carry. We
# will see this below.
prior = function(theta,alpha,beta){
	return(dbeta(theta,shape1=alpha,shape2=beta,log=T))
}

#the post_prob function: the posterior is proportional to the likelihood times the prior. Or
# the sum of the log(likelihood) and the log(prior). Thus, the function returns the proportional
# posterior density of theta given the data and the prior
post_prob = function(theta,al,be,alpha,beta){
	return(llike(theta,al,be)+prior(theta,alpha,beta))
}

#to estimate the posterior probability distribution of theta, we will use MCMC.
# this problem is simple, and could easily be solved analytically. But, the goal of this
# exercise is to also demonstrate the inner functions of an MCMC using a Metropolis-Hastings
# algorithm.
mcmc = function(data,steps,alpha,beta,prop_window){
	#observed number of heads
	al = sum(data==1)
	#observed number of tails
	be = sum(data==0)
	
	#creating some place to store some values
	#Accepted simulated theta values
	Theta = numeric(steps)
	#Posterior probability of each of accepted thetas
	Post = numeric(steps)
	
	#initiate the chain
	cur_theta = runif(1)
	cur_post = post_prob(theta=cur_theta,al,be,alpha,beta)
	Theta[1] = cur_theta
	Post[1] = cur_post
	
	#run the chain
	for(i in 2:steps){
		#propose a new theta
		new_theta = runif(1,min=max(cur_theta-prop_window,0),max=min(cur_theta+prop_window,1))
		new_post = post_prob(theta=new_theta,al,be,alpha,beta)
		
		#if it is better than the current theta, accept it
		# if not, accept it with probability rho
		rho = runif(1)
		if(exp(new_post-cur_post) > rho){
			cur_theta = new_theta
			cur_post = new_post
		}
		
		#save the values of the chain
		Theta[i] = cur_theta
		Post[i] = cur_post
	}
	return(data.frame(theta=Theta,post=Post))
}


#generate some data
dat=replicate(n=100,rbinom(n=1,size=1,prob=0.8))
dat

#run the mcmc
#parameters
alpha=25
beta=25
chain=mcmc(dat,10000,alpha=alpha,beta=beta,prop_window=0.5)

################
#four plot panel
#plot prior
par(mar=c(5,6,4,2)+0.1,mfrow=c(2,2))
plot(seq(0,1,0.01),dbeta(seq(0,1,0.01),25,25),col='dark green',type='l',ylim=c(0,8),xlab=expression(theta),ylab="Density",lwd=2,cex.lab=2,cex.axis=1.5)
title('Prior')

#plot like
plot(seq(0,1,0.01),exp(sapply(seq(0,1,0.01),llike,sum(dat==1),sum(dat==0))),col='navy blue',xlab=expression(theta),ylab="Density",lwd=2,cex.lab=2,cex.axis=1.5,type='l')
title("Likelihood")

#plot post
plot(chain[,1],type='l',ylab=expression(theta),xlab="Iterations",lwd=2,cex.lab=2,cex.axis=1.5,ylim=c(0,1))
title("MCMC")

hist(chain[,1],freq=F,xlab=expression(theta),ylab="Density",lwd=2,cex.lab=2,cex.axis=1.5,xlim=c(0,1),main="Posterior")

lines(seq(0,1,0.01),sapply(seq(0,1,0.01),function(theta) dbeta(theta,sum(dat==1)+alpha,sum(dat==0)+beta)),col='red',lwd=2)

#################
#comparative plot
par(mar=c(5,6,4,2)+0.1,mfrow=c(1,1))
plot(seq(0,1,0.01),dbeta(seq(0,1,0.01),alpha,beta),col='dark green',type='l',ylim=c(0,10),xlab=expression(theta),ylab="Density",lwd=2,cex.lab=2,cex.axis=1.5)

lines(seq(0,1,0.01),exp(sapply(seq(0,1,0.01),llike,sum(dat==1),sum(dat==0))),col='navy blue',lwd=2,lty=2)

lines(seq(0,1,0.01),sapply(seq(0,1,0.01),function(theta) dbeta(theta,sum(dat==1)+alpha,sum(dat==0)+beta)),col='red',lwd=2,lty=3)
legend(x=0,y=6,legend=c('prior','likelihood','posterior'),lty=c(1,2,3),col=c('dark green','navy blue','red'),lwd=2)
par(mar=c(5, 4, 4, 2) + 0.1)

##################################################################
#flat prior

#generate some data
dat=replicate(n=10,rbinom(n=1,size=1,prob=0.8))

#run the mcmc
#parameters
alpha=1
beta=1
chain=mcmc(dat,10000,alpha=alpha,beta=beta,prop_window=0.5)

################
#four plot panel
#plot prior
par(mar=c(5,6,4,2)+0.1,mfrow=c(2,2))
plot(seq(0,1,0.01),dbeta(seq(0,1,0.01),alpha,beta),col='dark green',type='l',ylim=c(0,8),xlab=expression(theta),ylab="Density",lwd=2,cex.lab=2,cex.axis=1.5)
title('Prior')

#plot like
plot(seq(0,1,0.01),exp(sapply(seq(0,1,0.01),llike,sum(dat==1),sum(dat==0))),col='navy blue',xlab=expression(theta),ylab="Density",lwd=2,cex.lab=2,cex.axis=1.5,type='l')
title("Likelihood")

#plot post
plot(chain[,1],type='l',ylab=expression(theta),xlab="Iterations",lwd=2,cex.lab=2,cex.axis=1.5,ylim=c(0,1))
title("MCMC")

hist(chain[,1],freq=F,xlab=expression(theta),ylab="Density",lwd=2,cex.lab=2,cex.axis=1.5,xlim=c(0,1),main="Posterior")

lines(seq(0,1,0.01),sapply(seq(0,1,0.01),function(theta) dbeta(theta,sum(dat==1)+alpha,sum(dat==0)+beta)),col='red',lwd=2)

#################
#comparative plot
par(mar=c(5,6,4,2)+0.1,mfrow=c(1,1))
plot(seq(0,1,0.01),dbeta(seq(0,1,0.01),alpha,beta),col='dark green',type='l',ylim=c(0,10),xlab=expression(theta),ylab="Density",lwd=2,cex.lab=2,cex.axis=1.5)

lines(seq(0,1,0.01),exp(sapply(seq(0,1,0.01),llike,sum(dat==1),sum(dat==0))),col='navy blue',lwd=2,lty=2)

lines(seq(0,1,0.01),sapply(seq(0,1,0.01),function(theta) dbeta(theta,sum(dat==1)+alpha,sum(dat==0)+beta)),col='red',lwd=2,lty=3)
legend(x=0,y=6,legend=c('prior','likelihood','posterior'),lty=c(1,2,3),col=c('dark green','navy blue','red'),lwd=2)
par(mar=c(5, 4, 4, 2) + 0.1)