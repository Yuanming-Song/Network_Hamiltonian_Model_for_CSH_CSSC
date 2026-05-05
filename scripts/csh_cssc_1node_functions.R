#Functions for estimation and evaluation of network Hamiltonian models for 1-node
#CSH/CSSC networks.  These are used, in particular, by the csh_cssc_1node_estimation.R function.
#
#Written by Carter T. Butts
#Modified by Yuanming Song
#Last updated 9/2/25 by CTB
#

#Load required libraries
library(sna)
library(ergm)
library(ergm.components)
library(ergm.multi)
library(parallel)


#Initial definitions (needed globally)
#resnames<-c("CSH","CSSC")  # Residue types


#Estimation function using a multi-compositional version of Yin and Butts pooling, with
#a version of SA adapted from Snijders (2002).  Composition is determined by size and the
#distribution of the "Resname" attribute within the input graphs.
#
#Arguments:
#  terms - a vector of ergm() model terms
#  train - a network.list object to be used for training (parameter estimation)
#  test - a network.list object to be used tor testing (model selection)
#  thin - number of thinning iterations to take during network simulation
#  dmax - maximum degree (max coordination number); this can be given per-species, though
#   at present only the max value is used
#  N1 - optionally, the number of draws for the initial variance estimate in the SA
#   algorithm; otherwise the heuristic of Snijders (2002) is used
#  N3 - the number of draws to use for final variance/covariance estimation
#  a1 - learning rate for the SA algorithm
#  phase.1.nr - logical; should a Newton-Rapheson step be used for initialization at
#    the end of phase 1?  This may speed up convergence if the initial estimate is close
#    to the MLE, but can also put the algorithm in a poor initial condition otherwise.
#  subphases - number of phases to use for SA (per Snijders (2002)).  Each phase is twice
#    as long as the one before it, so increasing this will greatly increase compute time;
#    however, more subphases can improve convergence, especially if the initial condition is
#    poor.
#  seed - optionally, a random number seed to use
#  theta.init - optionally, an initial coefficient vector; by default, an adjusted pooled
#    MPLE is used
#  verbose - logical; print progress messages?
#
#Return value:
#  A list containing the fitted model, together with estimated deviance on the test data.
#
#Details:
#  This function fits a pooled ERGM to the data in train, using a multi-composition extension of the 
#Yin-Butts pooling method.  Specifically, the input networks are split into compositional classes
#based on their "Resname" attributes, the mean sufficient statistics are computed for the networks
#in each class, and the group size weighted-mean is used as the vector of target statistics.
#A modified version of Snijders (2002) is used to perform estimation, where one graph at each update
#is produced from each compositional class, and the group size weighted mean deviation is used to
#produce the parameter update; where there is a single class, this simply becomes an SA-based
#implementation of standard Yin-Butts pooled maximum likelihood estimation.  An obvious benefit of
#this scheme is that it is quite fast, since one need draw only one graph per compositional class
#at each update, regardless of how many training graphs there are; after preprocessing is
#completed, larger data sizes do not increase compute time.  The SA algorithm (mostly) follows
#the Snijders (2002) scheme, and settings may be chosen accordingly.  Note that all models are
#fitted with a Krivitsky offset (so the input networks need not be the same size).
#  The initial condition for the pooled ERGM is obtained by a size-weighted average MPLE.  Given
#that some coefficients are not estimable within certain compositional settings, some basic
#sanity checking is done (e.g., removing non-finite estimates).  This is not fantastic, but it
#seems adequate for purpose.  The test deviance is computed using bridge sampling on the networks
#in test, given the fitted model.  Assuming this to be exchangable data from the same source as
#the training data, this can be used directly for model selection in the same manner as the AIC
#or AICc (with no additional complexity penalties), since the latter are analyical approximations to the 
#expected deviance on held-out data.
#
#Note that the edge coefficient reported is the base binding energy over -kT, and does not include
#either (1) the correction for the Krivitsky offset (-log(N)), nor the (2) correction for vibrational
#degrees of freedom (-1).  So, to actually simulate an ERGM with these coefficients, one must
#use theta[1]-1-log(N).  Using the resported edge coefficient without adjustments will result
#in heartache.  This adjustment is automatically used by the simulation and GOF methods asociated
#with models from estPool, so these coefficients can be used as-is for that purpose.  The reason
#for working with this modified form is that the network Hamiltonian parameters (phi) are
#then simply phi=-theta/(kT).  
#
ergmPool<-function(terms, train, test, thin=10*network.size(train[[1]])^2, dmax=network.size(train[[1]])-1, N1=NULL, N3=1000, a1=0.1, phase.1.nr=FALSE, subphases=4, seed=NULL, theta.init=NULL, verbose=TRUE){
  if(!is.null(seed))
    set.seed(seed)
  #First, classify the network types (here, exploiting the fact that there's only one covariate)
  if(verbose) cat("Classifying network types...\n")
  compvec<-sapply(train,function(z){table(z%v%"Resname")})
  compvec.test<-sapply(test,function(z){table(z%v%"Resname")})
  types<-unique(c(compvec,compvec.test))
  typeind<-match(compvec,types)
  typeind.test<-match(compvec.test,types)
  typecount<-tabulate(typeind,nbin=length(types))
  nv<-sapply(train,network.size)
  nv.test<-sapply(test,network.size)
  nt<-length(types)
  #Gather our statistics, while we may
  if(verbose) cat("Gathering statistics...\n")
  if(length(terms)==0)
    termstr<-""
  else
    termstr<-paste0("+",paste0(terms,collapse="+"))
  ostats<-summary(as.formula(paste0("train~edges",termstr))) #Training stats
  w<-(1/typecount[typeind])/nt  #Weighting vector
  tstats<-as.vector(crossprod(ostats,w))   #Target stats
  #Select one example from each type, finding the one closest to the mean
  if(verbose) cat("Choosing base nets...\n")
  selnets<-sapply(1:nt,function(z){sel<-which(typeind==z); tmean<-colMeans(ostats[sel,,drop=FALSE]); d2<-rowSums(sweep(ostats[sel,,drop=FALSE],2,tmean,"-")^2); sel[which.min(d2)]})
  basenets<-train[selnets]
  #Create our formula set (one per type)
  f<-list()
  for(i in 1:nt){
    f[[i]]<-as.formula(paste0("basenets[[",i,"]]~edges",termstr))
  }
  f.test<-as.formula(paste0("Networks(test) ~  N(~edges,offset=-log(sapply(test,network.size))) ", termstr))
  #Create our base simulation states - speeds things up considerably!
  sim.state<-list()
  for(i in 1:nt)
    sim.state[[i]]<-simulate(f[[i]], coef=NULL, nsim=0, control=control.simulate.formula(MCMC.interval=thin), constraints=~bd(maxout=max(dmax)), return.args="ergm_state")$object
  #Function to simulate stats across types
  simit<-function(theta,nsim=1){
    control<-control.simulate.formula(MCMC.interval=thin)
    control$MCMC.samplesize<-nsim
    #Start with the first one...
    co<-theta
    co[1]<-co[1]-1-log(nv[selnets[1]])
    stats<-ergm_MCMC_sample(sim.state[[1]], control=control, theta=co)$stats[[1]]
#    stats<-simulate(f[[1]], nsim=nsim, coef=c(theta[1]-1-log(nv[selnets[1]]),theta[-1]), control=control.simulate.formula(MCMC.interval=thin), output="stats")
    #Now increment, having established the size of the output
    for(i in 2:nt){
      co<-theta
      co[1]<-co[1]-1-log(nv[selnets[i]])
      stats<-stats+ergm_MCMC_sample(sim.state[[i]], control=control, theta=co)$stats[[1]]
#      stats<-stats+simulate(f[[i]], nsim=nsim, coef=c(theta[1]-1-log(nv[selnets[i]]),theta[-1]), control=control.simulate.formula(MCMC.interval=thin), output="stats")
    }
    #Return the normalized cross-type mean stats
    stats/nt
  }
  #Initialize with MPLE, averaged over types
  if(is.null(theta.init)){
    if(verbose) cat("Computing MPLEs...\n")
    theta<-vector()
    for(i in 1:nt)
      theta<-rbind(theta, coef(ergm(f[[i]], estimate="MPLE", constraints=~bd(maxout=max(dmax))), eval.loglik=FALSE))
    if(verbose) cat("Raw theta:\n")
    if(verbose) print(theta)
    theta[,1]<-theta[,1]+1+log(nv[selnets])
    theta[!is.finite(theta)]<-NA              #This can happen b/c of composition effects
    theta<-colMeans(theta,na.rm=TRUE)         #We only want to average well-defined coefs
    theta[!is.finite(theta)]<-0               #Ensure that we have well-behaved starting points
    if(verbose) cat("Adjusted theta:\n")
    if(verbose) print(theta)
  }else
    theta<-theta.init
  p<-length(theta)
  #Begin SA using an adapted version of Snijders (2002) - phase 1
  if(verbose) cat("SA phase 1...\n")
  if(is.null(N1))
    N1<-7+3*p
  ss<-simit(theta=theta,nsim=N1)
  mu<-colMeans(ss)
  D<-var(ss)
  diag(D)<-diag(D)+1e-1  #Mild regularization
  D0<-diag(D)
  if(phase.1.nr)
    theta<-as.vector(theta - pmax(pmin(a1*MASS::ginv(D)%*%(mu-tstats),5),-5))  #Initial NR step
  if(verbose) print(theta)
  #Now for phase 2
  if(verbose) cat("SA phase 2...\n")
  an<-a1
  for(k in 1:subphases){  #Walk through the subphases
    if(verbose) cat("\tSubphase",k,"\n")
    #Create the subphase step counts, if needed
    N2m<-2^(4*(k-1)/3)*(7+p)
    N2p<-N2m+200
    #Iterate
    flag<-FALSE
    i<-1
    z<-rep(0,p)
    mtheta<-rep(0,p)
    while((i<=N2m)||((!flag)&&(i<=N2p))){ #Stop if flag, or too many steps
      #Simulate the next draw
      if(verbose) print(theta)
      ss<-as.vector(simit(theta=theta,nsim=1))
      #Update
      zold<-z
      z<-ss-tstats
      if(sum(z*zold)<0)  #Early stopping condition
        flag<-TRUE
      theta<-theta-an*z/D0  #Local update
      mtheta<-mtheta+theta
      i<-i+1
    }
    theta<-mtheta/(i-1)  #Final value is average over the subphase
    if(verbose) print(theta)
    an<-an/2       #Decrement the learning rate for the next subphase
  }
  #Finally, phase 3 - estimate the variance/covariance matrix
  if(verbose) cat("SA phase 3...\n")
  sim<-simit(theta=theta,nsim=N3)
  scov<-var(sim)
  qcov<-MASS::ginv(scov)/length(train)  #Will be the inverse at the MLE; correct for N
  #Put the pieces together
  fit<-list(coef=theta, cov=qcov, se=diag(qcov)^0.5, ss.cov=scov, ss.mean=colMeans(sim), ss.target=tstats, terms=terms, dmax=dmax, train.sample=list(w=w, types=types, type.ind=typeind, type.count=typecount))
  #Now, we need to compute the deviance of the test data
  if(verbose) cat("Computing test deviance...\n")
  co<-theta
  co[1]<-co[1]-1  #Size is corrected for w/an offset in the formula
  dev<-ergm.bridge.0.llk(f.test, constraints=~bd(maxout=max(dmax)), coef=co, control=control.ergm.bridge(MCMC.interval=thin, MCMC.burnin=thin), llkonly=FALSE)
  fit$deviance.test<- -2*dev$llk
  fit$deviance.test.se<- 2*sqrt(dev$vcov.llr)
  #Return everything
  if(verbose) cat("Complete.\n")
  class(fit)<-"ergmPool"
  fit
}


#Coef method for ergmPool objects
coef.ergmPool<-function(object, ...){
  object$coef
}


#Print method for ergmPool objects
print.ergmPool<-function(x, ...){
  cat("Compositionally Pooled ERGM\n\nCoefs:\n")  
  print(x$coef)
}


#Summary method for ergmPool objects, and the associated print method
summary.ergmPool<-function(object, ...){
  class(object)<-c("summary.ergmPool",class(object))
  object
}

print.summary.ergmPool<-function(x, ...){
  cat("\nCompositionally Pooled ERGM\n\n")
  cat("Composition distribution:\n")
  cnt<-x$train.sample$type.count
  names(cnt)<-sapply(x$train.sample$types,paste0)
  print(cnt)
  cat("\nParameter Estimates:\n")  
  tab<-cbind(coef(x), x$se, coef(x)/x$se, 2*(1-pnorm(abs(coef(x)/x$se))))
  colnames(tab)<-c("Estimate", "Std.Err", "Z value", "Pr(>|z|)")
  printCoefmat(tab)
  cat("Test deviance",x$deviance.test,"(null deviance 0)\n\n")
}


#Special simulation command for models fit using estPool.  It requires an ergmPool object, 
#as well as one or more observed graphs to be used as templates (i.e., to supply vertex
#covariates and such).  
#
#Arguments:
#  object - ergmPool object from which draws should be taken
#  nsim - number of draws to take
#  seed - optionally, RNG seed to set
#  ... - additional arguments (currently ignored)
#  dat - list of base networks to use for simulation (one draw is taken for each); should be
#        a network.list, though a single network can be passed if desired.
#  thin - number of thinning iterations; by default, this is set to a multiple of
#         the number of networks times the square of the network size
#  verbose - logical; print progress messages?
#
#Return value:  A list of network.list objects, each of which corresponds to an element of dat;
#   if only a single input network was given, a single network.list will be returned.  Each
#   network.list will have length equal to nsim.
#
simulate.ergmPool<-function(object, nsim = 1, seed = NULL, ..., dat, thin=NULL, verbose = FALSE){
  #Set seed if desired
  if(!is.null(seed))
    set.seed(seed)
  #Set up the formula
  if(verbose) cat("Extracting the model...\n")
  if(!inherits(dat,"network.list")){            #Verify that the base nets are in the right form
    if(inherits(dat,"network")){
      dat<-list(dat)
      class(dat)<-"network.list"
    }else if(inherits(dat,"list")){
        class(dat)<-"network.list"
    }else
      stop("dat must be a network.list.\n")
  }
  if(length(object$terms)>0){
    f<-as.formula(paste0("Networks(dat) ~ N(~edges,offset=-log(sapply(dat,network.size))) + ", paste0(object$terms,collapse="+")))
  }else{
    f<-as.formula("Networks(dat) ~ N(~edges,offset=-log(sapply(dat,network.size)))")
  }
  #Simulate
  if(verbose) cat("Simulating draws for the selected base networks...\n")
  if(is.null(thin))                             #Ensure that thinning iterations are set
    thin<-length(dat)*network.size(dat[[1]])^2*10
  co<-coef(object)
  co[1]<-co[1]-1      #Need to adjust for bond vibrations
  sim<-ergm:::simulate.formula(f, nsim=nsim, coef=co, control=control.simulate.formula(MCMC.interval=thin), constraints=~bd(maxout=max(object$dmax)))
  if(nsim==1)
    sim<-list(sim)
  #Extract the results
  nets<-vector(mode="list",length=length(dat))
  for(i in 1:length(nets)){
    #Extract the draws from the ith input graph
    g.rep<-lapply(sim,function(z){z%s%which(z%v%".NetworkID"==i)})
    class(g.rep)<-"network.list"
    nets[[i]]<-g.rep
  }
  #Return the networks
  if(length(dat)==1)
    nets[[1]]
  else
    nets
}


#Model adequacy checking function (special case of gof()) for ergmPool, with support for checking
#versus held-out data.
#
#Arguments:
#  object - fitted model to check
#  dat - a vector of observed graphs against which to check
#  nets - sample size to take from dat
#  reps - number of simulation draws to take from fit
#  thin - thinning iterations 
#  verbose - logical; print progress messages?
#  ... - other parameters (ignored)
#
#Return value:
#  An object of class "gof.ergmPool", containing among other things  quantiles and z scores for
#  each of the observed networks versus the simulated networks.
#
#Details:
#  Currently, the statistics computed are the degree distribution, the ESP distribution, and
#  several component statistics (isolate count, dimer count, component count, and sum of squared
#  component sizes).  The quantile of each observed graph with respect to the simulated graphs on
#  each statistic is returned, as is the corresponding z-score.  In addition, the observed values,
#  and the means, SDs, and 95% simulation intervalus for the simulated samples are returned for
#  each target network.
#
#  The return value is also assigned a class ("gof.ergmPool"), facilitating downstream plot or
#  summary methods.
#
gof.ergmPool<-function(object, ..., dat, nets=25, reps=250, thin=nets*network.size(dat[[1]])^2*5, verbose=TRUE){
  #Draw sample and classify the network types (here, exploiting the fact that there's only one 
  #covariate)
  samp<-sample(1:length(dat),min(length(dat),nets)) #Subsample at random
  if(verbose) cat("Classifying network types...\n")
  compvec<-sapply(dat[samp],function(z){table(z%v%"Resname")})
  types<-unique(compvec)
  typeind<-match(compvec,types)
  typecount<-tabulate(typeind,nbin=length(types))
  if(verbose) cat("Extracting the model...\n")
  if(length(object$terms)>0){
    f<-as.formula(paste0("Networks(dat[samp]) ~ N(~edges,offset=-log(sapply(dat[samp],network.size))) + ", paste0(object$terms,collapse="+")))
  }else{
    f<-as.formula("Networks(dat[samp]) ~ N(~edges,offset=-log(sapply(dat[samp],network.size)))")
  }
  co<-coef(object)
  co[1]<-co[1]-1   #Need to correct for vibrational degrees of freedom
  names(co)[1]<-"N(1)~edges" #This is nuts, but ergm now needs it...
  if(verbose) cat("Simulating draws for the sampled networks...\n")
  print(co)
  sim<-simulate(f, nsim=reps, coef=co, control=control.simulate.formula(MCMC.interval=thin), constraints=~bd(maxout=max(object$dmax)))
  if(verbose) cat("Computing statistics...\n")
  maxd<-15
  maxe<-15
  o.deg<-matrix(0,nrow=nets,ncol=maxd+1)  #Observed
  o.esp<-matrix(0,nrow=nets,ncol=maxe+1)
  o.comp<-matrix(0,nrow=nets,ncol=4)
  m.deg<-matrix(0,nrow=nets,ncol=maxd+1)  #Mean sim
  m.esp<-matrix(0,nrow=nets,ncol=maxe+1)
  m.comp<-matrix(0,nrow=nets,ncol=4)
  sd.deg<-matrix(0,nrow=nets,ncol=maxd+1)  #SD sim
  sd.esp<-matrix(0,nrow=nets,ncol=maxe+1)
  sd.comp<-matrix(0,nrow=nets,ncol=4)
  q975.deg<-matrix(0,nrow=nets,ncol=maxd+1)  #97.5% sim quantile
  q975.esp<-matrix(0,nrow=nets,ncol=maxe+1)
  q975.comp<-matrix(0,nrow=nets,ncol=4)
  q025.deg<-matrix(0,nrow=nets,ncol=maxd+1)  #2.5% sim quantile
  q025.esp<-matrix(0,nrow=nets,ncol=maxe+1)
  q025.comp<-matrix(0,nrow=nets,ncol=4)
  q.deg<-matrix(0,nrow=nets,ncol=maxd+1)  #Quantile of obs in simulations
  q.esp<-matrix(0,nrow=nets,ncol=maxe+1)
  q.comp<-matrix(0,nrow=nets,ncol=4)
  z.deg<-matrix(0,nrow=nets,ncol=maxd+1)  #Z score of obs in simulations
  z.esp<-matrix(0,nrow=nets,ncol=maxe+1)
  z.comp<-matrix(0,nrow=nets,ncol=4)
  for(i in 1:nets){
    #Extract the networks
    g.rep<-lapply(sim,function(z){z%s%which(z%v%".NetworkID"==i)})
    class(g.rep)<-"network.list"
    #Degree
    obs<-summary(dat[[samp[i]]]~degree(0:maxd))
    if(i==1){
      snam<-names(obs)
      colnames(o.deg)<-snam
      colnames(m.deg)<-snam
      colnames(sd.deg)<-snam
      colnames(q975.deg)<-snam
      colnames(q025.deg)<-snam
      colnames(q.deg)<-snam
      colnames(z.deg)<-snam
    }
    sstats<-summary(g.rep~degree(0:maxd))
    o.deg[i,]<-obs
    m.deg[i,]<-colMeans(sstats)
    sd.deg[i,]<-apply(sstats,2,sd)
    q975.deg[i,]<-apply(sstats,2,quantile,0.975)
    q025.deg[i,]<-apply(sstats,2,quantile,0.025)
    q.deg[i,]<-colMeans(sweep(sstats,2,obs,"<="))
    z.deg[i,]<-(obs-m.deg[i,])/(sd.deg[i,]+1e-8)
    #ESP
    obs<-summary(dat[[samp[i]]]~esp(0:maxe))
    if(i==1){
      snam<-names(obs)
      colnames(o.esp)<-snam
      colnames(m.esp)<-snam
      colnames(sd.esp)<-snam
      colnames(q975.esp)<-snam
      colnames(q025.esp)<-snam
      colnames(q.esp)<-snam
      colnames(z.esp)<-snam
    }
    sstats<-summary(g.rep~esp(0:maxe))
    o.esp[i,]<-obs
    m.esp[i,]<-colMeans(sstats)
    sd.esp[i,]<-apply(sstats,2,sd)
    q975.esp[i,]<-apply(sstats,2,quantile,0.975)
    q025.esp[i,]<-apply(sstats,2,quantile,0.025)
    q.esp[i,]<-colMeans(sweep(sstats,2,obs,"<="))
    z.esp[i,]<-(obs-m.esp[i,])/(sd.esp[i,]+1e-8)
    #Components
    obs<-summary(dat[[samp[i]]]~isolates+dimers+components+compsizesum(pow=2))
    if(i==1){
      snam<-names(obs)
      colnames(o.comp)<-snam
      colnames(m.comp)<-snam
      colnames(sd.comp)<-snam
      colnames(q975.comp)<-snam
      colnames(q025.comp)<-snam
      colnames(q.comp)<-snam
      colnames(z.comp)<-snam
    }
    sstats<-summary(g.rep~isolates+dimers+components+compsizesum(pow=2))
    o.comp[i,]<-obs
    m.comp[i,]<-colMeans(sstats)
    sd.comp[i,]<-apply(sstats,2,sd)
    q975.comp[i,]<-apply(sstats,2,quantile,0.975)
    q025.comp[i,]<-apply(sstats,2,quantile,0.025)
    q.comp[i,]<-colMeans(sweep(sstats,2,obs,"<="))
    z.comp[i,]<-(obs-m.comp[i,])/(sd.comp[i,]+1e-8)
  }
  out<-list(degree.obs=o.deg, degree.mean=m.deg, degree.sd=sd.deg, degree.q025=q025.deg, degree.q975=q975.deg, degree.q=q.deg, degree.z=z.deg, comp.obs=o.comp, comp.mean=m.comp, comp.sd=sd.comp, comp.q025=q025.comp, comp.q975=q975.comp, comp.q=q.comp, comp.z=z.comp, esp.obs=o.esp, esp.mean=m.esp, esp.sd=sd.esp, esp.q025=q025.esp, esp.q975=q975.esp, esp.q=q.esp, esp.z=z.esp, types=types, type.count=typecount, type.ind=typeind)
  class(out)<-"gof.ergmPool"
  out
}


#Summary method for gof.ergmPool objects.  Crude, but perhaps helpful.
#
#Currently, we compute the following summaries:
#  - 95% simulation interval coverage by statistic, and the p-value for a one-sided binomial
#    test (over test networks) that the coverage probability is 0.95 (vs. the hypothesis that
#    the coverage probability is <0.95).
#  - The mean Z-scores, by statistic, across networks.
#  - An approximate "Chebychev-Fisher" test for excessive deviations from the predictive mean.
#    This is very crude, but intended to be a very robust test.  We start by observing that the
#    chance of observing a Z-value whose absolute value is greater than |z| is <= 1/z^2 (by
#    Chebyshev's inequality), no matter how the test statistic is distributed.  We can thus 
#    take this to be a very conservative approximate p-value for a two-sided test of the hypothesis
#    that the expected statistic for the observed graph was equal to the expectation under the
#    simulated model.  Given that we have multiple independent observed graphs, we then use 
#    Fisher's method to combine the p-values (i.e., computing X^2=sum_i -2 log p_i, and taking
#    X^2 to be chi-squared under the null with 2m degrees of freedom).  Since our p-values are
#    greater than or equal to (and usually much greater than) the "true" p-values we would otherwise
#    obtain, our X^2 approximation is biased down, and the test will be very conservative (i.e.,
#    loath to reject).  But we also have a lot of graphs, so in fact this still ends up having
#    quite a lot of power (as one can demonstrate to oneself with numerical experiments).  Why
#    not just look at the simulated quantiles, instead?  One can, and this is certainly a more
#    refined test, but it is not really helpful to ask if there is any detectable difference 
#    from the simulated distribution, because we know that the model is by intent an approximation.
#    More subtly, the variable discreteness of the quantiles in practice makes those sorts of
#    analyses quite tricky in practice.  It may thus be somewhat more useful to specifically look
#    at deviations from the mean, as a heuristic, and using a fairly conservative test.  It is 
#    intendedly heuristic, but then, no one should use p-values as anything other than heuristics
#    anyway.
#  
summary.gof.ergmPool<-function(object, ...){
 m<-NROW(object$degree.obs)  #Number of graphs/samples
 #Create summaries for quantile coverage
 cov95.deg<-apply((object$degree.obs>=object$degree.q025)&(object$degree.obs<=object$degree.q975),2,mean)
 cov95.deg.p<-pbinom((1-cov95.deg)*m,size=m,prob=0.05,lower.tail=FALSE)
 cov95.comp<-apply((object$comp.obs>=object$comp.q025)&(object$comp.obs<=object$comp.q975),2,mean)
 cov95.comp.p<-pbinom((1-cov95.comp)*m,size=m,prob=0.05,lower.tail=FALSE)
 cov95.esp<-apply((object$esp.obs>=object$esp.q025)&(object$esp.obs<=object$esp.q975),2,mean)
 cov95.esp.p<-pbinom((1-cov95.esp)*m,size=m,prob=0.05,lower.tail=FALSE)
 #Create summaries for Z-scores
 meanZ.deg<-colMeans(object$degree.z)
 meanZ.comp<-colMeans(object$comp.z)
 meanZ.esp<-colMeans(object$esp.z)
 chebfish.deg.x2<-colSums(-2*log(pmin(1/object$degree.z^2,1)))
 chebfish.deg.p<-1-pchisq(chebfish.deg.x2,df=2*m)
 chebfish.comp.x2<-colSums(-2*log(pmin(1/object$comp.z^2,1)))
 chebfish.comp.p<-1-pchisq(chebfish.comp.x2,df=2*m)
 chebfish.esp.x2<-colSums(-2*log(pmin(1/object$esp.z^2,1)))
 chebfish.esp.p<-1-pchisq(chebfish.esp.x2,df=2*m)
 #Save the output
 out<-list(coverage.95.deg=cov95.deg, coverage.95.deg.pval=cov95.deg.p, coverage.95.comp=cov95.comp, coverage.95.comp.pval=cov95.comp.p, coverage.95.esp=cov95.esp, coverage.95.esp.pval=cov95.esp.p, mean.z.deg=meanZ.deg, mean.z.comp=meanZ.comp, mean.z.esp=meanZ.esp, chebfish.deg.chisq=chebfish.deg.x2, chebfish.deg.pval=chebfish.deg.p, chebfish.comp.chisq=chebfish.comp.x2, chebfish.comp.pval=chebfish.comp.p, chebfish.esp.chisq=chebfish.esp.x2, chebfish.esp.pval=chebfish.esp.p)
 object$summary<-out
 class(object)<-"summary.gof.ergmPool"
 object
}

#Simple print method for summary.mygof
print.summary.gof.ergmPool<-function(x, digits=3, ...){
  cat("Multi-compositional Model Adequacy Summary\n\n")
  #Degree
  tab<-cbind(x$summary$mean.z.deg, x$summary$coverage.95.deg, x$summary$coverage.95.deg.pval, x$summary$chebfish.deg.chisq, x$summary$chebfish.deg.pval)
  colnames(tab)<-c("Mean Z", "95% Cov", "Pr(Cov<=C)", "Chisq", "Pr(X^2>=x^2)")
  cat("Degree:\n")
  print(tab,digits=digits)
  #ESP
  tab<-cbind(x$summary$mean.z.esp, x$summary$coverage.95.esp, x$summary$coverage.95.esp.pval, x$summary$chebfish.esp.chisq, x$summary$chebfish.esp.pval)
  colnames(tab)<-c("Mean Z", "95% Cov", "Pr(Cov<=C)", "Chisq", "Pr(X^2>=x^2)")
  cat("\nEdgewise Shared Partners:\n")
  print(tab,digits=digits)
  #Components
  tab<-cbind(x$summary$mean.z.comp, x$summary$coverage.95.comp, x$summary$coverage.95.comp.pval, x$summary$chebfish.comp.chisq, x$summary$chebfish.comp.pval)
  colnames(tab)<-c("Mean Z", "95% Cov", "Pr(Cov<=C)", "Chisq", "Pr(X^2>=x^2)")
  cat("\nComponents:\n")
  print(tab,digits=digits)
}


#Simple plot method for gof.ergmPool
plot.gof.ergmPool<-function(x, ...){
  #Set up plotting parameters
  op<-par(no.readonly=TRUE)
  on.exit(par(op))
  #Plot distributions by type, exploiting homogeneity of simulations by type
  nt<-length(x$types)
  par(mfrow=n2mfrow(nt))
  oask<-devAskNewPage(ask=TRUE)
  on.exit(devAskNewPage(ask=oask),add=TRUE)
  for(i in 1:nt){  #Degree
    main<-paste0("Degree Distribution, ",x$types[i])
    ylab<-"Vertices"
    xlab<-""
    sel<-x$type.ind==i      #Get the graphs of type i
    bso<-x$degree.obs[sel,,drop=FALSE]                #Obs values
    q025<-colMeans(x$degree.q025[sel,,drop=FALSE])    #Estimated 0.025 quantile
    q975<-colMeans(x$degree.q975[sel,,drop=FALSE])    #Estimated 0.975 quantile
    m<-colMeans(x$degree.mean[sel,,drop=FALSE])       #Estimated mean
    vals<-1:NCOL(x$degree.mean)                       #Index values
    xnam<-colnames(x$degree.obs)                      #Stat names
    xl<-range(vals)                                   #Plot limits
    yl<-range(c(as.vector(bso),q025,q975))
    plot(1,1,type="n", axes=FALSE, xlim=xl, ylim=yl, ylab=ylab, xlab=xlab, main=main, font.lab=2)
    axis(2,font=2)
    axis(1,font=2,las=3,at=vals,labels=xnam)
    polygon(c(vals,rev(vals)), c(q025,rev(q975)), border=NULL, col=rgb(1,0,0,alpha=0.05))
    for(j in 1:NROW(bso)){     #Plot the observed values
      points(vals,bso[j,], col=1)
    }
    lines(vals, colMeans(bso), lty=2, col=1, lwd=2)
    lines(vals, m, lwd=3, col=2)
    lines(vals, q025, lwd=2, col=2)
    lines(vals, q975, lwd=2, col=2)
  }
  par(mfrow=n2mfrow(nt))
  for(i in 1:nt){  #Edgewise shared partners
    main<-paste0("ESP Distribution, ",x$types[i])
    ylab<-"Edges"
    xlab<-""
    sel<-x$type.ind==i      #Get the graphs of type i
    bso<-x$esp.obs[sel,,drop=FALSE]                #Obs values
    q025<-colMeans(x$esp.q025[sel,,drop=FALSE])    #Estimated 0.025 quantile
    q975<-colMeans(x$esp.q975[sel,,drop=FALSE])    #Estimated 0.975 quantile
    m<-colMeans(x$esp.mean[sel,,drop=FALSE])       #Estimated mean
    vals<-1:NCOL(x$esp.mean)                       #Index values
    xnam<-colnames(x$esp.obs)                      #Stat names
    xl<-range(vals)                                   #Plot limits
    yl<-range(c(as.vector(bso),q025,q975))
    plot(1,1,type="n", axes=FALSE, xlim=xl, ylim=yl, ylab=ylab, xlab=xlab, main=main, font.lab=2)
    axis(2,font=2)
    axis(1,font=2,las=3,at=vals,labels=xnam)
    polygon(c(vals,rev(vals)), c(q025,rev(q975)), border=NULL, col=rgb(1,0,0,alpha=0.05))
    for(j in 1:NROW(bso)){     #Plot the observed values
      points(vals,bso[j,], col=1)
    }
    lines(vals, colMeans(bso), lty=2, col=1, lwd=2)
    lines(vals, m, lwd=3, col=2)
    lines(vals, q025, lwd=2, col=2)
    lines(vals, q975, lwd=2, col=2)
  }
  par(mfrow=n2mfrow(nt))
  for(i in 1:nt){  #Components
    main<-paste0("Component Statistics, ",x$types[i])
    ylab<-"Value"
    xlab<-""
    sel<-x$type.ind==i      #Get the graphs of type i
    bso<-x$comp.obs[sel,,drop=FALSE]                #Obs values
    q025<-colMeans(x$comp.q025[sel,,drop=FALSE])    #Estimated 0.025 quantile
    q975<-colMeans(x$comp.q975[sel,,drop=FALSE])    #Estimated 0.975 quantile
    m<-colMeans(x$comp.mean[sel,,drop=FALSE])       #Estimated mean
    vals<-1:NCOL(x$comp.mean)                       #Index values
    xnam<-colnames(x$comp.obs)                      #Stat names
    xl<-range(vals)                                   #Plot limits
    yl<-range(c(as.vector(bso),q025,q975))
    yvals<-round(seq(from=yl[1],to=yl[2],length=5))
    plot(1,1,type="n", axes=FALSE, xlim=xl, ylim=log1p(yl), ylab=ylab, xlab=xlab, main=main, font.lab=2)
    axis(2,font=2,at=log1p(yvals),labels=yvals)
    axis(1,font=2,las=3,at=vals,labels=xnam)
    polygon(c(vals,rev(vals)), log1p(c(q025,rev(q975))), border=NULL, col=rgb(1,0,0,alpha=0.05))
    for(j in 1:NROW(bso)){     #Plot the observed values
      points(vals,log1p(bso[j,]), col=1)
    }
    lines(vals, log1p(colMeans(bso)), lty=2, col=1, lwd=2)
    lines(vals, log1p(m), lwd=3, col=2)
    lines(vals, log1p(q025), lwd=2, col=2)
    lines(vals, log1p(q975), lwd=2, col=2)
  }
}
