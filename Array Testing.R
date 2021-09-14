############################
#####  Load Libraries  #####
############################
library(truncnorm)
library(MCMCpack)
library(BayesLogit)
library(R.utils)
library(coga)
#library(pracma)

#######################
#####  Functions  #####
#######################

###Binary Link Function for Disease Status
g.link <-function(x){
  res = exp(x)/(1 + exp(x))
  return(res)
}

###Biomarker Distribution Function
f <-function(zeta,theta0,theta1,Y){ 
  if(Y==0){
    res = dgamma(zeta,theta0[1],theta0[2])
  }
  if(Y==1){
    res = dgamma(zeta,theta1[1],theta1[2])
  }
  return(res)
}

###Sample Biomarker Measurement Error Variance
sigma.w.samp<-function(){
  res = rinvgamma(1,J/2 + a.sigma.w,sum((W.tilde.g - W)^2)/2 + b.sigma.w)
  return(res)
}

###Sample the Individual Biomarker Levels
zeta.tilde.samp <-function(zeta.ind){
  zeta.tilde.i.prop = rnorm(1,zeta.tilde.g[zeta.ind],sd.zeta.tilde[zeta.ind])
  zeta.tilde.prop = zeta.tilde.g
  zeta.tilde.prop[zeta.ind] = zeta.tilde.i.prop
  W.tilde.prop = W.tilde.g
  for(z.i in pools.containing.ind[[zeta.ind]]){
    W.tilde.prop[z.i] = mean(zeta.tilde.prop[mathcalP[[z.i]]])
  }
  den.g = sum(dnorm(W[pools.containing.ind[[zeta.ind]]], mean = W.tilde.g[pools.containing.ind[[zeta.ind]]],sqrt(sigma.w.g),log = TRUE))
  den.g = den.g + log(f(zeta.tilde.g[zeta.ind],theta.0.g,theta.1.g,Y.tilde.g[zeta.ind]))
  den.g = den.g + log(dnorm(zeta.tilde.i.prop,zeta.tilde.g[zeta.ind],sd.zeta.tilde[zeta.ind]))
  den.prop = sum(dnorm(W[pools.containing.ind[[zeta.ind]]], W.tilde.prop[pools.containing.ind[[zeta.ind]]],sqrt(sigma.w.g),log = TRUE))
  den.prop = den.prop + log(f(zeta.tilde.prop[zeta.ind],theta.0.g,theta.1.g,Y.tilde.g[zeta.ind]))
  den.prop = den.prop + log(dnorm(zeta.tilde.g[zeta.ind],zeta.tilde.i.prop, sd.zeta.tilde[zeta.ind]))
  A = min(1,exp(den.prop-den.g))
  if(is.na(A)){
    A = 0
  }
  r = rbinom(1,1,A)
  zeta.tilde.i.res = r*zeta.tilde.i.prop + (1-r)*zeta.tilde.g[zeta.ind]
  W.tilde.res = r*W.tilde.prop + (1-r)*W.tilde.g
  if(is.na(r)){
    print(den.g)
    print(den.prop)
    print(zeta.tilde.i.prop)
  }
  return(list(zeta.out = zeta.tilde.i.res,r.out = r,W.tilde.out = W.tilde.res))
}

###Sample the Latent true individual statuses
Y.tilde.samp <-function(){
  p0 = f(zeta.tilde.g,theta.0.g,theta.1.g,0)
  p0 = p0*(1-g.link(X%*%beta.g))
  p1 = f(zeta.tilde.g,theta.0.g,theta.1.g,1)
  p1 = p1*g.link(X%*%beta.g)
  y.out = rbinom(N,1,(p1/(p0+p1)))
  return(y.out)
}

###Sample the Poyla Gamma random variable for the population regression model
omega.samp<-function(){
  res = rpg(N,1,X%*%beta.g)
  return(res)
}

###Sample the population regression coefficients
beta.samp<-function(){
  sigma = solve(t(X)%*%Omega.g%*%X + Phi.inv.g)
  mean = sigma%*%t(X)%*%h.g
  res = mvrnorm(1,mean,sigma)
  return(res)
}

###Sample the biomarker distribution parameters for the negative biomarker distribution
#Note: We should be able to change from normal to gamma by redefing the f function and leaving this function unchanged
theta.0.samp <-function(theta.i){
  theta.prop <-rtruncnorm(1,theta.0.g[theta.i],sd =sd.theta.0[theta.i],a = theta.0.prior.lb[theta.i], b = theta.0.prior.ub[theta.i])
  theta.0.prop <- theta.0.g
  theta.0.prop[theta.i] = theta.prop
  den.g = sum(log(f(zeta.tilde.g[Y.tilde.g ==0],theta.0.g,theta.1.g,0)))
  #den.g = den.g + log(dtruncnorm(theta.0.g[theta.i],mean = theta.0.prior.mean[theta.i], sd = theta.0.prior.sd[theta.i],a = theta.0.prior.lb[theta.i], b = theta.0.prior.ub[theta.i]))
  den.g = den.g + log(dgamma(theta.0.g[theta.i],theta.0.prior.alpha[theta.i], theta.0.prior.beta[theta.i]))
  den.g = den.g + log(dtruncnorm(theta.prop,theta.0.g[theta.i],sd.theta.0[theta.i],a = theta.0.prior.lb[theta.i], b = theta.0.prior.ub[theta.i]))
  den.p = sum(log(f(zeta.tilde.g[Y.tilde.g ==0],theta.0.prop,theta.1.g,0)))
  #den.p = den.p + log(dtruncnorm(theta.0.prop[theta.i],mean = theta.0.prior.mean[theta.i], sd = theta.0.prior.sd[theta.i],a = theta.0.prior.lb[theta.i], b = theta.0.prior.ub[theta.i]))
  den.p = den.p + log(dgamma(theta.0.prop[theta.i],theta.0.prior.alpha[theta.i],theta.0.prior.beta[theta.i]))
  den.p = den.p + log(dtruncnorm(theta.0.g[theta.i],theta.0.prop[theta.i],sd.theta.0[theta.i],a = theta.0.prior.lb[theta.i], b = theta.0.prior.ub[theta.i]))
  A = min(1,exp(den.p-den.g))
  if(is.na(A)){
    A = 0
  }
  r = rbinom(1,1,A)
  theta.out = r*theta.prop + (1-r)*theta.0.g[theta.i]
  return(list(theta.out = theta.out, r.out = r))
}

###Sample the biomarker distribution parameters for the positive biomarker distribution
#Note: We should be able to change from normal to gamma by redefing the f function and leaving this function unchanged
theta.1.samp <-function(theta.i){
  theta.prop <-rtruncnorm(1,theta.1.g[theta.i],sd =sd.theta.1[theta.i],a = theta.1.prior.lb[theta.i], b = theta.1.prior.ub[theta.i])
  theta.1.prop <- theta.1.g
  theta.1.prop[theta.i] = theta.prop
  den.g = sum(log(f(zeta.tilde.g[Y.tilde.g ==1],theta.0.g,theta.1.g,1)))
  den.g = den.g + log(dgamma(theta.1.g[theta.i],theta.1.prior.alpha[theta.i],theta.1.prior.beta[theta.i]))
  den.g = den.g + log(dtruncnorm(theta.prop,theta.1.g[theta.i],sd.theta.1[theta.i],a = theta.1.prior.lb[theta.i], b = theta.1.prior.ub[theta.i]))
  den.p = sum(log(f(zeta.tilde.g[Y.tilde.g ==1],theta.0.g,theta.1.prop,1)))
  den.p = den.p + log(dgamma(theta.1.prop[theta.i],theta.1.prior.alpha[theta.i],theta.1.prior.beta[theta.i]))
  den.p = den.p + log(dtruncnorm(theta.1.g[theta.i],theta.1.prop[theta.i],sd.theta.1[theta.i],a = theta.1.prior.lb[theta.i], b = theta.1.prior.ub[theta.i]))
  A = min(1,exp(den.p-den.g))
  if(is.na(A)){
    A = 0
  }
  r = rbinom(1,1,A)
  theta.out = r*theta.prop + (1-r)*theta.1.g[theta.i]
  return(list(theta.out = theta.out, r.out = r))
}

w.con.den.S_p<-function(x,w,c.s,theta.0.mean,theta.0.sd,sig.w){
  res = pcoga(w,shape = rep(theta.0.mean,c.s),rate = rep(theta.0.sd,c.s)*c.s)
  res = res*dnorm(x,0,sqrt(sig.w))
  return(res)
}

w.con.den.S_e<-function(x,w,c.s,theta.0.mean,theta.0.sd, theta.1.mean,theta.1.sd,sig.w){
  res = pcoga(w,shape = c(rep(theta.0.mean,c.s-1),theta.1.mean), rate = c(rep(theta.0.sd,c.s-1),theta.1.sd)*c.s)
  res = res*dnorm(x,0,sqrt(sig.w))
  return(res)
}

S_p <- function(w,c.s,theta.0.mean,theta.0.sd,sig.w){
  res = integrate(w.con.den.S_p,lower = -Inf, upper = Inf,c.s = c.s,w = w,theta.0.mean=theta.0.mean,theta.0.sd = theta.0.sd, sig.w = sig.w )
  return(res$value)
}

###Sensitivity of a pool of size c with threshold of w, conditional on the pool containing 1 positive
S_e1 <- function(w,c.s,theta.0.mean,theta.0.sd,theta.1.mean,theta.1.sd,sig.w){
  res = 1-integrate(w.con.den.S_e,lower = -Inf, upper = Inf,c.s = c.s,w = w,theta.0.mean=theta.0.mean,theta.0.sd = theta.0.sd,theta.1.mean = theta.1.mean, theta.1.sd = theta.1.sd, sig.w = sig.w )$value
  return(res)
}

###Youden.index to be optimized to find optimal threshold
Youden.index <-function(w,fixed.par){
  c.s = fixed.par[1]
  theta.0.mean = fixed.par[2]
  theta.0.sd = fixed.par[3]
  theta.1.mean = fixed.par[4]
  theta.1.sd = fixed.par[5]
  sig.w = fixed.par[6]
  res = S_p(w,c.s,theta.0.mean,theta.0.sd,sig.w)
  res = res + S_e1(w,c.s,theta.0.mean,theta.0.sd,theta.1.mean,theta.1.sd,sig.w)
  res = res -1
  return(-res)
}

#############################
#####  Data Generation  #####
#############################
N = 900
X = cbind(1,rnorm(N),rbinom(N,1,0.5))
Q = dim(X)[2]
beta.true = c(-5,2,1)
eta.true = X%*%beta.true
p.true = g.link(eta.true)
Y.tilde.true = rbinom(N,1,p.true)
theta.0.true = c(2.5,0.5)
theta.1.true = c(80,2)
zeta.tilde.true = rgamma(N,theta.0.true[1]*(1-Y.tilde.true) + theta.1.true[1]*Y.tilde.true, theta.0.true[2]*(1-Y.tilde.true) + theta.1.true[2]*Y.tilde.true)
plot.seq = seq(0,50,length.out = 1000)
plot(plot.seq,dgamma(plot.seq,theta.0.true[1],theta.0.true[2]),type = 'l',ylim = c(0,0.5))
points(plot.seq,dgamma(plot.seq,theta.1.true[1],theta.1.true[2]),type = 'l',col = 'red')

###Create Pools for Array Testing  #####
c.size = 3
n.array = N/(c.size^2)
J = n.array*2*c.size
mathcalP = list()
ct = 1
for(k in 1:n.array){
  mat.k = matrix(seq(((k-1)*(c.size^2)+1),by = 1, length.out = (c.size^2)),c.size,c.size)
  for(r in 1:c.size){
    mathcalP[[ct]] = mat.k[r,]
    ct = ct + 1
    mathcalP[[ct]] = mat.k[,r]
    ct = ct+1
  }
}

pools.containing.ind = list()
for(i in 1:N){
  pool.list = c()
  for(j in 1:J){
    if(is.element(i,mathcalP[[j]])){
      pool.list = c(pool.list,j)
    }
  }
  pools.containing.ind[[i]] = pool.list
  print(i)
}

W.tilde.true = rep(NA,J)
for(j in 1:J){
  W.tilde.true[j] = mean(zeta.tilde.true[mathcalP[[j]]])
}
sigma.w.true = 0.005^2
W = rnorm(J,W.tilde.true,sqrt(sigma.w.true))
plot(W)
plot(W.tilde.true)

#####################################
#####  Specify Hyperparameters  #####
#####################################
N.cc = 100
N.cc.r = 1
Y.0 = rgamma(N.cc,theta.0.true[1],theta.0.true[2]) + rnorm(N.cc,0,sqrt(sigma.w.true))
mu.0 = mean(Y.0)
var.0 = var(Y.0) 

Y.1 = rgamma(N.cc,theta.1.true[1],theta.1.true[2]) + rnorm(N.cc,0,sqrt(sigma.w.true))
mu.1 = mean(Y.1)
var.1 = var(Y.1) 

theta.0.prior.mean = c(mu.0^2/var.0,mu.0/var.0)
theta.0.prior.var= 0.1*theta.0.prior.mean
theta.0.prior.alpha = theta.0.prior.mean^2/theta.0.prior.var
theta.0.prior.beta = theta.0.prior.mean/theta.0.prior.var
theta.1.prior.mean = c(mu.1^2/var.1,mu.1/var.1)
theta.1.prior.var= 0.1*theta.1.prior.mean
theta.1.prior.alpha = theta.1.prior.mean^2/theta.1.prior.var
theta.1.prior.beta = theta.1.prior.mean/theta.1.prior.var
theta.0.prior.lb = c(0,0)
theta.0.prior.ub = c(Inf,Inf)
theta.1.prior.lb = c(0,0)
theta.1.prior.ub = c(Inf,Inf)
a.sigma.w = 3
b.sigma.w = 3

############################
#####  The First MCMC  #####
############################
G = 30000
burn = 20000
zeta.tilde.g = rep(0.5,N)
W.tilde.g = rep(NA,J)
for(j in 1:J){
  W.tilde.g[j] = mean(zeta.tilde.g[mathcalP[[j]]])
}
sigma.w.g = sigma.w.true
Y.tilde.g = rep(0,N)
theta.0.g = theta.0.true
theta.1.g = theta.1.true
beta.g = beta.true
Phi.inv.g = solve(100*diag(Q))
omega.g = rep(1,N)
Omega.g = diag(omega.g)
h.g = (Y.tilde.g -1/2)

###Tuning
acc.zeta.tilde = rep(0,N)
sd.zeta.tilde = rep(0.1,N)
acc.theta.0 = rep(0,length(theta.0.g))
sd.theta.0 = rep(0.5,length(theta.0.g))
acc.theta.1 = rep(0,length(theta.1.g))
sd.theta.1 = rep(0.5,length(theta.1.g))

###Storage 
W.tilde.array = array(NA,c(J,G))
sigma.w.array = array(NA,c(1,G))
zeta.tilde.array = array(NA,c(N,G))
Y.tilde.array = array(NA,c(N,G))
theta.0.array = array(NA,c(length(theta.0.g),G))
theta.1.array = array(NA,c(length(theta.1.g),G))
beta.array = array(NA,c(Q,G))
omega.array = array(NA,c(N,G))

for(g in 1:G){
  ###Sample zeta.tilde
  for(i in 1:N){
    zeta.vals = zeta.tilde.samp(i)
    zeta.tilde.g[i] = zeta.vals$zeta.out
    acc.zeta.tilde[i] = acc.zeta.tilde[i] + zeta.vals$r.out
    W.tilde.g = zeta.vals$W.tilde.out
    if(is.na(zeta.tilde.g[i])){
      break
    }
  }
  if(is.na(zeta.tilde.g[i])){
    break
  }
  ###Sample sigma.w
  sigma.w.g = sigma.w.samp()
  ###Sample Y.tilde
  Y.tilde.g = Y.tilde.samp()
  h.g = (Y.tilde.g -1/2)
  ###Sample omega
  omega.g = omega.samp()
  Omega.g = diag(omega.g)
  h.g = (Y.tilde.g -1/2)
  ###Sample beta
  beta.g = beta.samp()
  ###Sample the biomarker distribution parameters
  for(i in 1:length(theta.0.g)){
    theta.vals = theta.0.samp(i)
    theta.0.g[i] = theta.vals$theta.out
    acc.theta.0[i] = acc.theta.0[i] + theta.vals$r.out
  }
  for(i in 1:length(theta.1.g)){
    theta.vals = theta.1.samp(i)
    theta.1.g[i] = theta.vals$theta.out
    acc.theta.1[i] = acc.theta.1[i] + theta.vals$r.out
  }
  W.tilde.array[,g] = W.tilde.g
  sigma.w.array[,g] = sigma.w.g
  zeta.tilde.array[,g] = zeta.tilde.g
  Y.tilde.array[,g] = Y.tilde.g
  theta.0.array[,g] = theta.0.g
  theta.1.array[,g] = theta.1.g
  beta.array[,g] = beta.g
  omega.array[,g] = omega.g
  if(g%%100 == 0 & g < burn){
    acc.zeta.tilde = acc.zeta.tilde/100
    ind = acc.zeta.tilde < 0.3
    sd.zeta.tilde[ind] = sd.zeta.tilde[ind]*0.9
    ind = acc.zeta.tilde > 0.7
    sd.zeta.tilde[ind] = sd.zeta.tilde[ind]*1.1
    acc.zeta.tilde = rep(0,N)
    acc.theta.0 = acc.theta.0/100
    ind = acc.theta.0 < 0.3
    sd.theta.0[ind] = sd.theta.0[ind]*0.9
    ind = acc.theta.0 > 0.7
    sd.theta.0[ind] = sd.theta.0[ind]*1.1
    acc.theta.0 = rep(0,length(theta.0.g))
    acc.theta.1 = acc.theta.1/100
    ind = acc.theta.1 < 0.3
    sd.theta.1[ind] = sd.theta.1[ind]*0.9
    ind = acc.theta.1 > 0.7
    sd.theta.1[ind] = sd.theta.1[ind]*1.1
    acc.theta.1 = rep(0,length(theta.0.g))
  }
  print(g)
}

###Generate Posterior Threshold Sample
ind.threshold.array = array(NA,c(1,G))
pool.threshold.array = array(NA,c(1,G))
#for(g in burn:G){
theta.0.g = apply(theta.0.array[,burn:G],1,mean)
theta.1.g = apply(theta.1.array[,burn:G],1,mean)
sigma.w.g = mean(sigma.w.array[,burn:G])
fixed.par.input = c(c.s = 1,theta.0.mean = theta.0.g[1], theta.0.sd = theta.0.g[2], theta.1.mean = theta.1.g[1], theta.1.sd = theta.1.g[2], sig.w = sigma.w.g)
pool.optim <- optimize(f = Youden.index, lower = 0, upper = 30, fixed.par = fixed.par.input)
ind.threshold.array[,burn:G] = pool.optim$minimum
fixed.par.input = c(c = c.size,theta.0.mean = theta.0.g[1], theta.0.sd = theta.0.g[2], theta.1.mean = theta.1.g[1], theta.1.sd = theta.1.g[2], sig.w = sigma.w.g)
pool.optim <- optimize(f = Youden.index, lower = 0, upper = 30, fixed.par = fixed.par.input)
pool.threshold.array[,burn:G] = pool.optim$minimum
#}

###Calculate True Optimal Thresholds for Individual Testing
fixed.par.input = c(c = 1,theta.0.mean = theta.0.true[1], theta.0.sd = theta.0.true[2], theta.1.mean = theta.1.true[1], theta.1.sd = theta.1.true[2], sig.w = sigma.w.true)
pool.optim <- optimize(f = Youden.index, lower = 0, upper = 30, fixed.par = fixed.par.input)
plot.seq = seq(0.001,30,length.out = 100)
plot.out = c()
for(i in 1:length(plot.seq)){
  plot.out = c(plot.out, Youden.index(plot.seq[i],fixed.par = fixed.par.input  ))
}
plot(plot.seq,plot.out,typ = 'l')
ind.threshold.true = pool.optim$minimum
abline(v = ind.threshold.true)
plot(plot.seq,dcoga(plot.seq, rep(theta.0.true[1],1),rep(theta.0.true[2],1)*1),type = 'l',ylim = c(0,1))
points(plot.seq,dcoga(plot.seq,c(rep(theta.0.true[1],1-1),theta.1.true[1]), c(rep(theta.0.true[2],1-1),theta.1.true[2])*1),type = 'l',col = 'red')
abline(v = ind.threshold.true)

###Calculate True Optimal Thresholds for Pool of size c.size
fixed.par.input = c(c = c.size,theta.0.mean = theta.0.true[1], theta.0.sd = theta.0.true[2], theta.1.mean = theta.1.true[1], theta.1.sd = theta.1.true[2], sig.w = sigma.w.true)
pool.optim <- optimize(f = Youden.index, lower = 0, upper = 30, fixed.par = fixed.par.input)
plot.seq = seq(-0.01,30,length.out = 100)
plot.out = c()
for(i in 1:length(plot.seq)){
  plot.out = c(plot.out, Youden.index(plot.seq[i],fixed.par = fixed.par.input  ))
}
plot(plot.seq,plot.out,typ = 'l')
pool.threshold.true = pool.optim$minimum
abline(v = pool.threshold.true)
plot(plot.seq,dcoga(plot.seq, rep(theta.0.true[1],c.size),rep(theta.0.true[2],c.size)*c.size),type = 'l')
points(plot.seq,dcoga(plot.seq,c(rep(theta.0.true[1],c.size-1),theta.1.true[1]), c(rep(theta.0.true[2],c.size-1),theta.1.true[2])*c.size),type = 'l',col = 'red')
abline(v = pool.threshold.true)

########################################################
#####  Calculate Point Estimators from First MCMC  #####
########################################################

beta.hat = apply(beta.array[,burn:G],1,mean)
beta.sd = apply(beta.array[,burn:G],1,sd)
beta.lower = apply(beta.array[,burn:G],1,quantile,prob = 0.025)
beta.upper = apply(beta.array[,burn:G],1,quantile,prob = 0.975)
beta.cov = beta.lower <= beta.true & beta.true <= beta.upper

sigma.w.hat = mean(sigma.w.array[,burn:G])
sigma.w.sd = sd(sigma.w.array[,burn:G])
sigma.w.lower = quantile(sigma.w.array[,burn:G],prob = 0.025)
sigma.w.upper = quantile(sigma.w.array[,burn:G],prob = 0.975)
sigma.w.cov = sigma.w.lower <= sigma.w.true & sigma.w.true <= sigma.w.upper

zeta.tilde.hat = apply(zeta.tilde.array[,burn:G],1,mean)
zeta.tilde.sd = apply(zeta.tilde.array[,burn:G],1,sd)
zeta.tilde.lower = apply(zeta.tilde.array[,burn:G],1,quantile,prob = 0.025)
zeta.tilde.upper = apply(zeta.tilde.array[,burn:G],1,quantile,prob = 0.975)
zeta.tilde.cov = zeta.tilde.lower <= zeta.tilde.true & zeta.tilde.true <= zeta.tilde.upper

W.tilde.hat = apply(W.tilde.array[,burn:G],1,mean)
W.tilde.sd = apply(W.tilde.array[,burn:G],1,sd)
W.tilde.lower = apply(W.tilde.array[,burn:G],1,quantile,prob = 0.025)
W.tilde.upper = apply(W.tilde.array[,burn:G],1,quantile,prob = 0.975)
W.tilde.cov = W.tilde.lower <= W.tilde.true & W.tilde.true <= W.tilde.upper

Y.tilde.hat = apply(Y.tilde.array[,burn:G],1,mean)
Y.tilde.sd = apply(Y.tilde.array[,burn:G],1,sd)
Y.tilde.lower = apply(Y.tilde.array[,burn:G],1,quantile,prob = 0.025)
Y.tilde.upper = apply(Y.tilde.array[,burn:G],1,quantile,prob = 0.975)
Y.tilde.cov = Y.tilde.lower <= Y.tilde.true & Y.tilde.true <= Y.tilde.upper

theta.0.hat = apply(theta.0.array[,burn:G],1,mean)
theta.0.sd = apply(theta.0.array[,burn:G],1,sd)
theta.0.lower = apply(theta.0.array[,burn:G],1,quantile,prob = 0.025)
theta.0.upper = apply(theta.0.array[,burn:G],1,quantile,prob = 0.975)
theta.0.cov = theta.0.lower <= theta.0.true & theta.0.true <= theta.0.upper

theta.1.hat = apply(theta.1.array[,burn:G],1,mean)
theta.1.sd = apply(theta.1.array[,burn:G],1,sd)
theta.1.lower = apply(theta.1.array[,burn:G],1,quantile,prob = 0.025)
theta.1.upper = apply(theta.1.array[,burn:G],1,quantile,prob = 0.975)
theta.1.cov = theta.1.lower <= theta.1.true & theta.1.true <= theta.1.upper

ind.threshold.hat = mean(ind.threshold.array[,burn:G])
ind.threshold.sd = sd(ind.threshold.array[,burn:G])
ind.threshold.lower = quantile(ind.threshold.array[,burn:G],prob = 0.025)
ind.threshold.upper = quantile(ind.threshold.array[,burn:G],prob = 0.975)
ind.threshold.cov = ind.threshold.lower <= ind.threshold.true & ind.threshold.true <= ind.threshold.upper

pool.threshold.hat = mean(pool.threshold.array[,burn:G])
pool.threshold.sd = sd(pool.threshold.array[,burn:G])
pool.threshold.lower = quantile(pool.threshold.array[,burn:G],prob = 0.025)
pool.threshold.upper = quantile(pool.threshold.array[,burn:G],prob = 0.975)
pool.threshold.cov = pool.threshold.lower <= pool.threshold.true & pool.threshold.true <= pool.threshold.upper

##############################################
#####  Assess Convergence of First MCMC  #####
##############################################
plot(beta.array[1,])
abline(h = beta.true[1])
plot(beta.array[2,])
abline(h = beta.true[2])
plot(beta.array[3,])
abline(h = beta.true[3])

plot(sigma.w.array[1,])
abline(h = sigma.w.true)

plot(theta.0.array[1,])
abline(h = theta.0.true[1])
plot(theta.0.array[2,])
abline(h = theta.0.true[2])
plot(theta.0.array[1,]/theta.0.array[2,])
abline(h = theta.0.true[1]/theta.0.true[2])
plot(theta.0.array[1,]/theta.0.array[2,]^2)
abline(h = theta.0.true[1]/theta.0.true[2]^2)
plot(theta.1.array[1,])
abline(h = theta.1.true[1])
plot(theta.1.array[2,])
abline(h = theta.1.true[2])
plot(theta.1.array[1,]/theta.1.array[2,])
abline(h = theta.1.true[1]/theta.1.true[2])
plot(theta.1.array[1,]/theta.1.array[2,]^2)
abline(h = theta.1.true[1]/theta.1.true[2]^2)


plot(zeta.tilde.array[1,])
abline(h=zeta.tilde.true[1])
plot(zeta.tilde.array[100,])
abline(h=zeta.tilde.true[100])
plot(zeta.tilde.array[200,])
abline(h=zeta.tilde.true[200])
plot(zeta.tilde.array[300,])
abline(h=zeta.tilde.true[300])

plot(W.tilde.array[1,])
abline(h = W.tilde.true[1])
plot(W.tilde.array[2,])
abline(h = W.tilde.true[2])
plot(W.tilde.array[3,])
abline(h = W.tilde.true[3])
plot(W.tilde.array[4,])
abline(h = W.tilde.true[4])

plot(ind.threshold.array[1,])
abline(h = ind.threshold.hat)
plot(pool.threshold.array[1,])
abline(h = pool.threshold.hat)

####################################################################
#####  Indentify Individuals to Retest Based on the First MCMC  ####
####################################################################

###Classify Pools
pool.status = W > pool.threshold.hat
pool.status.true = rep(0,J)
for(j in 1:J){
  ind.in.pool = mathcalP[[j]]
  if(sum(Y.tilde.true[ind.in.pool])>0){
    pool.status.true[j] = 1
  }
}

###Identify Individuals to Retest
mathcalP.retest = list()
retest.running.list = c()
ct = 1
for(i in 1:N){
  pools.i = pools.containing.ind[[i]]
  if(sum(pool.status[pools.i])>=2){
    mathcalP.retest[[ct]] = i
    retest.running.list = c(retest.running.list,i)
    ct = ct +1
  }
}
individuals.restested = unique(unlist(mathcalP.retest))
N.retest = length(individuals.restested)
J.retest = N.retest

###Retest the individuals
Y.tilde.true.retest = Y.tilde.true[individuals.restested]
zeta.tilde.true.retest = zeta.tilde.true[individuals.restested]
W.tilde.true.retest = zeta.tilde.true.retest
W.retest = rnorm(J.retest,W.tilde.true.retest,sqrt(sigma.w.true))

###Classify Individuals Based on First Round of Testing
ind.status.orginal = rep(0,N)
for(i in 1:N){
  if(is.element(i,individuals.restested)){
    ind = which(individuals.restested == i)
    if(W.retest[ind]>ind.threshold.hat){
      ind.status.orginal[i] = 1
    }
  }
}

####################################
#####  Set up the Second MCMC  #####
####################################

###Store Values from the First Round of Testing
W.original = W
mathcalP.original = mathcalP
J.original = J
pools.containing.ind.original = pools.containing.ind
W.tilde.array.original = W.tilde.array
sigma.w.array.original = sigma.w.array
zeta.tilde.array.original = zeta.tilde.array
Y.tilde.array.original = Y.tilde.array
theta.0.array.original = theta.0.array 
theta.1.array.original = theta.1.array
beta.array.original = beta.array
omega.array.original = omega.array
ind.threshold.array.original = ind.threshold.array
pool.threshold.array.original = pool.threshold.array

###Redefine variables to incorporate the new data
W = c(W.original,W.retest)
J = J.original + J.retest
ct = J.original +1
if(J.retest>0){
  for(j in 1:J.retest){
    mathcalP[[ct]] = mathcalP.retest[[j]]
    ct = ct + 1
  }
}

###Determine which pools contain each individual
pools.containing.ind.retest = list()
for(i in 1:N){
  pool.list = c()
  for(j in 1:J){
    if(is.element(i,mathcalP[[j]])){
      pool.list = c(pool.list,j)
    }
  }
  pools.containing.ind.retest[[i]] = pool.list
  print(i)
}
pools.containing.ind = pools.containing.ind.retest

###Initialize the MCMC
W.tilde.g = rep(NA,J)
for(j in 1:J){
  W.tilde.g[j] = mean(zeta.tilde.g[mathcalP[[j]]])
}

###Tuning
acc.zeta.tilde = rep(0,N)
sd.zeta.tilde = rep(0.1,N)
acc.theta.0 = rep(0,length(theta.0.g))
sd.theta.0 = rep(0.5,length(theta.0.g))
acc.theta.1 = rep(0,length(theta.1.g))
sd.theta.1 = rep(0.5,length(theta.1.g))

###Storage 
W.tilde.array = array(NA,c(J,G))
sigma.w.array = array(NA,c(1,G))
zeta.tilde.array = array(NA,c(N,G))
Y.tilde.array = array(NA,c(N,G))
theta.0.array = array(NA,c(length(theta.0.g),G))
theta.1.array = array(NA,c(length(theta.1.g),G))
tau.array = array(NA,c(Q,G))
v.array = array(NA,c(Q,G))
phi.array = array(NA,c(Q,G))
beta.array = array(NA,c(Q,G))
omega.array = array(NA,c(N,G))

for(g in 1:G){
  ###Sample zeta.tilde
  for(i in 1:N){
    zeta.vals = zeta.tilde.samp(i)
    zeta.tilde.g[i] = zeta.vals$zeta.out
    acc.zeta.tilde[i] = acc.zeta.tilde[i] + zeta.vals$r.out
    W.tilde.g = zeta.vals$W.tilde.out
    if(is.na(zeta.tilde.g[i])){
      break
    }
  }
  if(is.na(zeta.tilde.g[i])){
    break
  }
  ###Sample sigma.w
  sigma.w.g = sigma.w.samp()
  ###Sample Y.tilde
  Y.tilde.g = Y.tilde.samp()
  h.g = (Y.tilde.g -1/2)
  ###Sample omega
  omega.g = omega.samp()
  Omega.g = diag(omega.g)
  h.g = (Y.tilde.g -1/2)
  ###Sample beta
  beta.g = beta.samp()
  ###Sample the biomarker distribution parameters
  for(i in 1:length(theta.0.g)){
    theta.vals = theta.0.samp(i)
    theta.0.g[i] = theta.vals$theta.out
    acc.theta.0[i] = acc.theta.0[i] + theta.vals$r.out
  }
  for(i in 1:length(theta.1.g)){
    theta.vals = theta.1.samp(i)
    theta.1.g[i] = theta.vals$theta.out
    acc.theta.1[i] = acc.theta.1[i] + theta.vals$r.out
  }
  W.tilde.array[,g] = W.tilde.g
  sigma.w.array[,g] = sigma.w.g
  zeta.tilde.array[,g] = zeta.tilde.g
  Y.tilde.array[,g] = Y.tilde.g
  theta.0.array[,g] = theta.0.g
  theta.1.array[,g] = theta.1.g
  beta.array[,g] = beta.g
  omega.array[,g] = omega.g
  if(g%%100 == 0 & g < burn){
    acc.zeta.tilde = acc.zeta.tilde/100
    ind = acc.zeta.tilde < 0.3
    sd.zeta.tilde[ind] = sd.zeta.tilde[ind]*0.9
    ind = acc.zeta.tilde > 0.7
    sd.zeta.tilde[ind] = sd.zeta.tilde[ind]*1.1
    acc.zeta.tilde = rep(0,N)
    acc.theta.0 = acc.theta.0/100
    ind = acc.theta.0 < 0.3
    sd.theta.0[ind] = sd.theta.0[ind]*0.9
    ind = acc.theta.0 > 0.7
    sd.theta.0[ind] = sd.theta.0[ind]*1.1
    acc.theta.0 = rep(0,length(theta.0.g))
    acc.theta.1 = acc.theta.1/100
    ind = acc.theta.1 < 0.3
    sd.theta.1[ind] = sd.theta.1[ind]*0.9
    ind = acc.theta.1 > 0.7
    sd.theta.1[ind] = sd.theta.1[ind]*1.1
    acc.theta.1 = rep(0,length(theta.0.g))
  }
  print(g)
}

###Generate Posterior Threshold Sample
ind.threshold.array = array(NA,c(1,G))
pool.threshold.array = array(NA,c(1,G))
#for(g in burn:G){
theta.0.g = apply(theta.0.array[,burn:G],1,mean)
theta.1.g = apply(theta.1.array[,burn:G],1,mean)
sigma.w.g = mean(sigma.w.array[,burn:G])
fixed.par.input = c(c.s = 1,theta.0.mean = theta.0.g[1], theta.0.sd = theta.0.g[2], theta.1.mean = theta.1.g[1], theta.1.sd = theta.1.g[2], sig.w = sigma.w.g)
pool.optim <- optimize(f = Youden.index, lower = 0, upper = 30, fixed.par = fixed.par.input)
ind.threshold.array[,burn:G] = pool.optim$minimum
fixed.par.input = c(c = c.size,theta.0.mean = theta.0.g[1], theta.0.sd = theta.0.g[2], theta.1.mean = theta.1.g[1], theta.1.sd = theta.1.g[2], sig.w = sigma.w.g)
pool.optim <- optimize(f = Youden.index, lower = 0, upper = 30, fixed.par = fixed.par.input)
pool.threshold.array[,burn:G] = pool.optim$minimum
#}

#########################################################
#####  Calculate Point Estimators from Second MCMC  #####
#########################################################

beta.hat.2 = apply(beta.array[,burn:G],1,mean)
beta.sd.2 = apply(beta.array[,burn:G],1,sd)
beta.lower.2 = apply(beta.array[,burn:G],1,quantile,prob = 0.025)
beta.upper.2 = apply(beta.array[,burn:G],1,quantile,prob = 0.975)
beta.cov.2 = beta.lower.2 <= beta.true & beta.true <= beta.upper.2

sigma.w.hat.2 = mean(sigma.w.array[,burn:G])
sigma.w.sd.2 = sd(sigma.w.array[,burn:G])
sigma.w.lower.2 = quantile(sigma.w.array[,burn:G],prob = 0.025)
sigma.w.upper.2 = quantile(sigma.w.array[,burn:G],prob = 0.975)
sigma.w.cov.2 = sigma.w.lower.2 <= sigma.w.true & sigma.w.true <= sigma.w.upper.2

zeta.tilde.hat.2 = apply(zeta.tilde.array[,burn:G],1,mean)
zeta.tilde.sd.2 = apply(zeta.tilde.array[,burn:G],1,sd)
zeta.tilde.lower.2 = apply(zeta.tilde.array[,burn:G],1,quantile,prob = 0.025)
zeta.tilde.upper.2 = apply(zeta.tilde.array[,burn:G],1,quantile,prob = 0.975)
zeta.tilde.cov.2 = zeta.tilde.lower.2 <= zeta.tilde.true & zeta.tilde.true <= zeta.tilde.upper.2

W.tilde.hat.2 = apply(W.tilde.array[,burn:G],1,mean)
W.tilde.sd.2 = apply(W.tilde.array[,burn:G],1,sd)
W.tilde.lower.2 = apply(W.tilde.array[,burn:G],1,quantile,prob = 0.025)
W.tilde.upper.2 = apply(W.tilde.array[,burn:G],1,quantile,prob = 0.975)
W.tilde.true.2 = c(W.tilde.true,W.tilde.true.retest)
W.tilde.cov.2 = W.tilde.lower.2 <= W.tilde.true.2 & W.tilde.true.2 <= W.tilde.upper.2

Y.tilde.hat.2 = apply(Y.tilde.array[,burn:G],1,mean)
Y.tilde.sd.2 = apply(Y.tilde.array[,burn:G],1,sd)
Y.tilde.lower.2 = apply(Y.tilde.array[,burn:G],1,quantile,prob = 0.025)
Y.tilde.upper.2 = apply(Y.tilde.array[,burn:G],1,quantile,prob = 0.975)
Y.tilde.cov.2 = Y.tilde.lower.2 <= Y.tilde.true & Y.tilde.true <= Y.tilde.upper.2

theta.0.hat.2 = apply(theta.0.array[,burn:G],1,mean)
theta.0.sd.2 = apply(theta.0.array[,burn:G],1,sd)
theta.0.lower.2 = apply(theta.0.array[,burn:G],1,quantile,prob = 0.025)
theta.0.upper.2 = apply(theta.0.array[,burn:G],1,quantile,prob = 0.975)
theta.0.cov.2 = theta.0.lower.2 <= theta.0.true & theta.0.true <= theta.0.upper.2

theta.1.hat.2 = apply(theta.1.array[,burn:G],1,mean)
theta.1.sd.2 = apply(theta.1.array[,burn:G],1,sd)
theta.1.lower.2 = apply(theta.1.array[,burn:G],1,quantile,prob = 0.025)
theta.1.upper.2 = apply(theta.1.array[,burn:G],1,quantile,prob = 0.975)
theta.1.cov.2 = theta.1.lower.2 <= theta.1.true & theta.1.true <= theta.1.upper.2

ind.threshold.hat.2 = mean(ind.threshold.array[,burn:G])
ind.threshold.sd.2 = sd(ind.threshold.array[,burn:G])
ind.threshold.lower.2 = quantile(ind.threshold.array[,burn:G],prob = 0.025)
ind.threshold.upper.2 = quantile(ind.threshold.array[,burn:G],prob = 0.975)
ind.threshold.cov.2 = ind.threshold.lower.2 <= ind.threshold.true & ind.threshold.true <= ind.threshold.upper.2

pool.threshold.hat.2 = mean(pool.threshold.array[,burn:G])
pool.threshold.sd.2 = sd(pool.threshold.array[,burn:G])
pool.threshold.lower.2 = quantile(pool.threshold.array[,burn:G],prob = 0.025)
pool.threshold.upper.2 = quantile(pool.threshold.array[,burn:G],prob = 0.975)
pool.threshold.cov.2 = pool.threshold.lower.2 <= pool.threshold.true & pool.threshold.true <= pool.threshold.upper.2

###############################################
#####  Assess Convergence of Second MCMC  #####
###############################################
plot(beta.array[1,])
abline(h = beta.true[1])
plot(beta.array[2,])
abline(h = beta.true[2])
plot(beta.array[3,])
abline(h = beta.true[3])

plot(sigma.w.array[1,])
abline(h = sigma.w.true)

plot(theta.0.array[1,])
abline(h = theta.0.true[1])
plot(theta.0.array[2,])
abline(h = theta.0.true[2])
plot(theta.0.array[1,]/theta.0.array[2,])
abline(h = theta.0.true[1]/theta.0.true[2])
plot(theta.0.array[1,]/theta.0.array[2,]^2)
abline(h = theta.0.true[1]/theta.0.true[2]^2)
plot(theta.1.array[1,])
abline(h = theta.1.true[1])
plot(theta.1.array[2,])
abline(h = theta.1.true[2])
plot(theta.1.array[1,]/theta.1.array[2,])
abline(h = theta.1.true[1]/theta.1.true[2])
plot(theta.1.array[1,]/theta.1.array[2,]^2)
abline(h = theta.1.true[1]/theta.1.true[2]^2)

plot(zeta.tilde.array[1,])
abline(h=zeta.tilde.true[1])
plot(zeta.tilde.array[100,])
abline(h=zeta.tilde.true[100])
plot(zeta.tilde.array[200,])
abline(h=zeta.tilde.true[200])
plot(zeta.tilde.array[300,])
abline(h=zeta.tilde.true[300])

plot(W.tilde.array[1,])
abline(h = W.tilde.true.2[1])
plot(W.tilde.array[2,])
abline(h = W.tilde.true.2[2])
plot(W.tilde.array[3,])
abline(h = W.tilde.true.2[3])
plot(W.tilde.array[4,])
abline(h = W.tilde.true.2[4])

plot(ind.threshold.array[1,])
abline(h = ind.threshold.hat)
plot(pool.threshold.array[1,])
abline(h = pool.threshold.hat)

#####################################################################
#####  Indentify Individuals to Retest Based on the Second MCMC  ####
#####################################################################

###Classify Pools
pool.status = W > pool.threshold.hat.2

###Identify Individuals to Retest
mathcalP.retest.2 = list()
ct = 1

for(i in 1:N){
  pools.i = pools.containing.ind[[i]]
  if(sum(pool.status[pools.i])>=2 &  !is.element(i,retest.running.list) ){
    mathcalP.retest.2[[ct]] = i
    retest.running.list = c(retest.running.list,i)
    ct = ct +1
  }
}

individuals.restested.2 = unique(unlist(mathcalP.retest.2))
individuals.restested.2 = individuals.restested.2[!is.element(individuals.restested.2,individuals.restested)]
N.retest.2 = length(individuals.restested.2)
J.retest.2 = N.retest.2

###Retest the individuals
zeta.tilde.true.retest.2 = zeta.tilde.true[individuals.restested.2]
W.tilde.true.retest.2 = zeta.tilde.true.retest.2
W.retest.2 = rnorm(J.retest.2,W.tilde.true.retest.2,sqrt(sigma.w.true))
all.individuals.retested = c(individuals.restested,individuals.restested.2)
W.retest.all = c(W.retest,W.retest.2)

###Classify Individuals Based on Second Round of Testing
ind.status.2 = rep(0,N)
ct = 1
for(i in 1:N){
  if(is.element(i,all.individuals.retested)){
    ind = which(all.individuals.retested == i)
    if(W.retest.all[ind]>ind.threshold.hat.2){
      ind.status.2[i] = 1
    }
    ct = ct + 1
  }
}

beta.out = rbind(beta.true,beta.hat,beta.hat -beta.true,beta.sd,beta.cov)
beta.out = rbind(beta.out,beta.true,beta.hat.2,beta.hat.2 -beta.true,beta.sd.2,beta.cov.2)
write.table(beta.out,'beta.txt',row.names = FALSE,col.names = FALSE)

sigma.w.out = rbind(sigma.w.true,sigma.w.hat,sigma.w.hat -sigma.w.true,sigma.w.sd,sigma.w.cov)
sigma.w.out = rbind(sigma.w.out,sigma.w.true,sigma.w.hat.2,sigma.w.hat.2 -sigma.w.true,sigma.w.sd.2,sigma.w.cov.2)
write.table(sigma.w.out,'sigmaw.txt',row.names = FALSE,col.names = FALSE)

zeta.tilde.out = rbind(zeta.tilde.true,zeta.tilde.hat,zeta.tilde.hat -zeta.tilde.true,zeta.tilde.sd,zeta.tilde.cov)
zeta.tilde.out = rbind(zeta.tilde.out,zeta.tilde.true,zeta.tilde.hat.2,zeta.tilde.hat.2 -zeta.tilde.true,zeta.tilde.sd.2,zeta.tilde.cov.2)
write.table(zeta.tilde.out,'Ztilde.txt',row.names = FALSE,col.names = FALSE)

W.tilde.out = rbind(W.tilde.true,W.tilde.hat,W.tilde.hat -W.tilde.true,W.tilde.sd,W.tilde.cov)
write.table(W.tilde.out,'Wtilde1.txt',row.names = FALSE,col.names = FALSE)

W.tilde.true.2 = c(W.tilde.true,W.tilde.true.retest)
length.add = J.original + N - length(W.tilde.true.2)
W.tilde.out = rbind(c(W.tilde.true.2,rep(-99,length.add)),c(W.tilde.hat.2,rep(-99,length.add)),c(W.tilde.hat.2 -W.tilde.true.2,rep(-99,length.add)),c(W.tilde.sd.2,rep(-99,length.add)),c(W.tilde.cov.2,rep(-99,length.add)))
write.table(W.tilde.out,'Wtilde2.txt',row.names = FALSE,col.names = FALSE)

length.add = J.original + N - length(W.retest.all)
W.tilde.out = c(W.retest.all,rep(-99,length.add))
write.table(W.tilde.out,'Wtilde3.txt',row.names = FALSE,col.names = FALSE)

Y.tilde.out = rbind(Y.tilde.true,Y.tilde.hat,Y.tilde.hat -Y.tilde.true,Y.tilde.sd,Y.tilde.cov)
Y.tilde.out = rbind(Y.tilde.out,Y.tilde.true,Y.tilde.hat.2,Y.tilde.hat.2 -Y.tilde.true,Y.tilde.sd.2,Y.tilde.cov.2)
write.table(Y.tilde.out,'Ytilde.txt',row.names = FALSE,col.names = FALSE)

theta.out = rbind(c(theta.0.true,theta.1.true),c(theta.0.hat,theta.1.hat),c(theta.0.true,theta.1.true)-c(theta.0.hat,theta.1.hat),c(theta.0.sd,theta.1.sd),c(theta.0.cov,theta.1.cov))
theta.out = rbind(theta.out,c(theta.0.true,theta.1.true),c(theta.0.hat.2,theta.1.hat.2),c(theta.0.true,theta.1.true)-c(theta.0.hat.2,theta.1.hat.2),c(theta.0.sd.2,theta.1.sd.2),c(theta.0.cov.2,theta.1.cov.2))
write.table(theta.out,'theta.txt',row.names = FALSE,col.names = FALSE)

threshold.out = rbind(c(ind.threshold.true,pool.threshold.true),c(ind.threshold.hat,pool.threshold.hat),c(ind.threshold.true,pool.threshold.true)-c(ind.threshold.hat,pool.threshold.hat),c(ind.threshold.sd,pool.threshold.sd),c(ind.threshold.cov,pool.threshold.cov))
threshold.out = rbind(threshold.out,c(ind.threshold.true,pool.threshold.true),c(ind.threshold.hat.2,pool.threshold.hat.2),c(ind.threshold.true,pool.threshold.true)-c(ind.threshold.hat.2,pool.threshold.hat.2),c(ind.threshold.sd.2,pool.threshold.sd.2),c(ind.threshold.cov.2,pool.threshold.cov.2))
write.table(threshold.out,'threshold.txt',row.names = FALSE,col.names = FALSE)

classification.out = rbind(Y.tilde.true,ind.status.orginal,ind.status.2)
write.table(classification.out,'classification.txt',row.names = FALSE,col.names = FALSE)

