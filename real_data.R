rm(list=ls())

library(circular)
library(CircStats)
library(splines2)
library(numDeriv)

library(parallel)
library(doParallel)
library(foreach)

# l2 norm
l2 = function(v){
  sqrt(sum(v^2))
}

#projected data
x_pr = function(beta, x){
  as.vector(x %*% (beta / l2(beta)))
}

beta = function(phi) c(1/sqrt(1+t(phi)%*%phi), as.numeric(1/sqrt(1+t(phi)%*%phi))*phi)

jac.beta = function(phi){
  a1 = -as.numeric(1/(1+t(phi)%*%phi)^(3/2))*phi
  a2 = as.numeric(1/sqrt(1+t(phi)%*%phi))*diag(1,2) - as.numeric(1/(1+t(phi)%*%phi)^(3/2))*(phi%*%t(phi))
  return(matrix(rbind(a1,a2), 3, 2, byrow = T))
}

# Shortest arc-length distance
cdist = function(a, b){
  delta = abs(a - b)     
  pmin(delta, 2*pi - delta)
}


# Original covariates x0
# Centred covariates x
x = apply(x0, 2, function(x) x-mean(x)) 


######### Fit Fisher-Lee model ###########

fl = suppressWarnings(lm.circular(y, x, init=rnorm(ncol(x),0, 0.3), type='c-l'))



##########  Fit ZIvM 1 model ###########


n.knots = 5
d.f = n.knots + 2
int_b = function(x, df = d.f, degree = 2) {
  ibs_x = ibs(x, df, degree = degree)
  ibs_0 = predict(ibs_x, 0)
  sweep(ibs_x, 2, ibs_0, "-")
}

E = function(x, y, phi, gamma, a, mu){
  b = int_b(x_pr(beta(phi), x))
  n = length(y)
  S = 0
  for(i in 1:n){
    S = S + ( (1-a[i])*cos(y[i]-mu-2*atan(sum(gamma*b[i,]))) )
  }
  return(-S)
}

loglik = function(theta, x, y){
  n = length(y)
  b = int_b(x_pr(beta(theta[4:5]), x))
  mu = rep(0, n)
  for(i in 1:n){
    mu[i] = theta[2]+2*atan(sum(theta[-c(1:5)]*b[i,]))
  }
  sum(log(theta[1]*dvm(y, 0, 1000) + (1-theta[1])*dvm(y, mu, theta[3]) ))
}

lc.em = function(x, y, init, tol = 1e-3, max.iter = 200){
  n = length(y)
  theta.0 = init
  diff = 1
  iter = 0
  while(diff > tol && iter <= max.iter){
    
    theta.1 = theta.0
    
    b = int_b(x_pr(beta(theta.0[4:5]), x))
    
    mu = rep(0,n)
    for(i in 1:n){
      mu[i] = theta.0[2]+2*atan(sum(theta.0[-c(1:5)]*b[i,]))
    }
    
    a = (theta.0[1]*dvm(y,0,1000))/(theta.0[1]*dvm(y,0,1000)+(1-theta.0[1])*dvm(y,mu,theta.0[3]))
    
    gamma_new = optim(theta.0[-c(1:5)], function(gamma)E(x, y, theta.0[4:5], gamma, a, theta.0[2]), method = "BFGS", control=list(maxit=2000))$par
    phi_new = optim(theta.0[4:5], function(phi)E(x, y, phi, gamma_new, a, theta.0[2]), method = "BFGS", control=list(maxit=2000))$par
    
    #E_opt = E2_opt(x, y, theta.0[4:5], theta.0[-c(1:5)], a, theta.0[2])
    #gamma_new = E_opt[-c(1:2)]
    #phi_new = E_opt[1:2]
    
    p_new = mean(a)
    
    b = int_b(x_pr(beta(phi_new), x))
    
    cs = rep(0, n)
    sn = rep(0, n)
    for(i in 1:n){
      cs[i] = cos(y[i]-2*atan(sum(gamma_new*b[i,])))
      sn[i] = sin(y[i]-2*atan(sum(gamma_new*b[i,])))
    }
    
    C = weighted.mean(cs, 1-a)
    D = weighted.mean(sn, 1-a)
    
    mu_new = atan2(D, C) 
    kappa_new = A1inv(sqrt(C^2 + D^2))
    theta.0 = c(p_new, mu_new, kappa_new, phi_new, gamma_new)
    #diff = l2(theta.0 - theta.1)
    diff = l2(theta.0[-c(1:3)] - theta.1[-c(1:3)])
    iter = iter + 1
  }
  if (diff <= tol) {
    return(theta.0)   # success: return estimate
  } else {
    return(NULL)      # fail: max.iter reached without convergence
  }
  
}


zivm1 = function(x, y, n.knots = 3, m = 200){
  
  #Fit Fisher-Lee model
  fl = suppressWarnings(lm.circular(y, x, init = rnorm(ncol(x), 0, 0.3), type='c-l'))
  
  #d.f = n.knots + 2
  #n = length(y)
  
  # Initial values
  p0 = sum(y==0)/length(y)
  #p.init = runif(m, max(0, p0 - 0.1), min(1, p0 + 0.1))
  #mu.init = runif(m, as.numeric(fl$mu)-0.3, as.numeric(fl$mu)+0.3)
  #kappa.init = rexp(m, 1/fl$kappa)
  beta0 = fl$coefficients
  phi0 = beta0[-1]/beta0[1]
  #phi.init = MASS::mvrnorm(m, mu = phi0, Sigma = 0.1^2 * diag(ncol(x)-1))
  phi.init = matrix(phi0, nrow = m, ncol = length(phi0), byrow = TRUE)
  gamma.init = matrix(0, nrow = m, ncol = d.f)
  for(i in 1:m){
    gamma.init[i,] = rnorm(d.f, 0, 0.3)
  }
  #init.em = unname(cbind(p.init, mu.init, kappa.init, phi.init, gamma.init))
  init.em = cbind(p0, as.numeric(fl$mu), fl$kappa, phi.init, gamma.init)
  
  # Parallel computation 
  ncores = detectCores()
  cl = makeCluster(ncores-1)
  clusterExport(cl,ls(globalenv()),envir = globalenv())
  registerDoParallel(cl)
  t1 = Sys.time()
  res = foreach(i = 1:m, .combine=rbind, .packages = c("splines2","CircStats", "circular"))%dopar%{
    #cat(i,"\t")
    theta.em = lc.em(x, y, init.em[i,])
    ll= loglik(theta.em, x, y)
    c(theta.em, ll)
  }
  stopCluster(cl)
  t2 = Sys.time()
  t2 - t1
  
  ll = res[, ncol(res)]
  theta.em = res[, -ncol(res)]
  theta_hat = theta.em[which.max(ll),]
  
  aic = 2*length(theta_hat) - 2*loglik(theta_hat, x, y)
  bic = length(theta_hat)*log(n) - 2*loglik(theta_hat, x, y)
  est = c(theta_hat[1:3], beta(theta_hat[4:5]), theta_hat[-c(1:5)])
  cov = solve(hessian(function(theta)-loglik(theta, x, y),theta_hat))
  se = sqrt(diag(cov))
  cov.beta = jac.beta(theta_hat[4:5])%*%cov[4:5,4:5]%*%t(jac.beta(theta_hat[4:5]))
  se.beta = sqrt(diag(cov.beta))
  se_est=c(se[1:3], se.beta, se[-c(1:5)])
  b = est[4:6]
  g = est[-c(1:6)]
  
  # Estimated mean direction
  B_0 = int_b(x_pr(b, x))
  f_0 = rep(0, nrow(x))
  for(i in 1:nrow(x)){
    f_0[i] = sum(g*B_0[i,])
  }
  mu_ZIVM1 = (est[2] + 2*atan(f_0)) %% (2*pi)
  
  # Sum of residuals (to be computed for non-zero responses and corresponding estimated mean directions)
  id = (y != 0)
  e1 = sum(cdist(y[id], mu_ZIVM1[id]))
  
  #CRPS
  n.sim = 1e6
  set.seed(1)
  sim = runif(n.sim, 0, 2*pi)
  sim1 = runif(n.sim, 0, 2*pi)
  sim2 = runif(n.sim, 0, 2*pi)
  crp11 = sapply(1:n, function(i){
    2*pi * mean(cdist(sim, y[i]) * (est[1]*(CircStats::dvm(sim,0,1000)) + (1-est[1])*(CircStats::dvm(sim,mu_ZIVM1[i],est[3]))))
  })
  crp12 = sapply(1:n, function(i){
    4*pi^2 * mean(cdist(sim1, sim2) *
                    (est[1]*(CircStats::dvm(sim1,0,1000)) + (1-est[1])*(CircStats::dvm(sim1,mu_ZIVM1[i],est[3]))) *
                    (est[1]*(CircStats::dvm(sim2,0,1000)) + (1-est[1])*(CircStats::dvm(sim2,mu_ZIVM1[i],est[3]))))
  })
  crps1 = crp11-0.5*crp12
  
  result = list(est[1:6], se_est[1:6], aic, bic, e1, mean(crps1))
  names(result)=c("Estimate","Std.Error","AIC","BIC","Sum_Residual","Avg.CRPS")
  names(result[[1]]) = c("p", "mu", "kappa","beta_1","beta_2","beta_3")
  names(result[[2]]) = c("p", "mu", "kappa","beta_1","beta_2","beta_3")
  return(lapply(result,round,4))
}


t1 = Sys.time()
fit1=zivm1(x,y)
t2 = Sys.time()
t2-t1
fit1



#N=1
fit1
$Estimate
p      mu   kappa  beta_1  beta_2  beta_3 
0.2914  0.4532  2.9292  0.3675 -0.8287  0.4222 

$Std.Error
p     mu  kappa beta_1 beta_2 beta_3 
0.0677 0.1865 0.5846 0.0429 0.0440 0.0809 

$AIC
[1] 68.1498

$BIC
[1] 84.0617

$Sum_Residual
[1] 17.9988

$Avg.CRPS
[1] 0.309



#N=2
$Estimate
p      mu   kappa  beta_1  beta_2  beta_3 
0.2994  0.4745  2.9916  0.5119 -0.8176  0.2635 

$Std.Error
p     mu  kappa beta_1 beta_2 beta_3 
0.0673 0.2039 0.5973 0.0375 0.0185 0.0492 

$AIC
[1] 68.7756

$BIC
[1] 86.6764

$Sum_Residual
[1] 18.6639

$Avg.CRPS
[1] 0.3226



#N=3
$Estimate
p      mu   kappa  beta_1  beta_2  beta_3 
0.2867  0.6658  3.2714  0.2059 -0.8483  0.4879 

$Std.Error
p     mu  kappa beta_1 beta_2 beta_3 
0.0681 0.2188 0.6609 0.0258 0.0575 0.0661 

$AIC
[1] 66.6917

$BIC
[1] 86.5815

$Sum_Residual
[1] 17.449

$Avg.CRPS
[1] 0.2889



#N=4
> fit1
$Estimate
p      mu   kappa  beta_1  beta_2  beta_3 
0.2881  0.8264  3.4449  0.2826 -0.8594  0.4261 

$Std.Error
p     mu  kappa beta_1 beta_2 beta_3 
0.0679 0.2414 0.7005 0.0220 0.0565 0.0548 

$AIC
[1] 66.0506

$BIC
[1] 87.9295

$Sum_Residual
[1] 17.2091

$Avg.CRPS
[1] 0.2824


#N=5
> fit1
$Estimate
p      mu   kappa  beta_1  beta_2  beta_3 
0.2840  0.7511  3.5351  0.2410 -0.8665  0.4372 

$Std.Error
p     mu  kappa beta_1 beta_2 beta_3 
0.0683 0.2332 0.7217 0.0228 0.0425 0.0522 

$AIC
[1] 66.8805

$BIC
[1] 90.7483

$Sum_Residual
[1] 16.9735

$Avg.CRPS
[1] 0.274



#######################################################################


n.knots = 4
d.f = n.knots + 2
int_b = function(x, df = d.f, degree = 2) {
  ibs_x = ibs(x, df, degree = degree)
  ibs_0 = predict(ibs_x, 0)
  sweep(ibs_x, 2, ibs_0, "-")
}

# Generate initial values
m = 200
#set.seed(1)
p0 = sum(y==0)/length(y)
#p.init = runif(m, max(0, p0 - 0.1), min(1, p0 + 0.1))
#mu.init = runif(m, as.numeric(fl$mu)-0.3, as.numeric(fl$mu)+0.3)
#kappa.init = rexp(m, 1/fl$kappa)
beta0 = fl$coefficients
phi0 = beta0[-1]/beta0[1]
#phi.init = MASS::mvrnorm(m, mu = phi0, Sigma = 0.1^2 * diag(ncol(x)-1))
phi.init = matrix(phi0, nrow = m, ncol = length(phi0), byrow = TRUE)
gamma.init = matrix(0, nrow = m, ncol = d.f)
for(i in 1:m){
  gamma.init[i,] = rnorm(d.f, 0, 0.3)
}
#init.em = unname(cbind(p.init, mu.init, kappa.init, phi.init, gamma.init))
init.em = cbind(p0, as.numeric(fl$mu), fl$kappa, phi.init, gamma.init)


# Parallel computation 
ncores = detectCores()
cl = makeCluster(ncores-1)
registerDoParallel(cl)
t1 = Sys.time()
res = foreach(i = 1:m, .combine=rbind, .packages = c("splines2","CircStats", "circular"),.errorhandling = "remove")%dopar%{
  #cat(i,"\t")
  theta.em = lc.em(x, y, init.em[i,])
  ll= loglik(theta.em, x, y)
  c(theta.em, ll)
}
stopCluster(cl)
t2 = Sys.time()
t2 - t1

#nrow(res)

ll = res[, ncol(res)]
theta.em = res[, -ncol(res)]
theta_hat = theta.em[which.max(ll),]
aic = 2*length(theta_hat) - 2*loglik(theta_hat, x, y)
bic = length(theta_hat)*log(length(y)) - 2*loglik(theta_hat, x, y)
est = c(theta_hat[1:3], beta(theta_hat[4:5]), theta_hat[-c(1:5)])
cov = solve(hessian(function(theta)-loglik(theta, x, y),theta_hat))
se = sqrt(diag(cov))
cov.beta = jac.beta(theta_hat[4:5])%*%cov[4:5,4:5]%*%t(jac.beta(theta_hat[4:5]))
se.beta = sqrt(diag(cov.beta))
se_est=c(se[1:3], se.beta, se[-c(1:5)])
b = est[4:6]
g = est[-c(1:6)]

# Estimated mean direction
B_0 = int_b(x_pr(b, x))
f_0 = rep(0, nrow(x))
for(i in 1:nrow(x)){
  f_0[i] = sum(g*B_0[i,])
}
mu_ZIVM1 = (est[2] + 2*atan(f_0)) %% (2*pi)

# Sum of residuals (to be computed for non-zero responses and corresponding estimated mean directions)
e1 = sum(cdist(y[1:37], mu_ZIVM1[1:37]))

#CRPS
n.sim = 1e6
set.seed(1)
sim = runif(n.sim, 0, 2*pi)
sim1 = runif(n.sim, 0, 2*pi)
sim2 = runif(n.sim, 0, 2*pi)
crp11 = sapply(1:n, function(i){
  2*pi * mean(cdist(sim, y[i]) * (est[1]*(CircStats::dvm(sim,0,1000)) + (1-est[1])*(CircStats::dvm(sim,mu_ZIVM1[i],est[3]))))
})
crp12 = sapply(1:n, function(i){
  4*pi^2 * mean(cdist(sim1, sim2) *
                  (est[1]*(CircStats::dvm(sim1,0,1000)) + (1-est[1])*(CircStats::dvm(sim1,mu_ZIVM1[i],est[3]))) *
                  (est[1]*(CircStats::dvm(sim2,0,1000)) + (1-est[1])*(CircStats::dvm(sim2,mu_ZIVM1[i],est[3]))))
})
crps1 = crp11-0.5*crp12

result = list(est[1:6],se_est[1:6],aic,bic,e1,mean(crps1))
names(result)=c("Estimate","Std.Error","AIC","BIC","Sum_Residual","Avg.CRPS")
names(result[[1]]) = c("p", "mu", "kappa","beta_1","beta_2","beta_3")
names(result[[2]]) = c("p", "mu", "kappa","beta_1","beta_2","beta_3")
lapply(result,round,4)


######################################################




# Plotting estimated function
N = 100
X1 = X2 = X3 = seq(-1, 1, length = N)
X = cbind(X1, X2, X3)
#b = est[4:6]
B = int_b(x_pr(b, X))
#g = est[-c(1:6)]

V = cov[-c(1:5),-c(1:5)]
var.h = conf.lower = conf.upper = rep(0, N)
for(i in 1:N){
  var.h[i] = t(B[i,])%*%V%*%B[i,]
  conf.lower[i] = t(g)%*%B[i,] - 1.96*sqrt(var.h[i])
  conf.upper[i] = t(g)%*%B[i,] + 1.96*sqrt(var.h[i])
}
f_1 = rep(0, N)
for(i in 1:N){
  f_1[i] = sum(g*B[i,])
}
plot(x_pr(b, X), f_1, type = "l", lwd=2, xlab=expression(beta^T * X), ylab=expression(h(beta^T * X)), ylim = c(-0.28,0.32))
lines(x_pr(b, X), conf.lower, lty=2)
lines(x_pr(b, X), conf.upper, lty=2)
h = function(x) x
curve(h, add = TRUE, lty=6, lwd=2)
legend("topleft", legend=c("ZIvM 1", "Confidence band for ZIvM 1","Fisher-Lee"), lty=c(1,2,6), lwd=c(2,1,2),cex=1.1,y.intersp=0.5)   




u1 = seq(min(x_pr(b, x)), max(x_pr(b, x)), length = N)
B = int_b(u1)
fhat = as.vector(B %*% g)
var.h = rowSums((B %*% V) * B)
lower = fhat - 1.96 * sqrt(var.h)
upper = fhat + 1.96 * sqrt(var.h)

#N=3
par(mar = c(5, 4.5, 4, 2) + 0.1)
plot(u1, fhat, type="l", lwd=2,xlab=expression(hat(beta)^T * x), ylab=expression(hat(h)(hat(beta)^T * x)), ylim = c(-1.75,2.25))
lines(u1, lower, lty=2)
lines(u1, upper, lty=2)

#N=4
par(mar = c(5, 4.5, 4, 2) + 0.1)
plot(u1, fhat, type="l", lwd=2,xlab=expression(hat(beta)^T * x), ylab=expression(hat(h)(hat(beta)^T * x)), ylim = c(-7.7,1.6))
lines(u1, lower, lty=2)
lines(u1, upper, lty=2)
legend("bottomleft", legend=c("ZIvM 1", "Confidence band for ZIvM 1","Fisher-Lee"), lty=c(1,2,6), lwd=c(2,1,2),cex=1.1,y.intersp=0.5)   


#mu.hat = (est[2] + 2 * atan(fhat)) %% (2*pi)
#plot(u1, mu.hat,type="l")



# Plotting Donut-plots of ZIvM 1 and Fisher-Lee

mu_fl = (as.numeric(fl$mu) + 2*atan(as.vector(x %*% fl$coefficients))) %% (2*pi) # Estimated mean direction in Fisher-Lee model

plot_circ_diag = function(y, y_hat) {
  stopifnot(length(y) == length(y_hat))
  
  delta = (y_hat - y) %% (2*pi) 
  
  # radius
  r = 1 + cos(delta)
  
  x = r * cos(y_hat)
  y_cart = r * sin(y_hat)
  
  # setup plot
  plot(0,0, type="n", asp=1, xlim=c(-2.5,2.5), ylim=c(-2.5,2.5),
       axes=FALSE, xlab="", ylab="")
  symbols(rep(0,2), rep(0,2), circles=c(1,2), inches=FALSE, add=TRUE, lwd=1)
  
  # center dot
  points(0,0, pch=19, col="black")
  
  # split into clockwise vs anticlockwise
  cw  = delta > pi   # clockwise → solid
  acw = delta > 0 & delta <= pi   # anticlockwise → open
  
  # plot points
  points(x[cw],  y_cart[cw],  pch=19, col="blue",cex=2)   # filled
  points(x[acw], y_cart[acw], pch=1,  col="blue",cex=2)   # open
}

par(mfrow = c(1,2),mar = c(3,2,2,1),oma = c(0,0,0,0))    # shrink subplot margins
plot_circ_diag(y[1:37], mu_fl[1:37])
plot_circ_diag(y[1:37], mu_ZIVM1[1:37])






############# Fit ZIvM 2 model ###############


z = cbind(1, x)

p = function(delta, x){
  eta = as.vector(x %*% delta)
  1/(1+exp(-eta))
}

w = function(delta, x){
  diag(p(delta, x)*(1-p(delta, x)))
}

lc.em2 = function(x,z, y, init, tol = 1e-3, max.iter = 200){
  n = length(y)
  theta.0 = init
  
  b = int_b(x_pr(beta(theta.0[7:8]), x))
  mu = rep(0,n)
  for(i in 1:n){
    mu[i] = theta.0[5]+2*atan(sum(theta.0[-c(1:8)]*b[i,]))
  }
  a = (p(theta.0[1:4],z)*dvm(y,0,1000))/(p(theta.0[1:4],z)*dvm(y,0,1000)+(1-p(theta.0[1:4],z))*dvm(y,mu,theta.0[6]))
  
  diff = 1
  iter = 0
  while(diff > tol && iter <= max.iter){
    theta.1 = theta.0
    
    b = int_b(x_pr(beta(theta.0[7:8]), x))
    
    mu = rep(0,n)
    for(i in 1:n){
      mu[i] = theta.0[5]+2*atan(sum(theta.0[-c(1:8)]*b[i,]))
    }
    
    a = (p(theta.0[1:4],z)*dvm(y,0,1000))/(p(theta.0[1:4],z)*dvm(y,0,1000)+(1-p(theta.0[1:4],z))*dvm(y,mu,theta.0[6]))
    
    gamma_new = optim(theta.0[-c(1:8)], function(gamma)E(x, y, theta.0[7:8], gamma, a, theta.0[5]), method = "BFGS", control=list(maxit=2000))$par
    phi_new = optim(theta.0[7:8], function(phi)E(x, y, phi, gamma_new, a, theta.0[5]), method = "BFGS", control=list(maxit=2000))$par
    
    #E_opt = E2_opt(x, y, theta.0[7:8], theta.0[-c(1:8)], a, theta.0[5])
    #gamma_new = E_opt[-c(1:2)]
    #phi_new = E_opt[1:2]
    
    delta_new = theta.0[1:4] + solve(t(z) %*% w(theta.0[1:4],z) %*% z) %*% t(z) %*% (a-p(theta.0[1:4],z))
    #delta_new = theta.0[1:4] + ginv(t(z) %*% w(theta.0[1:4],z) %*% z) %*% t(z) %*% (a-p(theta.0[1:4],z))
    delta_new = as.numeric(delta_new)
    
    b = int_b(x_pr(beta(phi_new), x))
    
    cs = rep(0, n)
    sn = rep(0, n)
    for(i in 1:n){
      cs[i] = cos(y[i]-2*atan(sum(gamma_new*b[i,])))
      sn[i] = sin(y[i]-2*atan(sum(gamma_new*b[i,])))
    }
    
    C = weighted.mean(cs, 1-a)
    D = weighted.mean(sn, 1-a)
    
    mu_new = atan2(D, C) 
    kappa_new = A1inv(sqrt(C^2 + D^2))
    theta.0 = c(delta_new, mu_new, kappa_new, phi_new, gamma_new)
    diff = l2(theta.0 - theta.1)
    iter = iter + 1
  }
  #return(theta.0)
  if (diff <= tol) {
    return(theta.0)   # success: return estimate
  } else {
    return(NULL)      # fail: max.iter reached without convergence
  }
}

loglik2 = function(theta, x,z, y){
  n = length(y)
  b = int_b(x_pr(beta(theta[7:8]), x))
  mu = rep(0, n)
  for(i in 1:n){
    mu[i] = theta[5]+2*atan(sum(theta[-c(1:8)]*b[i,]))
  }
  sum(log(p(theta[1:4],z)*dvm(y, 0, 1000) + (1-p(theta[1:4],z))*dvm(y, mu, theta[6]) ))
}


# Generate initial values
m = 200
#set.seed(1)
delta.init = matrix(0, nrow = m, ncol = ncol(z))
mu.init = runif(m, as.numeric(fl$mu)-0.3, as.numeric(fl$mu)+0.3)
kappa.init = rexp(m, 1/fl$kappa)
beta0 = fl$coefficients
phi0 = beta0[-1]/beta0[1]
phi.init = MASS::mvrnorm(m, mu = phi0, Sigma = 0.1^2 * diag(ncol(x)-1))
gamma.init = matrix(0, nrow = m, ncol = d.f)
for(i in 1:m){
  delta.init[i,] = rnorm(ncol(z), 0, 0.3)
  gamma.init[i,] = rnorm(d.f, 0, 0.3)
}

init.em = unname(cbind(delta.init, mu.init, kappa.init, phi.init, gamma.init))

ncores = detectCores()
cl = makeCluster(ncores-1)
registerDoParallel(cl)
t1 = Sys.time()
res1 = foreach(i = 1:m, .combine=rbind, .packages = c("splines2","CircStats", "circular"),.errorhandling = "remove")%dopar%{
  #cat(i,"\t")
  theta.em = lc.em2(x,z, y, init.em[i,])
  ll= loglik2(theta.em, x,z, y)
  c(theta.em, ll)
}
stopCluster(cl)
t2 = Sys.time()
t2 - t1

nrow(res1)/m

ll = res1[, ncol(res1)]
theta.em = res1[, -ncol(res1)]
theta_hat1 = theta.em[which.max(ll),]

fit1 = 2*length(theta_hat1) - 2*loglik2(theta_hat1, x,z, y)

bic1 = length(theta_hat1)*log(length(y)) - 2*loglik2(theta_hat1, x,z, y)

est1 = c(theta_hat1[1:6], beta(theta_hat1[7:8]), theta_hat1[-c(1:8)])

cov1 = solve(hessian(function(theta)-loglik2(theta, x,z, y),theta_hat1), method.args = list(eps = 1e-6))
se1 = sqrt(diag(cov1))

cov1.beta = jac.beta(theta_hat1[7:8])%*%cov1[7:8,7:8]%*%t(jac.beta(theta_hat1[7:8]))
se1.beta = sqrt(diag(cov1.beta))

N = 100
X1 = X2 = seq(-1, 1, length = N)
X = cbind(X1, X2)

b1 = est1[6:7]
B1 = int_b(x_pr(b1, X))
g1 = est1[-c(1:7)]
f_1 = rep(0, N)
for(i in 1:N){
  f_1[i] = sum(g1*B1[i,])
}
plot(x_pr(b1, X), f_1, type = "l")


r1 = data.frame(cbind(est1, c(se1[1:5], se1.beta, se1[-c(1:6)])))
colnames(r1) = c("Estimate", "Std. Error")
rownames(r1)[4:7] = c("mu", "kappa","beta_1","beta_2")

result1 = list(r1, fit1, bic1)
names(result1[[2]])="AIC"
names(result1[[3]])="BIC"

round(result1[[1]], 4)
round(result1[[2]], 4)
round(result1[[3]], 4)


B_1 = int_b(x_pr(b1, x))
f_1 = rep(0, nrow(x))
for(i in 1:nrow(x)){
  f_1[i] = sum(g1*B_1[i,])
}
mu1_ZIVM = (est1[5] + 2*atan(f_1)) %% (2*pi)





