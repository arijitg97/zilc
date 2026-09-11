library(circular)
library(CircStats)
library(splines2)
library(numDeriv)

library(parallel)
library(doParallel)
library(foreach)

x_centre = function(x) x - mean(x)
x = apply(x0, 2, x_centre)

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

E.vmsp = function(x, y, phi, gamma, mu){
  b = int_b(x_pr(beta(phi), x))
  n = length(y)
  S = 0
  for(i in 1:n){
    S = S + cos(y[i]-mu-2*atan(sum(gamma*b[i,]))) 
  }
  return(-S)
}



lc.vmsp = function(x, y, init, tol = 1e-3, max.iter = 200){
  n = length(y)
  theta.0 = init
  diff = 1
  iter = 0
  
  while(diff > tol && iter <= max.iter){
    theta.1 = theta.0
    
    b = int_b(x_pr(beta(theta.0[3:4]), x))
    
    gamma_new = optim(theta.0[-c(1:4)], function(gamma)E.vmsp(x, y, theta.0[3:4], gamma, theta.0[1]), method = "BFGS", control=list(maxit=2000))$par
    phi_new = optim(theta.0[3:4], function(phi)E.vmsp(x, y, phi, gamma_new, theta.0[1]), method = "BFGS", control=list(maxit=2000))$par
    
    b = int_b(x_pr(beta(phi_new), x))
    
    cs = rep(0, n)
    sn = rep(0, n)
    for(i in 1:n){
      cs[i] = cos(y[i]-2*atan(sum(gamma_new*b[i,])))
      sn[i] = sin(y[i]-2*atan(sum(gamma_new*b[i,])))
    }
    
    C = mean(cs)
    D = mean(sn)
    
    mu_new = atan2(D, C) 
    kappa_new = A1inv(sqrt(C^2 + D^2))
    theta.0 = c(mu_new, kappa_new, phi_new, gamma_new)
    #diff = l2(theta.0 - theta.1)
    diff = l2(theta.0[-c(1:2)] - theta.1[-c(1:2)])
    iter = iter + 1
  }
  if (diff <= tol) {
    return(theta.0)   
  } else {
    return(NULL)     
  }
}

loglik.vmsp = function(theta, x, y){
  n = length(y)
  b = int_b(x_pr(beta(theta[3:4]), x))
  mu = rep(0, n)
  for(i in 1:n){
    mu[i] = theta[1]+2*atan(sum(theta[-c(1:4)]*b[i,]))
  }
  sum(log(dvm(y, mu, theta[2]) ))
}



########################################################

vmsp = function(x, y, m = 200){
  #Fit Fisher-Lee model
  fl = suppressWarnings(lm.circular(y, x, init = rnorm(ncol(x), 0, 0.3), type='c-l'))
  beta0 = fl$coefficients
  phi0 = beta0[-1]/beta0[1]
  phi.init = matrix(phi0, nrow = m, ncol = length(phi0), byrow = TRUE)
  gamma.init = matrix(0, nrow = m, ncol = d.f)
  for(i in 1:m){
    gamma.init[i,] = rnorm(d.f, 0, 0.3)
  }
  init = cbind(as.numeric(fl$mu), fl$kappa, phi.init, gamma.init)
  
  # Parallel computation 
  ncores = detectCores()
  cl = makeCluster(ncores-1)
  clusterExport(cl,ls(globalenv()),envir = globalenv())
  registerDoParallel(cl)
  res = foreach(i = 1:m, .combine=rbind, .packages = c("splines2","CircStats", "circular"))%dopar%{
    theta.m = lc.vmsp(x, y, init[i,])
    ll= loglik.vmsp(theta.m, x, y)
    c(theta.m, ll)
  }
  stopCluster(cl)
  
  ll = res[, ncol(res)]
  theta.m = res[, -ncol(res)]
  theta.hat = theta.m[which.max(ll),]
  
  aic = 2*length(theta.hat) - 2*loglik.vmsp(theta.hat, x, y)
  bic = length(theta.hat)*log(length(y)) - 2*loglik.vmsp(theta.hat, x, y)
  est = c(theta.hat[1:2], beta(theta.hat[3:4]), theta.hat[-c(1:4)])
  cov = solve(hessian(function(theta)-loglik.vmsp(theta, x, y),theta.hat))
  se = sqrt(diag(cov))
  cov.beta = jac.beta(theta.hat[3:4])%*%cov[3:4,3:4]%*%t(jac.beta(theta.hat[3:4]))
  se.beta = sqrt(diag(cov.beta))
  se_est=c(se[1:2], se.beta, se[-c(1:4)])
  b = est[3:5]
  g = est[-c(1:5)]
  
  # Estimated mean direction
  B_0 = int_b(x_pr(b, x))
  f_0 = rep(0, nrow(x))
  for(i in 1:nrow(x)){
    f_0[i] = sum(g*B_0[i,])
  }
  mu.vmsp = (est[1] + 2*atan(f_0)) %% (2*pi)
  
  # Sum of residuals (to be computed for non-zero responses and corresponding estimated mean directions)
  id = (y != 0)
  e1 = sum(cdist(y[id], mu.vmsp[id]))
  
  #CRPS
  n.sim = 1e6
  set.seed(1)
  sim = runif(n.sim, 0, 2*pi)
  sim1 = runif(n.sim, 0, 2*pi)
  sim2 = runif(n.sim, 0, 2*pi)
  crp11 = sapply(1:n, function(i){
    2*pi * mean(cdist(sim, y[i]) * (CircStats::dvm(sim, mu.vmsp[i], est[2])))
  })
  crp12 = sapply(1:n, function(i){
    4*pi^2 * mean(cdist(sim1, sim2) *
                    (CircStats::dvm(sim1, mu.vmsp[i], est[2])) *
                    (CircStats::dvm(sim2, mu.vmsp[i], est[2])))
  })
  crps1 = crp11-0.5*crp12
  
  result = list(est[1:5], se_est[1:5], aic, bic, e1, mean(crps1))
  names(result)=c("Estimate","Std.Error","AIC","BIC","Sum_Residual","Avg.CRPS")
  names(result[[1]]) = c("mu", "kappa","beta_1","beta_2","beta_3")
  names(result[[2]]) = c("mu", "kappa","beta_1","beta_2","beta_3")
  return(lapply(result,round,4))
}


n.knots = 3
d.f = n.knots + 2
int_b = function(x, df = d.f, degree = 2) {
  ibs_x = ibs(x, df, degree = degree)
  ibs_0 = predict(ibs_x, 0)
  sweep(ibs_x, 2, ibs_0, "-")
}

t1=Sys.time()
fit2 = vmsp(x,y)
t2=Sys.time()
t2-t1
fit2

plot_circ_diag(y[1:37], mu.vmsp[1:37])



V = cov[-c(1:4),-c(1:4)]
N = 100
u1 = seq(min(x_pr(b, x)), max(x_pr(b, x)), length = N)
B = int_b(u1)
fhat = as.vector(B %*% g)
var.h = rowSums((B %*% V) * B)
lower = fhat - 1.96 * sqrt(var.h)
upper = fhat + 1.96 * sqrt(var.h)

par(mar = c(5, 4.5, 4, 2) + 0.1)
plot(u1, fhat, type="l", lwd=2,xlab=expression(hat(beta)^T * x), ylab=expression(hat(h)(hat(beta)^T * x)), ylim = c(-1.55,1.6))
lines(u1, lower, lty=2)
lines(u1, upper, lty=2)
h = function(x) x
curve(h, add = TRUE, lty=6, lwd=2)
legend("topleft", legend=c("vMSp", "Confidence band for vMSp","Fisher-Lee"), lty=c(1,2,6), lwd=c(2,1,2),cex=1.1,y.intersp=0.5)   



####### LOOCV #########


n.knots = 2
d.f = n.knots + 2
int_b = function(x, df = d.f, degree = 2) {
  ibs_x = ibs(x, df, degree = degree)
  ibs_0 = predict(ibs_x, 0)
  sweep(ibs_x, 2, ibs_0, "-")
}
int1_b = function(x.new, x.ref, df = d.f, degree = 2){
  ibs_ref = ibs(x.ref, df = df, degree = degree)
  ibs_new = predict(ibs_ref, x.new)
  ibs_0   = predict(ibs_ref, 0)
  sweep(ibs_new, 2, ibs_0, "-")
}



n.sim = 1e5
set.seed(1)
sim = runif(n.sim, 0, 2*pi)
sim1 = runif(n.sim, 0, 2*pi)
sim2 = runif(n.sim, 0, 2*pi)

cv_vmsp = numeric(n)

t1 = Sys.time()
for(i in 1:n){
  
  test.idx = i
  train.idx = setdiff(1:n, i)
  
  # Split predictors and response
  x.train = x[train.idx, , drop = FALSE]
  x.test  = x[test.idx, , drop = FALSE]
  
  y.train = y[train.idx]
  y.test  = y[test.idx]
 
  fit = vmsp(x.train, y.train, m=200)
  b = fit[3:5]
  g = fit[-c(1:5)]
  
  xp.train = x_pr(b, x.train)
  xp.test  = x_pr(b, x.test)
  
  xp.test = pmin(pmax(xp.test, min(xp.train)),max(xp.train))
  
  if(length(unique(round(xp.train, 8))) <= 3){
    pred4 = NA
  } else {
  
  B_0 = try(int1_b(xp.test, xp.train))
  if(inherits(B_0, "try-error")){
    pred = NA
  } else {
    f_0 = sum(g*B_0[1, ])
    pred = (fit[1] + 2*atan(f_0)) %% (2*pi)
  }
  }
  
  if(all(is.na(pred))){
    cv_vmsp[i] = NA
  } else {
    crp1 = 2*pi * mean(cdist(sim, y.test) * (CircStats::dvm(sim, pred, fit[2])))
    crp2 = 4*pi^2 * mean(cdist(sim1, sim2) *(CircStats::dvm(sim1, pred, fit[2])) *(CircStats::dvm(sim2, pred, fit[2])))
    crp = crp1-0.5*crp2
    cv_vmsp[i] = mean(crp)
  }
}
t2 = Sys.time()
t2-t1




# Average CV error
round(mean(cv_vmsp), 4)






















