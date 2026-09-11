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


n.knots = 2
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
    return(theta.0)   
  } else {
    return(NULL)      
  }  
}


zivm1 = function(x, y, n.knots = 3, m = 200){
  
  #Fit Fisher-Lee model
  fl = suppressWarnings(lm.circular(y, x, init = rnorm(ncol(x), 0, 0.3), type='c-l'))

  # Initial values
  p0 = sum(y==0)/length(y)
  beta0 = fl$coefficients
  phi0 = beta0[-1]/beta0[1]
  phi.init = matrix(phi0, nrow = m, ncol = length(phi0), byrow = TRUE)
  gamma.init = matrix(0, nrow = m, ncol = d.f)
  for(i in 1:m){
    gamma.init[i,] = rnorm(d.f, 0, 0.3)
  }
  init.em = cbind(p0, as.numeric(fl$mu), fl$kappa, phi.init, gamma.init)
  
  # Parallel computation 
  ncores = detectCores()
  cl = makeCluster(ncores-1)
  clusterExport(cl,ls(globalenv()),envir = globalenv())
  registerDoParallel(cl)
  t1 = Sys.time()
  res = foreach(i = 1:m, .combine=rbind, .packages = c("splines2","CircStats", "circular"))%dopar%{
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
fit1 = zivm1(x,y)
t2 = Sys.time()
t2-t1
fit1


############################################


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


#N=4
par(mar = c(5, 4.5, 4, 2) + 0.1)
plot(u1, fhat, type="l", lwd=2,xlab=expression(hat(beta)^T * x), ylab=expression(hat(h)(hat(beta)^T * x)), ylim = c(-7.7,1.6))
lines(u1, lower, lty=2)
lines(u1, upper, lty=2)
legend("bottomleft", legend=c("ZIvM 1", "Confidence band for ZIvM 1","Fisher-Lee"), lty=c(1,2,6), lwd=c(2,1,2),cex=1.1,y.intersp=0.5)   



# Plotting Donut-plots of ZIvM 1 and Fisher-Lee

mu_fl = (as.numeric(fl$mu) + 2*atan(as.vector(x %*% fl$coefficients))) %% (2*pi) # Estimated mean direction in Fisher-Lee model

plot_circ_diag = function(y, y_hat) {
  stopifnot(length(y) == length(y_hat))
  delta = (y_hat - y) %% (2*pi) 
  r = 1 + cos(delta)
  x = r * cos(y_hat)
  y_cart = r * sin(y_hat)
  plot(0,0, type="n", asp=1, xlim=c(-2.5,2.5), ylim=c(-2.5,2.5), axes=FALSE, xlab="", ylab="")
  symbols(rep(0,2), rep(0,2), circles=c(1,2), inches=FALSE, add=TRUE, lwd=1)
  points(0,0, pch=19, col="black")
  
  # Split into clockwise vs anticlockwise
  cw  = delta > pi   # clockwise → solid
  acw = delta > 0 & delta <= pi   # anticlockwise → open
  
  points(x[cw],  y_cart[cw],  pch=19, col="blue",cex=2)   # Filled
  points(x[acw], y_cart[acw], pch=1,  col="blue",cex=2)   # Open
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
    delta_new = theta.0[1:4] + solve(t(z) %*% w(theta.0[1:4],z) %*% z) %*% t(z) %*% (a-p(theta.0[1:4],z))
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
  if (diff <= tol) {
    return(theta.0)   
  } else {
    return(NULL)      
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







