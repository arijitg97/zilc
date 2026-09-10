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

# Projected data
x_pr = function(beta, x){
  as.vector(x %*% (beta / l2(beta)))
}

# Beta reparameterization
beta = function(phi) c(1/sqrt(1+t(phi)%*%phi), as.numeric(1/sqrt(1+t(phi)%*%phi))*phi)

# Jacobian matrix of beta with respect to phi
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

p = function(delta, x){
  eta = as.vector(x %*% delta)
  1/(1+exp(-eta))
}

w = function(delta, x){
  diag(p(delta, x)*(1-p(delta, x)))
}

########## ZIvM 2 ###########

# Objective function for updating (phi,gamma) in ZIvM models
E.zivm = function(x, y, phi, gamma, a, mu){
  b = int_b(x_pr(beta(phi), x))
  delta = y-mu-2*atan(as.vector(b %*% gamma))
  return(-sum((1-a) * cos(delta)))
}

# Log-likelihood for ZIvM 2 model
loglik.zivm2 = function(theta, x, z, y){
  q1 = ncol(x)-1
  q2 = ncol(z)
  b = int_b(x_pr(beta(theta[(q2+3):(q2+2+q1)]), x))
  mu = theta[(q2+1)]+2*atan(as.vector(b %*% theta[-c(1:(q2+2+q1))]))
  sum(log(p(theta[1:q2],z)*dvm(y, 0, 1000) + (1-p(theta[1:q2],z))*dvm(y, mu, theta[(q2+2)]) ))
}

# Function returning estimates under the ZIvM 2 model
lc.zivm2 = function(x, z, y, init, tol = 1e-3, max.iter = 200){
  theta.0 = init  # initial for delta,mu,kappa,phi,gamma
  diff = 1
  iter = 0
  q1 = ncol(x)-1
  q2 = ncol(z)
  while(diff > tol && iter < max.iter){
    theta.1 = theta.0
    b = int_b(x_pr(beta(theta.0[(q2+3):(q2+2+q1)]), x))
    mu = theta.0[(q2+1)]+2*atan(as.vector(b %*% theta.0[-c(1:(q2+2+q1))]))
    alpha = p(theta.0[1:q2],z)
    a = (alpha*dvm(y,0,1000))/(alpha*dvm(y,0,1000)+(1-alpha)*dvm(y,mu,theta.0[(q2+2)]))
    init1 = theta.0[(q2+3):(q2+2+q1)]  # initial for phi
    init2 = theta.0[-c(1:(q2+2+q1))]  # initial for gamma
    gamma_new = optim(init2, function(gamma)E.zivm(x, y, init1, gamma, a, theta.0[(q2+1)]), method = "BFGS", control=list(maxit=2000))$par
    phi_new = optim(init1, function(phi)E.zivm(x, y, phi, gamma_new, a, theta.0[(q2+1)]), method = "BFGS", control=list(maxit=2000))$par
    delta_new = theta.0[1:q2] + solve(t(z) %*% w(theta.0[1:q2],z) %*% z) %*% t(z) %*% (a-p(theta.0[1:q2],z))
    delta_new = as.numeric(delta_new)
    b = int_b(x_pr(beta(phi_new), x))
    del = y-2*atan(as.vector(b %*% gamma_new))
    C = weighted.mean(cos(del), 1-a)
    D = weighted.mean(sin(del), 1-a)
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
}

# Function returning estimates under the ZIvM 2 model given multiple random initial values                    
zivm2.em = function(x, z, y, init){
  M = nrow(init)
  theta.em = matrix(0, nrow = M, ncol = ncol(init))
  ll = rep(0, M)
  for(i in 1:M){
    theta.em[i,] = lc.zivm2(x, z, y, init[i,])
    ll[i] = loglik.zivm2(theta.em[i,], x, z, y)
  }
  theta_hat = theta.em[which.max(ll),]
  return(theta_hat)
}                


########## ZIJP 2 #############


# Normalizing constant of Jones-Pewsey distribution
JPNCon = function(kappa, psi) {
  if (kappa < 0.001) {ncon = 1/(2*pi) ; return(ncon) }
  else {
    eps = 1e-6
    if (abs(psi) <= eps) { ncon = 1/(2*pi*I.0(kappa)) ; return(ncon) }
    else {
      intgrnd = function(x) { (cosh(kappa*psi)+sinh(kappa*psi)*cos(x))**(1/psi) }
      ncon = 1/integrate(intgrnd, lower=-pi, upper=pi)$value
      return(ncon) } }
}

# Density of Jones-Pewsey distribution
djp = function(theta, mu, kappa, psi, ncon) {
  if (kappa < 0.001) {pdfval = 1/(2*pi) ; return(pdfval)}
  else {
    eps = 1e-6
    if (abs(psi) <= eps) {
      pdfval = ncon*exp(kappa*cos(theta-mu)) ; return(pdfval) }
    else {
      pdfval = (cosh(kappa*psi)+sinh(kappa*psi)*cos(theta-mu))**(1/psi)
      pdfval = ncon*pdfval ; return(pdfval) } }
}

# Random number generation from Jones-Pewsey distribution
jpsim = function(n, mu, kappa, psi, ncon) {
  fmax = djp(mu, mu, kappa, psi, ncon) ; theta = 0
  for (j in 1:n) { stopgo = 0
  while (stopgo == 0) {
    u1 = runif(1, 0, 2*pi) ; pdfu1 = djp(u1, mu, kappa, psi, ncon)
    u2 = runif(1, 0, fmax)
    if (u2 <= pdfu1) { theta[j] = u1 ; stopgo = 1 }
  } }
  return(theta)
}

# Negative log-likelihood of Zero-inflated Jones-Pewsey (ZIJP) distribution
zijpnll = function(x, y, a, p, phi, gamma) {
  b = int_b(x_pr(beta(phi), x))
  mu = p[1] + 2*atan(as.vector(b %*% gamma))
  kappa = p[2] ; psi = p[3]
  if (abs(kappa*psi) > 10) return(Inf)
  else { ncon = JPNCon(kappa, psi)
  return(-sum((1-a)*log(djp(y, mu, kappa, psi, ncon)))) }  
}

# Objective function for updating (phi,gamma) in ZIJP models
E.zijp = function(x, y, phi, gamma, a, mu, kappa, psi){
  if(abs(kappa*psi) > 10) return(Inf)
  b = int_b(x_pr(beta(phi), x))
  delta = y-mu-2*atan(as.vector(b %*% gamma))
  eps = 1e-6
  if(abs(psi) <= eps){
    return(-sum((1-a) * cos(delta)))
  } else {
    S = sum ((1-a)*log(cosh(kappa*psi) + sinh(kappa*psi)*cos(delta)))
    return(-(1/psi)*S)
  }
}

# Log-likelihood for ZIJP 2 model
loglik.zijp2 = function(theta, x, z, y){
  q1 = ncol(x)-1
  q2 = ncol(z)
  b = int_b(x_pr(beta(theta[(q2+4):(q2+3+q1)]), x))
  mu = theta[(q2+1)]+2*atan(as.vector(b %*% theta[-c(1:(q2+3+q1))]))
  if(abs(theta[(q2+2)]*theta[(q2+3)]) > 10) return(-Inf)
  else {
    ncon = JPNCon(theta[(q2+2)], theta[(q2+3)])
    sum(log(p(theta[1:q2],z)*dvm(y, 0, 1000) + (1-p(theta[1:q2],z))*djp(y, mu, theta[(q2+2)], theta[(q2+3)], ncon)))
  }
}

# Function returning estimates under the ZIJP 2 model
lc.zijp2 = function(x, z, y, init, tol = 1e-3, max.iter = 200){
  theta.0 = init
  diff = 1
  iter = 0
  q1 = ncol(x)-1
  q2 = ncol(z)
  while(diff > tol && iter < max.iter){
    theta.1 = theta.0
    b = int_b(x_pr(beta(theta.0[(q2+4):(q2+3+q1)]), x))
    mu = theta.0[(q2+1)]+2*atan(as.vector(b %*% theta.0[-c(1:(q2+3+q1))]))
    ncon = JPNCon(theta.0[(q2+2)], theta.0[(q2+3)])
    alpha = p(theta.0[1:q2],z)
    a = (alpha*dvm(y,0,1000))/(alpha*dvm(y,0,1000)+(1-alpha)*djp(y,mu,theta.0[(q2+2)],theta.0[(q2+3)], ncon))
    init1 = theta.0[(q2+4):(q2+3+q1)]  # initial for phi
    init2 = theta.0[-c(1:(q2+3+q1))]   # initial for gamma
    gamma_new = optim(init2, function(gamma)E.zijp(x, y, init1, gamma, a, theta.0[(q2+1)], theta.0[(q2+2)], theta.0[(q2+3)]), method = "BFGS", control=list(maxit=2000))$par
    phi_new = optim(init1, function(phi)E.zijp(x, y, phi, gamma_new, a, theta.0[(q2+1)], theta.0[(q2+2)], theta.0[(q2+3)]), method = "BFGS", control=list(maxit=2000))$par
    delta_new = theta.0[1:q2] + solve(t(z) %*% w(theta.0[1:q2],z) %*% z) %*% t(z) %*% (a-p(theta.0[1:q2],z))
    delta_new = as.numeric(delta_new)
    out = optim(theta.0[(q2+1):(q2+3)], function(p) zijpnll(x, y, a, p, phi_new, gamma_new), method = "L-BFGS-B", lower = c(-pi, 0, -Inf), upper = c(pi, Inf, Inf), control=list(maxit=2000))
    mu_new = out$par[1] ; kappa_new = out$par[2] ; psi_new = out$par[3]
    theta.0 = c(delta_new, mu_new, kappa_new, psi_new, phi_new, gamma_new)
    diff = l2(theta.0 - theta.1)
    iter = iter + 1
  }
  if (diff <= tol) {
    return(theta.0)   
  } else {
    return(NULL)      
  }
}

# Function returning estimates under the ZIJP 2 model given multiple random initial values
zijp2.em = function(x, z, y, init){
  M = nrow(init)
  theta.em = matrix(0, nrow = M, ncol = ncol(init))
  ll = rep(0, M)
  for(i in 1:M){
    theta.em[i,] = lc.zijp2(x, z, y, init[i,])
    ll[i] = loglik.zijp2(theta.em[i,], x, z, y)
  }
  theta_hat = theta.em[which.max(ll),]
  return(theta_hat)
} 


####### Simulation #######


n.knots = 3               # Number of interior knots
d.f = n.knots + 2         # Degrees of freedom of the spline function

# Integrated B-spline basis function               
int_b = function(x, df = d.f, degree = 2) {
  ibs_x = ibs(x, df, degree = degree)
  ibs_0 = predict(ibs_x, 0)
  sweep(ibs_x, 2, ibs_0, "-")
}

# True values of parameters
delta0 = -1; delta1 = -0.5; delta2 = 0.5
delta = c(delta0, delta1, delta2)
mu0 = pi/3
kappa = 1.5
psi = -0.75
ncon = JPNCon(kappa, psi)
beta.0 = c(1, 1)/sqrt(2)
phi = beta.0[-1]/beta.0[1]

f0 = function(x) 1/(1+ exp(-3*(x - 0.2)))
f_0 = function(x) f0(x)-f0(0)             # True single-index function

# Random initial values for spline coefficients
m = 20

# Parallel computation of estimates over multiple replications under ZIJP 1 
ncores = detectCores()
t1 = Sys.time()
cl = makeCluster(ncores-1)
registerDoParallel(cl)
res = foreach(j= 1:510, .combine = rbind, .packages = c("circular","CircStats","splines2"), .errorhandling = "remove") %dopar% {
  # Data generation under ZIJP 2
  set.seed(j)
  n = 200
  x1 = rnorm(n, 0, 0.6)
  x2 = runif(n, -1.1, 1.1)
  x = cbind(x1, x2)
  U = runif(n)
  y = rep(0, n)
  mu = mu0 + 2*atan(as.vector(f_0(x %*% beta.0)))
  for (i in 1:n){
    if(U[i]< p){
      y[i] = rvm(1, 0, 1000)
    }
    else{
      y[i] = jpsim(1, mu[i], kappa, psi, ncon)
    }
  }
  gamma.init = matrix(rnorm(m * d.f, 0, 0.3), m, d.f)
  init1 = cbind(delta0, delta1, delta2, mu0, kappa, psi, phi, gamma.init)
  zijp2.em(x, z, y, init1)   # ZIJP 2 fit
  #init2 = cbind(delta0, delta1, delta2, mu0, kappa, phi, gamma.init)
  #zivm2.em(x, z, y, init2)  # ZIvM 2 fit
}
stopCluster(cl)
t2 = Sys.time()
t2 - t1

del = 501:nrow(res)         # Considering 500 replications
res.1 = res[-del,]
beta_hat = unname(t(sapply(res.1[,7], beta)))         # ZIJP 2
res1 = cbind(res.1[,1:6], beta_hat, res.1[,-c(1:7)])  # ZIJP 2
#beta_hat = unname(t(sapply(res.1[,6], beta)))         # ZIvM 2
#res1 = cbind(res.1[,1:5], beta_hat, res.1[,-c(1:6)])  # ZIvM 2

int1_b = function(x.new, x.ref, df = d.f, degree = 2){
  ibs_ref = ibs(x.ref, df = df, degree = degree)
  ibs_new = predict(ibs_ref, x.new)
  ibs_0   = predict(ibs_ref, 0)
  sweep(ibs_new, 2, ibs_0, "-")
}

#Equispaced points between -1 and 1 for plotting the average estimated single-index function               
N = 100
t = seq(-1, 1, length = N)

ncores = detectCores()
cl = makeCluster(ncores - 1)
registerDoParallel(cl)
X1 = foreach(k = 1:nrow(res1), .combine = rbind, .packages = "splines2") %dopar% {
  
  ## regenerate training covariates
  set.seed(k)
  n = 200
  x1 = rnorm(n, 0, 0.6)
  x2 = runif(n, -1.1, 1.1)
  x  = cbind(x1, x2)
  
  ## estimated beta
  b = res1[k, 7:8]  # ZIJP 2
  #b = res1[k, 6:7]  # ZIvM 2
  
  ## projected training values
  eta = x_pr(b, x)
  
  ## evaluate fitted basis on common grid
  B = int1_b(t, eta)
  
  ## estimated spline coefficients
  g = res1[k, -c(1:8)]  # ZIJP 2
  #g = res1[k, -c(1:7)]  # ZIvM 2
  
  ## fitted function
  hhat = as.vector(B %*% g)
  
  mu.true = (mu0 + 2 * atan(f_0(t))) %% (2*pi)
  mu.hat  = (res1[k,4] + 2 * atan(hhat)) %% (2*pi)
  mspe = mean(1 - cos(mu.hat - mu.true))            # Circular MSPE
  c(mspe, hhat)
}
stopCluster(cl)

cMSPE = mean(X1[,1])  # Average circular MSPE

# Plotting the average estimated single-index function with 95% pointwise confidence interval
mean.fun = colMeans(X1[,-1])
par(mar = c(5, 4.5, 4, 2) + 0.1)
plot(t, mean.fun, type = "l",xlab = expression(hat(beta)^T*x), ylab = expression(bar(h)(hat(beta)^T*x)), lwd=2, ylim=c(-0.83,1.03))
curve(f_0, add = T, col = "red", lwd = 2)
lower = apply(X1[, -1], 2, quantile, probs = 0.025)
upper = apply(X1[, -1], 2, quantile, probs = 0.975)
lines(t, lower, lty = 2)
lines(t, upper, lty = 2)

# Parameter estimates (with their standard deviations) and average circular MSPE
est = data.frame(cbind(apply(res1[,c(1:8)], 2, mean), apply(res1[,c(1:8)], 2, sd)))   # ZIJP 2
#est = data.frame(cbind(apply(res1[,c(1:7)], 2, mean), apply(res1[,c(1:7)], 2, sd)))   # ZIvM 2
est[4,] = c(circ.mean(res1[,4]), sd.circular(res1[,4]))   
colnames(est) = c("Estimate", "Std. Error")
rownames(est)[1:6] = c("delta0", "delta1", "delta2", "mu", "kappa", "psi")  # ZIJP 2
#rownames(est)[1:5] = c("delta0", "delta1", "delta2", "mu", "kappa")         # ZIvM 2
round(est, 4)
round(cMSPE, 4)


























































