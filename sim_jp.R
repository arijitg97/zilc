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



##########################################



JPNCon = function(kappa, psi) {
  if (kappa < 0.001) {ncon = 1/(2*pi) ; return(ncon) }
  else {
    eps=1e-6
    if (abs(psi) <= eps) { ncon = 1/(2*pi*I.0(kappa)) ; return(ncon) }
    else {
      intgrnd = function(x) { (cosh(kappa*psi)+sinh(kappa*psi)*cos(x))**(1/psi) }
      ncon = 1/integrate(intgrnd, lower=-pi, upper=pi)$value
      return(ncon) } }
}


djp = function(theta, mu, kappa, psi, ncon) {
  if (kappa < 0.001) {pdfval = 1/(2*pi) ; return(pdfval)}
  else {
    eps=1e-6
    if (abs(psi) <= eps) {
      pdfval = ncon*exp(kappa*cos(theta-mu)) ; return(pdfval) }
    else {
      pdfval = (cosh(kappa*psi)+sinh(kappa*psi)*cos(theta-mu))**(1/psi)
      pdfval = ncon*pdfval ; return(pdfval) } }
}


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






zijpnll = function(x, y, a, p, phi, gamma) {
  b = int_b(x_pr(beta(phi), x))
  mu = p[1] + 2*atan(as.vector(b %*% gamma))
  kappa = p[2] ; psi = p[3]
  if (abs(kappa*psi) > 10) return(Inf)
  else { ncon = JPNCon(kappa, psi)
  return(-sum((1-a)*log(djp(y, mu, kappa, psi, ncon)))) }  
}


E.zijp1 = function(x, y, phi, gamma, a, mu, kappa, psi){
  
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


loglik.zijp1 = function(theta, x, y){
  q = ncol(x)-1
  b = int_b(x_pr(beta(theta[5:(4+q)]), x))
  mu = theta[2]+2*atan(as.vector(b %*% theta[-c(1:(4+q))]))
  if(abs(theta[3]*theta[4]) > 10) return(-Inf)
  else {
    ncon = JPNCon(theta[3], theta[4])
    sum(log(theta[1]*dvm(y, 0, 1000) + (1-theta[1])*djp(y, mu, theta[3], theta[4], ncon)))
  }
}


lc.zijp1 = function(x, y, init, tol = 1e-3, max.iter = 200){
  
  theta.0 = init
  diff = 1
  iter = 0
  q = ncol(x)-1
  
  while(diff > tol && iter < max.iter){
    theta.1 = theta.0
    
    b = int_b(x_pr(beta(theta.0[5:(4+q)]), x))
    mu = theta.0[2]+2*atan(as.vector(b %*% theta.0[-c(1:(4+q))]))
    
    ncon = JPNCon(theta.0[3], theta.0[4])
    a = (theta.0[1]*dvm(y,0,1000))/(theta.0[1]*dvm(y,0,1000)+(1-theta.0[1])*djp(y,mu,theta.0[3],theta.0[4], ncon))
    
    gamma_new = optim(theta.0[-c(1:(4+q))], function(gamma)E.zijp1(x, y, theta.0[5:(4+q)], gamma, a, theta.0[2], theta.0[3], theta.0[4]), method = "BFGS", control=list(maxit=2000))$par
    phi_new = optim(theta.0[5:(4+q)], function(phi)E.zijp1(x, y, phi, gamma_new, a, theta.0[2], theta.0[3], theta.0[4]), method = "BFGS", control=list(maxit=2000))$par
    
    out = optim(theta.0[2:4], function(p) zijpnll(x, y, a, p, phi_new, gamma_new), method = "L-BFGS-B", lower = c(-pi, 0, -Inf), upper = c(pi, Inf, Inf), control=list(maxit=2000))
    mu_new = out$par[1] ; kappa_new = out$par[2] ; psi_new = out$par[3]
    
    p_new = mean(a)
    
    theta.0 = c(p_new, mu_new, kappa_new, psi_new, phi_new, gamma_new)
    diff = l2(theta.0 - theta.1)
    #diff = l2(theta.0[-c(1:3)] - theta.1[-c(1:3)])
    iter = iter + 1
    #cat(iter,diff,loglik.zijp1(theta.0, x, y),"\n")
    
  }
  #return(theta.0)
  if (diff <= tol) {
    return(theta.0)   
  } else {
    return(NULL)      
  }
}



####### Simulation #######


n.knots = 1
d.f = n.knots + 3
int_b = function(x, df = d.f, degree = 2) {
  ibs_x = ibs(x, df, degree = degree)
  ibs_0 = predict(ibs_x, 0)
  sweep(ibs_x, 2, ibs_0, "-")
}




lc.em = function(x, y, init){
  M = nrow(init)
  theta.em = matrix(0, nrow = M, ncol = 9)
  ll = rep(0, M)
  for(i in 1:M){
    theta.em[i,] = lc.zijp1(x, y, init[i,])
    ll[i] = loglik.zijp1(theta.em[i,], x, y)
  }
  theta_hat = theta.em[which.max(ll),]
  return(theta_hat)
}  


#true values
p = 0.1
mu0 = pi/4
kappa = 2.5
psi = 1.75
ncon = JPNCon(kappa, psi)
#beta.0 = c(1/sqrt(3), sqrt(2/3))
beta.0 = c(1/2, sqrt(3)/2)
phi = beta.0[-1]/beta.0[1]

m = 20
gamma.init = matrix(0, nrow = m, ncol = d.f)
for(i in 1:m){
  gamma.init[i,] = rnorm(d.f, 0, 0.3)
}
init1 = cbind(p, mu0, kappa, psi, phi, gamma.init)


ncores = detectCores()
t1 = Sys.time()
cl = makeCluster(ncores-1)
registerDoParallel(cl)
r = foreach(j= 1:1050, .combine = rbind, .packages = c("circular","CircStats","splines2"), .errorhandling = "remove") %dopar% {
  
  set.seed(j)
  n = 200
  x1 = rnorm(n, 1, 2)
  x2 = runif(n, -1, 1)
  x = cbind(x1, x2)
  U = runif(n)
  y = rep(0, n)
  mu = mu0 + 2*atan(as.vector(2*(x %*% beta.0)^2))
  #mu = mu0 + 2*atan(as.vector(x %*% beta.0))
  for (i in 1:n){
    #mu[i] = mu0 + 2*atan(sin(5*x1[i])/2)
    if(U[i]< p){
      y[i] = rvm(1, 0, 1000)
    }
    else{
      #y[i] = rvm(1, mu[i], kappa)
      y[i] = jpsim(1, mu[i], kappa, psi, ncon)
    }
  }
  
  lc.em(x, y, init1)
}
stopCluster(cl)
t2 = Sys.time()
t2 - t1


apply(r, 2, mean)




































