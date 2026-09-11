rm(list=ls())

library(circular)
library(CircStats)
library(splines2)
library(numDeriv)

library(parallel)
library(doParallel)
library(foreach)

d1 = read.csv("D:/data/cataract/snare.csv")
d2 = read.csv("D:/data/cataract/vectis.csv")
d3 = read.csv("D:/data/cataract/conv.csv")
#d4 = read.csv("D:/data/cataract/tors.csv")

#getwd()

x11 = na.omit(as.numeric(d1$X.34))[-6]
u1 = na.omit(as.numeric(d1$X.35))[-6]
x12 = sin((4*(pi/180)*u1)%% (2*pi))
x13 = cos((4*(pi/180)*u1)%% (2*pi))
#x12 = sin((pi/180)*u1)
#x13 = cos((pi/180)*u1)
y1 = na.omit(as.numeric(d1$X.43))

x21 = na.omit(as.numeric(d2$X.34))[-15]
u2 = na.omit(as.numeric(d2$X.35))[-15]
x22 = sin((4*(pi/180)*u2)%% (2*pi))
x23 = cos((4*(pi/180)*u2)%% (2*pi))
#x22 = sin((pi/180)*u2)
#x23 = cos((pi/180)*u2)
y2 = na.omit(as.numeric(d2$X.43))

x31 = na.omit(as.numeric(d3$X.34))[-c(7,13,17)]
u3 = na.omit(as.numeric(d3$X.35))[-c(7,13,17)]
x32 = sin((4*(pi/180)*u3)%% (2*pi))
x33 = cos((4*(pi/180)*u3)%% (2*pi))
#x32 = sin((pi/180)*u3)
#x33 = cos((pi/180)*u3)
y3 = na.omit(as.numeric(d3$X.43))


y1 = (4*(pi/180)*y1) %% (2*pi)
y2 = (4*(pi/180)*y2) %% (2*pi)
y3 = (4*(pi/180)*y3) %% (2*pi)

y = c(y1, y2, y3)  # axis of astigmatism in radians for PO 1 month

u = c(x11, x21, x31)  # intensity of astigmatism for PO 7 days
v = c(x12, x22, x32)  # sin(axis of astigmatism for PO 7 days) 
w = c(x13, x23, x33) 

x0 = cbind(u, v, w)    

n=nrow(x0)

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

E2_opt = function(x, y, phi.init, gamma.init, mu, max.iter = 100, tol = 1e-3){
  
  phi.0 = phi.init
  gamma.0 = gamma.init
  diff = 1
  iter = 0
  
  while(diff > tol && iter <= max.iter) {
    phi.1 = phi.0
    gamma.1 = gamma.0
    
    gamma_new = optim(gamma.0, function(gamma)E(x, y, phi.0, gamma, mu), method = "BFGS", control=list(maxit=2000))$par
    phi_new = optim(phi.0, function(phi)E(x, y, phi, gamma_new, mu), method = "BFGS", control=list(maxit=2000))$par
    
    phi.0 = phi_new
    gamma.0 = gamma_new
    diff = l2(c(phi.0, gamma.0) - c(phi.1, gamma.1))
    iter = iter + 1
  }
  #return(c(phi.0, gamma.0))
  if (diff <= tol) {
    return(c(phi.0, gamma.0))   
  } else {
    return(NULL)      
  }
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
    
    #E_opt = E2_opt(x, y, theta.0[3:4], theta.0[-c(1:4)], theta.0[1])
    #gamma_new = E_opt[-c(1:2)]
    #phi_new = E_opt[1:2]
    
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
  #return(theta.0)
  if (diff <= tol) {
    return(theta.0)   # success: return estimate
  } else {
    return(NULL)      # fail: max.iter reached without convergence
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






############### 1 knot ############################




n.knots = 1
d.f = n.knots + 3
int_b = function(x, df = d.f, degree = 2) {
  ibs_x = ibs(x, df, degree = degree)
  ibs_0 = predict(ibs_x, 0)
  sweep(ibs_x, 2, ibs_0, "-")
}


m = 100
set.seed(123)
mu.init = runif(m, -pi, pi)
kappa.init = rexp(m, 1/20)
#phi.init = runif(m, -1, 1)
#gamma.init = mvrnorm(m, mu = runif(d.f, -1, 1), Sigma = diag(rep(1, d.f)))
phi.init = matrix(0, nrow = m, ncol = ncol(x)-1)
gamma.init = matrix(0, nrow = m, ncol = d.f)
for(i in 1:m){
  #set.seed(i)
  #phi.init[i,] = runif(ncol(x)-1, -1, 1)
  phi.init[i,] = rnorm(ncol(x)-1, 0, 0.3)
  #gamma.init[i,] = runif(d.f, -1, 1)
  gamma.init[i,] = rnorm(d.f, 0, 0.3)
}


#init.em = cbind(p, mu0, kappa, atan2(0, 1), gamma.init)
init.em = unname(cbind(mu.init, kappa.init, phi.init, gamma.init))




ncores = detectCores()
cl = makeCluster(ncores-1)
registerDoParallel(cl)
t1 = Sys.time()
res1 = foreach(i = 1:m, .combine=rbind, .packages = c("splines2","CircStats", "circular"),.errorhandling = "remove")%dopar%{
  #cat(i,"\t")
  theta.em = lc.em2(x, y, init.em[i,])
  ll= loglik(theta.em, x, y)
  c(theta.em, ll)
}
stopCluster(cl)
t2 = Sys.time()
t2 - t1

nrow(res1)
nrow(res1)/m

ll = res1[, ncol(res1)]
theta.em = res1[, -ncol(res1)]
theta_hat1 = theta.em[which.max(ll),]
theta_hat1

max(ll)


fit1 = 2*length(theta_hat1) - 2*loglik(theta_hat1, x, y)
bic1 = length(theta_hat1)*log(length(y)) - 2*loglik(theta_hat1, x, y)
est1 = c(theta_hat1[1:2], beta(theta_hat1[3:4]), theta_hat1[-c(1:4)])


cov1 = solve(numDeriv::hessian(function(theta)-loglik(theta, x, y),theta_hat1))
#cov1 = solve(pracma::hessian(function(theta)-loglik(theta, x, y),theta_hat1))
se1 = sqrt(diag(cov1))


cov1.beta = jac.beta(theta_hat1[3:4])%*%cov1[3:4,3:4]%*%t(jac.beta(theta_hat1[3:4]))
se1.beta = sqrt(diag(cov1.beta))



N = 100
X1 = X2 = X3 = seq(-1, 1, length = N)
X = cbind(X1, X2, X3)

b1 = est1[3:5]
B1 = int_b(x_pr(b1, X))
g1 = est1[-c(1:5)]
f_1 = rep(0, N)
for(i in 1:N){
  f_1[i] = sum(g1*B1[i,])
}
plot(x_pr(b1, X), f_1, type = "l")


r1 = data.frame(cbind(est1, c(se1[1:2], se1.beta, se1[-c(1:4)])))
colnames(r1) = c("Estimate", "Std. Error")
rownames(r1)[1:5] = c("mu", "kappa","beta_1","beta_2","beta_3")

result1 = list(r1, fit1, bic1)
names(result1[[2]])="AIC"
names(result1[[3]])="BIC"

round(result1[[1]], 4)
round(result1[[2]], 4)
round(result1[[3]], 4)


#getwd()

write.csv(result1, file="N=1.csv")
#write.csv(theta_hat1, file="N=01.csv")



B_1 = int_b(x_pr(b1, x))
f_11 = rep(0, nrow(x))
for(i in 1:nrow(x)){
  f_11[i] = sum(g1*B_1[i,])
}

mu1_ZIVM = est1[1] + 2*atan(f_11)
round(sum(1-cos(y[1:37]-mu1_ZIVM[1:37])),4)


##########################################




n.knots = 2
d.f = n.knots + 3
int_b = function(x, df = d.f, degree = 2) {
  ibs_x = ibs(x, df, degree = degree)
  ibs_0 = predict(ibs_x, 0)
  sweep(ibs_x, 2, ibs_0, "-")
}

m = 100
set.seed(123)
mu.init = runif(m, -pi, pi)
kappa.init = rexp(m, 1/20)
#phi.init = runif(m, -1, 1)
#gamma.init = mvrnorm(m, mu = runif(d.f, -1, 1), Sigma = diag(rep(1, d.f)))
phi.init = matrix(0, nrow = m, ncol = ncol(x)-1)
gamma.init = matrix(0, nrow = m, ncol = d.f)
for(i in 1:m){
  #set.seed(i)
  #phi.init[i,] = runif(ncol(x)-1, -1, 1)
  phi.init[i,] = rnorm(ncol(x)-1, 0, 0.3)
  #gamma.init[i,] = runif(d.f, -1, 1)
  gamma.init[i,] = rnorm(d.f, 0, 0.3)
}


#init.em = cbind(p, mu0, kappa, atan2(0, 1), gamma.init)
init.em = unname(cbind(mu.init, kappa.init, phi.init, gamma.init))


ncores = detectCores()
cl = makeCluster(ncores-1)
registerDoParallel(cl)
t1 = Sys.time()
res2 = foreach(i = 1:m, .combine=rbind, .packages = c("splines2","CircStats", "circular"),.errorhandling = "remove")%dopar%{
  #cat(i,"\t")
  theta.em = lc.em2(x, y, init.em[i,])
  ll= loglik(theta.em, x, y)
  c(theta.em, ll)
}
stopCluster(cl)
t2 = Sys.time()
t2 - t1

nrow(res2)
nrow(res2)/m

ll = res2[, ncol(res2)]
theta.em = res2[, -ncol(res2)]
theta_hat2 = theta.em[which.max(ll),]
theta_hat2

max(ll)


fit2 = 2*length(theta_hat2) - 2*loglik(theta_hat2, x, y)
bic2 = length(theta_hat2)*log(length(y)) - 2*loglik(theta_hat2, x, y)
est2 = c(theta_hat2[1:2], beta(theta_hat2[3:4]), theta_hat2[-c(1:4)])


cov2 = solve(numDeriv::hessian(function(theta)-loglik(theta, x, y),theta_hat2))
#cov2 = solve(pracma::hessian(function(theta)-loglik(theta, x, y),theta_hat2))
se2 = sqrt(diag(cov2))


cov2.beta = jac.beta(theta_hat2[3:4])%*%cov2[3:4,3:4]%*%t(jac.beta(theta_hat2[3:4]))
se2.beta = sqrt(diag(cov2.beta))


N = 100
X1 = X2 = X3 = seq(-1, 1, length = N)
X = cbind(X1, X2, X3)
b2 = est2[3:5]
B2 = int_b(x_pr(b2, X))
g2 = est2[-c(1:5)]
f_2 = rep(0, N)
for(i in 1:N){
  f_2[i] = sum(g2*B2[i,])
}
plot(x_pr(b2, X), f_2, type = "l",lwd=2, xlab=expression(beta^T * X), ylab=expression(h(beta^T * X)))
lines(x_pr(b2, X), conf.lower, lty=2)
lines(x_pr(b2, X), conf.upper, lty=2)


V = cov2[-c(1:4),-c(1:4)]
#sqrt(diag(V9))

var.h = conf.lower = conf.upper = rep(0, N)
for(i in 1:N){
  var.h[i] = t(B2[i,])%*%V%*%B2[i,]
  conf.lower[i] = t(g2)%*%B2[i,] - 1.96*sqrt(var.h[i])
  conf.upper[i] = t(g2)%*%B2[i,] + 1.96*sqrt(var.h[i])
}




r2 = data.frame(cbind(est2, c(se2[1:2], se2.beta, se2[-c(1:4)])))
colnames(r2) = c("Estimate", "Std. Error")
rownames(r2)[1:5] = c("mu", "kappa","beta_1","beta_2","beta_3")

result2 = list(r2, fit2, bic2)
names(result2[[2]])="AIC"
names(result2[[3]])="BIC"

round(result2[[1]], 4)
round(result2[[2]], 4)
round(result2[[3]], 4)


#write.csv(result2, file="N=2.csv")


B_2 = int_b(x_pr(b2, x))
f_21 = rep(0, nrow(x))
for(i in 1:nrow(x)){
  f_21[i] = sum(g2*B_2[i,])
}
mu_ZIVM_2 = est2[1] + 2*atan(f_21)
round(sum(1-cos(y[1:37]-mu_ZIVM_2[1:37])),4)




##############################################################




n.knots = 3
d.f = n.knots + 3
int_b = function(x, df = d.f, degree = 2) {
  ibs_x = ibs(x, df, degree = degree)
  ibs_0 = predict(ibs_x, 0)
  sweep(ibs_x, 2, ibs_0, "-")
}

m = 100
set.seed(123)
mu.init = runif(m, -pi, pi)
kappa.init = rexp(m, 1/20)
#phi.init = runif(m, -1, 1)
#gamma.init = mvrnorm(m, mu = runif(d.f, -1, 1), Sigma = diag(rep(1, d.f)))
phi.init = matrix(0, nrow = m, ncol = ncol(x)-1)
gamma.init = matrix(0, nrow = m, ncol = d.f)
for(i in 1:m){
  #set.seed(i)
  #phi.init[i,] = runif(ncol(x)-1, -1, 1)
  phi.init[i,] = rnorm(ncol(x)-1, 0, 0.3)
  #gamma.init[i,] = runif(d.f, -1, 1)
  gamma.init[i,] = rnorm(d.f, 0, 0.3)
}

init.em = unname(cbind(mu.init, kappa.init, phi.init, gamma.init))

ncores = detectCores()
cl = makeCluster(ncores-1)
registerDoParallel(cl)
t1 = Sys.time()
res3 = foreach(i = 1:m, .combine=rbind, .packages = c("splines2","CircStats", "circular"),.errorhandling = "remove")%dopar%{
  #cat(i,"\t")
  theta.em = lc.em2(x, y, init.em[i,])
  ll= loglik(theta.em, x, y)
  c(theta.em, ll)
}
stopCluster(cl)
t2 = Sys.time()
t2 - t1

nrow(res3)
nrow(res3)/m

ll = res3[, ncol(res3)]
theta.em = res3[, -ncol(res3)]
theta_hat3 = theta.em[which.max(ll),]
theta_hat3


max(ll)

fit3 = 2*length(theta_hat3) - 2*loglik(theta_hat3, x, y)
bic3 = length(theta_hat3)*log(length(y)) - 2*loglik(theta_hat3, x, y)
est3 = c(theta_hat3[1:2], beta(theta_hat3[3:4]), theta_hat3[-c(1:4)])


cov3 = solve(numDeriv::hessian(function(theta)-loglik(theta, x, y),theta_hat3))
#cov3 = solve(pracma::hessian(function(theta)-loglik(theta, x, y),theta_hat1))
se3 = sqrt(diag(cov3))


cov3.beta = jac.beta(theta_hat3[3:4])%*%cov3[3:4,3:4]%*%t(jac.beta(theta_hat3[3:4]))
se3.beta = sqrt(diag(cov3.beta))



N = 100
X1 = X2 = X3 = seq(-1, 1, length = N)
X = cbind(X1, X2, X3)

b3 = est3[3:5]
B3 = int_b(x_pr(b3, X))
g3 = est3[-c(1:5)]
f_3 = rep(0, N)
for(i in 1:N){
  f_3[i] = sum(g3*B3[i,])
}
plot(x_pr(b3, X), f_3, type = "l")


r3 = data.frame(cbind(est3, c(se3[1:2], se3.beta, se3[-c(1:4)])))
colnames(r3) = c("Estimate", "Std. Error")
rownames(r3)[1:5] = c("mu", "kappa","beta_1","beta_2","beta_3")

result3 = list(r3, fit3, bic3)
names(result3[[2]])="AIC"
names(result3[[3]])="BIC"

round(result3[[1]], 4)
round(result3[[2]], 4)
round(result3[[3]], 4)


#getwd()

write.csv(result1, file="N=1.csv")
#write.csv(theta_hat1, file="N=01.csv")



B_3 = int_b(x_pr(b3, x))
f_31 = rep(0, nrow(x))
for(i in 1:nrow(x)){
  f_31[i] = sum(g3*B_3[i,])
}
mu3_ZIVM = (est3[1] + 2*atan(f_31)) %% (2*pi)
round(sum(1-cos(y[1:37]-mu3_ZIVM[1:37])),4)
#5.5578


round(sum(1-cos(y-mu3_ZIVM)),4)
#5.9765




########################################################


n.knots = 2
d.f = n.knots + 3
int_b = function(x, df = d.f, degree = 2) {
  ibs_x = ibs(x, df, degree = degree)
  ibs_0 = predict(ibs_x, 0)
  sweep(ibs_x, 2, ibs_0, "-")
}

m = 200
#set.seed(1)
mu.init = runif(m, 0.3296243-0.3, 0.3296243+0.3)
#mu0 = 0.3296243
kappa.init = rexp(m, 1/3.700447)
#kappa0 = 3.700447
#gamma.init = mvrnorm(m, mu = runif(d.f, -1, 1), Sigma = diag(rep(1, d.f)))

phi0=c(-3.019547 , 1.739724)
library(MASS)
phi.init = mvrnorm(m, mu = phi0, Sigma = 0.1^2 * diag(ncol(x)-1))
#phi.init = matrix(0, nrow = m, ncol = ncol(x)-1)
gamma.init = matrix(0, nrow = m, ncol = d.f)

for(i in 1:m){
  #set.seed(i)
  #phi.init[i,] = runif(ncol(x)-1, -1, 1)
  #phi.init[i,] = rnorm(ncol(x)-1, 0, 0.3)
  #gamma.init[i,] = runif(d.f, -1, 1)
  gamma.init[i,] = rnorm(d.f, 0, 0.3)
}

init.em = unname(cbind(mu.init, kappa.init, phi.init, gamma.init))

ncores = detectCores()
cl = makeCluster(ncores-1)
registerDoParallel(cl)
t1 = Sys.time()
res4 = foreach(i = 1:m, .combine=rbind, .packages = c("splines2","CircStats", "circular"),.errorhandling = "remove")%dopar%{
  #cat(i,"\t")
  theta.em = lc.em2(x, y, init.em[i,])
  ll= loglik(theta.em, x, y)
  c(theta.em, ll)
}
stopCluster(cl)
t2 = Sys.time()
t2 - t1

nrow(res4)
nrow(res4)/m

ll = res4[, ncol(res4)]
theta.em = res4[, -ncol(res4)]
theta_hat4 = theta.em[which.max(ll),]
theta_hat4


max(ll)

fit4 = 2*length(theta_hat4) - 2*loglik(theta_hat4, x, y)
bic4 = length(theta_hat4)*log(length(y)) - 2*loglik(theta_hat4, x, y)
est4 = c(theta_hat4[1:2], beta(theta_hat4[3:4]), theta_hat4[-c(1:4)])


cov4 = solve(numDeriv::hessian(function(theta)-loglik(theta, x, y),theta_hat4))
#cov4 = solve(pracma::hessian(function(theta)-loglik(theta, x, y),theta_hat4))
se4 = sqrt(diag(cov4))


cov4.beta = jac.beta(theta_hat4[3:4])%*%cov4[3:4,3:4]%*%t(jac.beta(theta_hat4[3:4]))
se4.beta = sqrt(diag(cov4.beta))



N = 100
X1 = X2 = X3 = seq(-1, 1, length = N)
X = cbind(X1, X2, X3)

b4 = est4[3:5]
B4 = int_b(x_pr(b4, X))
g4 = est4[-c(1:5)]
f_4 = rep(0, N)
for(i in 1:N){
  f_4[i] = sum(g4*B4[i,])
}
plot(x_pr(b4, X), f_4, type = "l")


r4 = data.frame(cbind(est4, c(se4[1:2], se4.beta, se4[-c(1:4)])))
colnames(r4) = c("Estimate", "Std. Error")
rownames(r4)[1:5] = c("mu", "kappa","beta_1","beta_2","beta_3")

result4 = list(r4, fit4, bic4)
names(result4[[2]])="AIC"
names(result4[[3]])="BIC"

round(result4[[1]], 4)
round(result4[[2]], 4)
round(result4[[3]], 4)


#getwd()

#write.csv(result1, file="N=4.csv")
#write.csv(theta_hat1, file="N=01.csv")



B_4 = int_b(x_pr(b4, x))
f_41 = rep(0, nrow(x))
for(i in 1:nrow(x)){
  f_41[i] = sum(g4*B_4[i,])
}
mu4_ZIVM = (est4[1] + 2*atan(f_41)) %% (2*pi)
#round(sum(1-cos(y[1:37]-mu4_ZIVM[1:37])),4)
# 5.5436


#round(sum(1-cos(y-mu4_ZIVM)),4)




########################################################




n.knots = 5
d.f = n.knots + 3
int_b = function(x, df = d.f, degree = 2) {
  ibs_x = ibs(x, df, degree = degree)
  ibs_0 = predict(ibs_x, 0)
  sweep(ibs_x, 2, ibs_0, "-")
}

m = 100
#set.seed(123)
mu.init = runif(m, -pi, pi)
kappa.init = rexp(m, 1/20)
#phi.init = runif(m, -1, 1)
#gamma.init = mvrnorm(m, mu = runif(d.f, -1, 1), Sigma = diag(rep(1, d.f)))
phi.init = matrix(0, nrow = m, ncol = ncol(x)-1)
gamma.init = matrix(0, nrow = m, ncol = d.f)
for(i in 1:m){
  #set.seed(i)
  #phi.init[i,] = runif(ncol(x)-1, -1, 1)
  phi.init[i,] = rnorm(ncol(x)-1, 0, 0.3)
  #gamma.init[i,] = runif(d.f, -1, 1)
  gamma.init[i,] = rnorm(d.f, 0, 0.3)
}

init.em = unname(cbind(mu.init, kappa.init, phi.init, gamma.init))

ncores = detectCores()
cl = makeCluster(ncores-1)
registerDoParallel(cl)
t1 = Sys.time()
res5 = foreach(i = 1:m, .combine=rbind, .packages = c("splines2","CircStats", "circular"),.errorhandling = "remove")%dopar%{
  #cat(i,"\t")
  theta.em = lc.em2(x, y, init.em[i,])
  ll= loglik(theta.em, x, y)
  c(theta.em, ll)
}
stopCluster(cl)
t2 = Sys.time()
t2 - t1

nrow(res5)
nrow(res5)/m

ll = res5[, ncol(res5)]
theta.em = res5[, -ncol(res5)]
theta_hat5 = theta.em[which.max(ll),]
theta_hat5


max(ll)

fit5 = 2*length(theta_hat5) - 2*loglik(theta_hat5, x, y)
bic5 = length(theta_hat5)*log(length(y)) - 2*loglik(theta_hat5, x, y)
est5 = c(theta_hat5[1:2], beta(theta_hat5[3:4]), theta_hat5[-c(1:4)])


cov5 = solve(numDeriv::hessian(function(theta)-loglik(theta, x, y),theta_hat5))
#cov5 = solve(pracma::hessian(function(theta)-loglik(theta, x, y),theta_hat5))
se5 = sqrt(diag(cov5))


cov5.beta = jac.beta(theta_hat5[3:4])%*%cov5[3:4,3:4]%*%t(jac.beta(theta_hat5[3:4]))
se5.beta = sqrt(diag(cov5.beta))



N = 100
X1 = X2 = X3 = seq(-1, 1, length = N)
X = cbind(X1, X2, X3)

b5 = est5[3:5]
B5 = int_b(x_pr(b5, X))
g5 = est5[-c(1:5)]
f_5 = rep(0, N)
for(i in 1:N){
  f_5[i] = sum(g5*B5[i,])
}
plot(x_pr(b5, X), f_5, type = "l")


r5 = data.frame(cbind(est5, c(se5[1:2], se5.beta, se5[-c(1:4)])))
colnames(r5) = c("Estimate", "Std. Error")
rownames(r5)[1:5] = c("mu", "kappa","beta_1","beta_2","beta_3")

result5 = list(r5, fit5, bic5)
names(result5[[2]])="AIC"
names(result5[[3]])="BIC"

round(result5[[1]], 4)
round(result5[[2]], 4)
round(result5[[3]], 4)


#getwd()

write.csv(result1, file="N=5.csv")
#write.csv(theta_hat1, file="N=01.csv")



B_5 = int_b(x_pr(b5, x))
f_51 = rep(0, nrow(x))
for(i in 1:nrow(x)){
  f_51[i] = sum(g5*B_5[i,])
}
mu5_ZIVM = est5[1] + 2*atan(f_51)
round(sum(1-cos(y[1:37]-mu5_ZIVM[1:37])),4)
# 5.5436


round(sum(1-cos(y-mu5_ZIVM)),4)






########################################################





cdist = function(a, b){
  delta = abs(a - b)
  pmin(delta, 2*pi - delta)
}
round(sum(cdist(y, mu4_ZIVM)), 4)
#18.3781

round(sum(cdist(y[1:37], mu4_ZIVM[1:37])), 4)
#16.7549


round(sum(cdist(y[1:37], mu3_ZIVM[1:37])), 4)
#16.5089


#residual vs predictors(projected along beta)
plot(x_pr(b3,x[1:37,]),1-cos(y[1:37]-mu3_ZIVM[1:37]))

#qq-plot of observed data vs predicted data
#qqplot(y,mu3_ZIVM)
#qqline(y)




##########################################

#N=4
n.sim = 1e6
set.seed(1)
sim = runif(n.sim, 0, 2*pi)
sim1 = runif(n.sim, 0, 2*pi)
sim2 = runif(n.sim, 0, 2*pi)

crp1 = sapply(1:n, function(i){
  2*pi * mean(cdist(sim, y[i]) * (CircStats::dvm(sim,mu4_ZIVM[i],est4[2])))
})
crp2 = sapply(1:n, function(i){
  4*pi^2 * mean(cdist(sim1, sim2) *
                  (CircStats::dvm(sim1,mu4_ZIVM[i],est4[2])) *
                  (CircStats::dvm(sim2,mu4_ZIVM[i],est4[2])))
})
crp = crp1-0.5*crp2
round(mean(crp),4)
#0.2601




#0.2611 ; N=4
#0.2665 ; N=3



########################################


x1 = runif(n,min(u),max(u)); x2=runif(n,min(v),max(v)); x3=runif(n,min(w),max(w))
x_hat = apply(cbind(x1,x2,x3), 2, x_centre)
B0 = int_b(x_pr(b3, x_hat))
f0 = rep(0, nrow(x_hat))
for(i in 1:nrow(x_hat)){
  f0[i] = sum(g3*B0[i,])
}
mu_hat = est3[1]+2*atan(f0)
y_hat = rep(0,n)
for(i in 1:n){
  y_hat[i] = rvm(1,mu_hat[i], est3[2])
}
qqplot(y, y_hat)


y_hat = rep(0,n)
for(i in 1:n){
  y_hat[i] = rvm(1,mu3_ZIVM[i], est3[2])
}
qqplot(y, y_hat)
abline(0,1)



r_in <- 1.25
r_out <- 2.5

# Cartesian coordinates
x_in <- r_in * cos(y_hat)
y_in <- r_in * sin(y_hat)
x_out <- r_out * cos(y)
y_out <- r_out * sin(y)

# Set up plot
plot(0,0, type="n", asp=1, xlim=c(-2.5,2.5), ylim=c(-2.5,2.5), axes=FALSE, xlab="", ylab="")
symbols(rep(0,2), rep(0,2), circles=c(r_in, r_out), inches=FALSE, add=TRUE, lwd=1)

# Add dot at the center
points(0, 0, pch=19, col="black", cex=1.5)

# Draw edges (connect inner to outer)
for(i in 1:n){
  segments(x_in[i], y_in[i], x_out[i], y_out[i], col="black")
}

# Add points on circles
points(x_in, y_in, pch=19, col="blue", cex = 1.5)
points(x_out, y_out, pch=19, col="blue", cex = 1.5)

###############################


plot_circ_diag <- function(y, y_hat) {
  stopifnot(length(y) == length(y_hat))
  
  delta <- (y_hat - y ) %% (2*pi) 
  
  # radius
  r <- 1 + cos(delta)
  
  # Cartesian coordinates (use predicted angle)
  x <- r * cos(y_hat)
  y_cart <- r * sin(y_hat)
  
  # setup plot
  plot(0,0, type="n", asp=1, xlim=c(-2.5,2.5), ylim=c(-2.5,2.5),
       axes=FALSE, xlab="", ylab="")
  symbols(rep(0,2), rep(0,2), circles=c(1,2), inches=FALSE, add=TRUE, lwd=1)
  
  # center dot
  points(0,0, pch=19, col="black")
  
  # split into clockwise vs anticlockwise
  cw  <- delta > pi   # clockwise → solid
  acw <- delta > 0 & delta <= pi   # anticlockwise → open
  
  # plot points
  points(x[cw],  y_cart[cw],  pch=19, col="blue")   # filled
  points(x[acw], y_cart[acw], pch=1,  col="blue")   # open
}

#plot_circ_diag(y,y_hat)
plot_circ_diag(y[1:37], mu4_ZIVM[1:37])


#################################

qqplot(y, mu4_ZIVM %% (2*pi))
abline(0,1)


# non-zero y
plot_circ_diag(y[1:37], mu4_ZIVM[1:37])
qqplot(y[1:37], mu4_ZIVM[1:37] %% (2*pi))
abline(0,1)

##############################




vmsp = function(x, y, m = 200){
  #Fit Fisher-Lee model
  fl = suppressWarnings(lm.circular(y, x, init = rnorm(ncol(x), 0, 0.3), type='c-l'))
  
  #d.f = n.knots + 3
  #n = length(y)
  
  # Initial values
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
  #init = unname(cbind(mu.init, kappa.init, phi.init, gamma.init))
  init = cbind(as.numeric(fl$mu), fl$kappa, phi.init, gamma.init)
  
  
  # Parallel computation 
  ncores = detectCores()
  cl = makeCluster(ncores-1)
  clusterExport(cl,ls(globalenv()),envir = globalenv())
  registerDoParallel(cl)
  #t1 = Sys.time()
  res = foreach(i = 1:m, .combine=rbind, .packages = c("splines2","CircStats", "circular"))%dopar%{
    #cat(i,"\t")
    theta.m = lc.vmsp(x, y, init[i,])
    ll= loglik.vmsp(theta.m, x, y)
    c(theta.m, ll)
  }
  stopCluster(cl)
  #t2 = Sys.time()
  #t2 - t1
  
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





#N=1
$Estimate
mu   kappa  beta_1  beta_2  beta_3 
0.3172  3.8490  0.2432 -0.7568  0.6068 

$Std.Error
mu  kappa beta_1 beta_2 beta_3 
0.1463 0.6640 0.0333 0.0438 0.0576 

$AIC
[1] 103.3984

$BIC
[1] 117.3213

$Sum_Residual
[1] 18.6585

$Avg.CRPS
[1] 0.2993




#N=2
$Estimate
mu   kappa  beta_1  beta_2  beta_3 
0.3698  3.8695  0.2289 -0.7949  0.5620 

$Std.Error
mu  kappa beta_1 beta_2 beta_3 
0.1634 0.6681 0.0388    NaN 0.0846 

$AIC
[1] 105.0497

$BIC
[1] 120.9615

$Sum_Residual
[1] 18.453

$Avg.CRPS
[1] 0.2975



#N=3
$Estimate
mu   kappa  beta_1  beta_2  beta_3 
0.5967  4.2951  0.1900 -0.8443  0.5011 

$Std.Error
mu  kappa beta_1 beta_2 beta_3 
0.1818 0.7529 0.0334 0.1152 0.0870 

$AIC
[1] 100.303

$BIC
[1] 118.2039

$Sum_Residual
[1] 17.3295

$Avg.CRPS
[1] 0.2807




#N=4
$Estimate
mu   kappa  beta_1  beta_2  beta_3 
0.6011  4.5021  0.1729 -0.8563  0.4867 

$Std.Error
mu  kappa beta_1 beta_2 beta_3 
0.1970 0.7942 0.0231 0.1098 0.0661 

$AIC
[1] 99.3166

$BIC
[1] 119.2064

$Sum_Residual
[1] 16.3896

$Avg.CRPS
[1] 0.2719



#N=5
$Estimate
mu   kappa  beta_1  beta_2  beta_3 
0.6930  4.7163  0.2249 -0.8615  0.4552 

$Std.Error
mu  kappa beta_1 beta_2 beta_3 
0.1998 0.8371 0.0274    NaN 0.0660 

$AIC
[1] 98.3976

$BIC
[1] 120.2764

$Sum_Residual
[1] 17.017

$Avg.CRPS
[1] 0.2654





####### LOOCV #########


n.knots = 1
d.f = n.knots + 3
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
#0.3776 (n.knots=4)


round(mean(cv_vmsp, na.rm = TRUE), 4)




#Updated

#0.3692 (n.knots = 0)
#0.4028 (n.knots = 1)
#0.3724 (n.knots = 2)
#0.4052 (n.knots = 3)
#0.3773 (n.knots = 4)
#0.39 (n.knots = 5)


















