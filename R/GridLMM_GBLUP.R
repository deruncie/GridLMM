library(GridLMM)
library(rrBLUP)
GridLMM_GBLUP = function(Knn,K_no,...) {
  
}


no = 500
nn = 100
p = 1000
Xo = matrix(rbinom(no*p,1,rep(.5,2)),no,p)
Xn = matrix(rbinom(nn*p,1,rep(.5,2)),nn,p)
Xom = colMeans(Xo)
Xo = sweep(Xo,2,Xom,'-')
Xn = sweep(Xn,2,Xom,'-')
Xn[1:50,1:(p/2)] = 0
b = rnorm(p)

h2 = .5
s2 = 1
uo = Xo %*% b
un = Xn %*% b

un = un/sd(uo) * sqrt(h2*s2)
uo = uo/sd(uo) * sqrt(h2*s2)
eo = rnorm(no,0,sqrt((1-h2)*s2))
yo = uo + eo
yo2 = uo + 2*eo
yo = yo-mean(yo)


Ko = tcrossprod(Xo)/p*2
# Kn = tcrossprod(Xn,rbind(Xn,Xo))/p*2
Kno = tcrossprod(Xn,Xo)/p*2
Knn = tcrossprod(Xn)/p*2
K = rbind(cbind(Ko,t(Kno)),cbind(Kno,Knn))

Kinv = solve(K)  
Ko_inv = solve(Ko)

data = data.frame(y = yo, ID = 1:no)
rownames(Ko) = colnames(Ko) = data$ID

res1 = GridLMM_posterior(y~1+(1|ID),data,relmat = list(ID = Ko))

n_h2 = nrow(res1$h2s_results)
Z = diag(1,no)
X = matrix(1,no)
XZ = cbind(X,Z)
ZtZ = diag(1,no)
Kinv = solve(Ko)
XtZtZX = crossprod(XZ)
inf_Ko_inv = cbind(0,rbind(0,Ko_inv))
h2_index = sapply(res1$h2s_solutions,function(x) x$h2s)
# i = which(h2_index == res1$h2s_results$ID[which.max(res1$h2s_results$posterior)])
x = which.max(res1$h2s_results$posterior)
moments = lapply(1:n_h2,function(x) {
  i = which(h2_index == res1$h2s_results$ID[x])
  mu = res1$h2s_solutions[[i]]$mu_star[1]
  h2 = res1$h2s_solutions[[i]]$h2s[1]
  if(h2 == 0) return(matrix(0,nr = nn,nc=2))
  # Co_inv = crossprod(cbind(X,Z))/(1-h2)
  # Co_inv[-1,-1] = Co_inv[-1,-1] + Ko_inv/h2
  Co_inv = XtZtZX/(1-h2) + inf_Ko_inv/h2
  Co = solve(Co_inv)
  # Co = solve(Ko_inv/h2 + diag(1/(1-h2),no))
  a_star = res1$h2s_solutions[[i]]$a_star
  b_star = res1$h2s_solutions[[i]]$b_star[1]
  # t((yo-mu)) %*% solve(h2*Ko + diag(1-h2,no)) %*% (yo-mu)/2
  res = sapply(1:nn,function(j) {
    # KjbtKjdinv = Kno[j,,drop=F] %*% Ko_inv
    # Kjinv_A = 1/(Knn[j,j] - KjbtKjdinv %*% t(Kno[j,,drop=F]))/(h2)
    # Kjinv_B = - t(KjbtKjdinv) %*% Kjinv_A
    KjbtKjdinv = c(0,Kno[j,,drop=F]) %*% inf_Ko_inv
    Kjinv_A = 1/(Knn[j,j] - KjbtKjdinv %*% c(0,Kno[j,,drop=F]))/(h2)
    Kjinv_B = - t(KjbtKjdinv) %*% Kjinv_A
    # Kjinv_D = Kinv_o +  t(KjbtKjdinv) %*% Kjinv_A %*% KjbtKjdinv
    # Kjinv_D_inv = Ko - crossprod(Kno[j,,drop=F])/Knn[j,j]
    C_a = 1/Kjinv_A + 1/Kjinv_A * t(Kjinv_B) %*% Co %*% Kjinv_B * 1/Kjinv_A
    C_b = -Co %*% Kjinv_B %*% (1/Kjinv_A)
    
    # u_star_hat = t(C_b) %*% (yo-mu)/(1-h2)
    u_star_hat = t(C_b) %*% crossprod(XZ,yo/(1-h2))
    var = b_star/a_star*C_a * 2*a_star/(2*a_star-2)
    # c(u_star_hat,var)
    c(u_star_hat,C_a)
  })
  t(res)
  # Cinv = solve(K)/h2
  # Cinv[c(1:no),c(1:no)] = Cinv[c(1:no),c(1:no)] + ZtZ/(1-h2)
  # C = solve(Cinv)
  C = Co
  # u_star = C %*% Z %*% (yo-mu)/(1-h2)
  # var = b_star/a_star*diag(Co) * 2*a_star/(2*a_star-2)
  a = C %*% crossprod(cbind(X,Z),yo/(1-h2))
  mu_star = a[1]
  u_star = a[-1]
  var = b_star/a_star*diag(Co)[-1] * 2*a_star/(2*a_star-2)
  return(cbind(u_star,var))
  
  # u_hat = Co %*% (yo-mu)/(1-h2)
  
    
  #   # 
  #   # C_A = solve(Kjinv_A - t(Kjinv_B) %*% Kjinv_D_inv %*% Kjinv_B)
  #   # C_Bt = -C_A %*% t(Kjinv_B) %*% Kjinv_D_inv
  #   # C_D = solve(Kjinv_D + Kjinv_B %*% Kjinv_A %*% t(Kjinv_B))
  #   # 
  #   Kj = rbind(cbind(Knn[j,j],Kno[j,,drop=F]),cbind(t(Kno[j,,drop=F]),Ko))
  #   Cinv = solve(Kj)/h2
  #   Cinv[-1,-1] = Cinv[-1,-1] + diag(1/(1-h2),no)
  #   C = solve(Cinv)
  #   # Kinv_o = Ko_inv - Ko_inv %*% t(Kno[j,,drop=F]) %*% solve(Knn[j,j] - Kno[j,,drop=F] %*% Ko_inv %*% t(Kno[j,,drop=F]), t(Ko_inv %*% t(Kno[j,,drop=F])))
  #   # Kinv_o_inv = Ko - crossprod(Kno[j,,drop=F])/Knn[j,j]
  #   # Kinv_n = 
  #   # C_nn = ()
  #   
  # })
  # Cnn = 1/()
  # C = solve(ZtZ/(1-h2) + Kinv/h2)
  # u_hat = C %*% Z %*% (data$y - mu)/(1-h2)
  # u_star_hat = Kno %*% Kinv %*% u_hat
  # u_star_hat2 = h2*Kno %*% Z %*% (data$y - mu)/(1-h2)
  # a_star = res1$h2s_solutions[[i]]$a_star
  # b_star = res1$h2s_solutions[[i]]$b_star[1]
  # lambda_inv = sapply(1:nn,function(j) {
  #   knoKinv = Kno[j,] %*% Kinv
  #   h2*(Knn[j,j] - knoKinv %*% t(Kno[j,,drop=F]) + knoKinv %*% C %*% t(knoKinv))
  # })
  # vars = a_star*lambda_inv/b_star * 2*a_star/(2*a_star-2)
  # # u_star_hat[is.na(u_star_hat)]
  # cbind(u_star_hat,vars)
})

mean = rowSums(sweep(do.call(cbind,lapply(moments,function(x) x[,1])),2,res1$h2s_results$posterior,'*'))
var = rowSums(sweep(do.call(cbind,lapply(moments,function(x) x[,1]^2+x[,2] - mean^2)),2,res1$h2s_results$posterior,'*'))


res2 = mixed.solve(y= c(yo,rep(NA,nn)), K = K,SE=T) #/sd(yo)
h2 = res2$Vu/(res2$Vu+res2$Ve)
mu = res2$beta[1]
res2b = mixed.solve(y=yo,K=Ko,SE=T)
h2 = res2b$Vu/(res2b$Vu+res2b$Ve)
mu = res2b$beta[1]
res2c = mixed.solve(y= c(yo2/sd(yo2),rep(NA,nn)), K = K,SE=T)

Vinv = solve(res2b$Vu*Ko+res2b$Ve*diag(1,no))
u_star2 = res2b$Vu*Ko %*% Vinv %*% (yo-res2b$beta[1])


library(brms)
# rownames(K) = colnames(K) = 1:nrow(K)
# res3 = brm(y ~ (1|ID),cov_ranef = list(ID = K),data = data)
sKo = svd(Ko)
data$qty = t(sKo$u) %*% data$y
data$qt1 = t(sKo$u) %*% matrix(1,nrow(data))
D = diag(sKo$d);rownames(D) = colnames(D) = data$ID
res3b = brm(qty~0+qt1+(1|ID),cov_ranef = list(ID = D),data=data,save_all_pars=T,chains=1)
res3b_s = as.matrix(res3b)
res3b_s = res3b_s[,grep('r_ID',colnames(res3b_s))]
res3b_s2 = res3b_s %*% t(sKo$u)
plot(mean,colMeans(res3b_s2));abline(0,1)
plot(mean,res2b$u);abline(0,1)
plot(res2b$u,colMeans(res3b_s2));abline(0,1)
plot(u_star,colMeans(res3b_s2));abline(0,1)
plot(var,apply(res3b_s2,2,var));abline(0,1)
plot(res2b$u.SE^2,apply(res3b_s2,2,var));abline(0,1)

res3c = mixed.solve(data$qty,K=D,X=matrix(data$qt1))
res3d = stan_lmer(qty~0+qt1+(1|ID),data=data,prior = normal(location=0,scale = sKo$d))


res2 = mixed.solve(y = c(yo,rep(NA,nn)),Z = rbind(Xo,Xn))
h2 = res2$Vu/(res2$Vu+res2$Ve)
mu = res2$beta[1]
plot(Xo %*% res2$u,u_hat);abline(0,1)
plot(Xn %*% res2$u,u_star_hat);abline(0,1)
plot(Xn %*% res2$u,mean);abline(0,1)