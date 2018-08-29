library(GridLMM)

set.seed(1)
# simulation with 1 RE, n=100
n = 4000
nG = 20
data = data.frame(Group = factor(rep(1:nG,each = n/nG)),cov = sample(c(-1,1),n,replace=T)) # cov is a covariate for all tests


h2 = 0.15
data$y = sqrt(h2)*rnorm(nG)[data$Group] + sqrt(1-h2)*rnorm(n)

library(lme4)
m1 = lmer(y~(1|Group),data)
vars = as.data.frame(VarCorr(m1))$vcov;vars/sum(vars)
x = model.matrix(~1,data)
# x = matrix(rnorm(n*20),n)
g = model.matrix(~0+Group,data)
chol_Vi_R = diag(1,n)
r = GridLMM_test_setTest(as.matrix(data$y),chol_Vi_R,g*sqrt(ncol(g)),x)
r$taus
r2 = GridLMM_test_setTest2(as.matrix(data$y),chol_Vi_R,g*sqrt(ncol(g)),x)
r2$taus

microbenchmark(GridLMM_test_setTest(as.matrix(data$y),chol_Vi_R,g*sqrt(ncol(g)),x),
               GridLMM_test_setTest2(as.matrix(data$y),chol_Vi_R,g*sqrt(ncol(g)),x),
               times = 10)

library(microbenchmark)
microbenchmark(lmer(y~(1|Group),data),GridLMM_test_setTest(as.matrix(data$y),chol_Vi_R,g*sqrt(ncol(g)),x),times=6)
tau = 0.611854
G = tcrossprod(g*sqrt(ncol(g)))
determinant(tau/ncol(g)*G + (1-tau)*diag(1,n))
R = chol(tau/ncol(g)*G + (1-tau)*diag(1,n))
ystd = solve(t(R),data$y)
xstd = solve(t(R),x)
e = resid(lm(ystd~0+xstd))
sum(e^2)


p = 1000
Z = model.matrix(~0+Group,data)
K = Z %*% t(Z)

Xg = matrix(sample(c(0,1),nG*p,replace=T),nrow = nG)
# Xg[,4:6] = Xg[,1:3]
rownames(Xg) = 1:nG
X = Z %*% Xg
X = matrix(sample(c(0,2),n*p,replace=T),nrow = n)
X = sweep(X,2,colMeans(X),'-')
data$ID = factor(1:nrow(data))
rownames(X) = data$ID
K = X %*% t(X) / p
rownames(K) = colnames(K) = data$ID

X1 = X
X = matrix(sample(c(0,1),n*p,replace=T),nrow = n)
X = sweep(X,2,colMeans(X),'-')
rownames(X) = data$ID

h2 = 0.5
prop_X = 0.1
cV = chol(h2*K + (1-h2)*diag(1,n)) # This K doesn't not include the testing X
K = K/2 + 1/2*X %*% t(X) / p # This K does include the testing X

beta = 1*c(1,0,0,rep(0,p-3))  # first SNP is real
beta_cov = 0
beta_gxe = c(1,0,1,rep(0,p-3))

g = X %*% beta + data$cov*(X %*% beta_gxe);g=g/sd(g)*1
e = t(cV) %*% rnorm(n);e=e/sd(e)
y = data$cov * beta_cov + sqrt(prop_X) * c(g) + sqrt(1-prop_X) * e
y = y/sd(y)
data$y = y



g0 = GridLMM_ML(y~cov+(1|Group),data,relmat = list(Group = tcrossprod(scale_SNPs(Xg,scaleX = F))/p))
g0$results
g0 = GridLMM_ML(y~cov+(1|ID),data,relmat = list(ID = K))
g0$results

g1 = GridLMM_GWAS(y~cov + (1|Group) + (0+cov|Group),~1,~0,data = data,X = Xg,X_ID = 'Group',
                  # h2_start = c(0.9,0.1),
                  algorithm = 'Full',
                  h2_step = 0.1, save_V_folder = 'V_folder',mc.cores=1)
proximal_markers = lapply(1:p,function(x) seq(max(1,x-10),min(p,x+10)))
proximal_markers[-c(1:2)] = lapply(proximal_markers[-c(1:2)],function(x) x[x>2])
g2 = GridLMM_GWAS(y~cov + (1|Group),~1,~0,data = data,X = Xg,X_ID = 'Group',h2_step = 0.1,mc.cores=1,algorithm = 'Full',
                  proximal_markers = proximal_markers)
g2b = GridLMM_GWAS(y~cov + (1|Group),~1,~0,data = data,X = Xg,X_ID = 'Group',proximal_markers = proximal_markers, 
                   V_setup = g2$setup$V_setup)
head(g1$results)
head(g2$results)
head(g2b$results)

g3 = GridLMM_GWAS(y~cov + (1|ID),~1,~0,data = data,X = X,X_ID = 'ID',h2_step = 0.01,mc.cores=1,algorithm = 'Fast',
                  method='ML',centerX=T,relmat = list(ID = K))
g4 = GridLMM_GWAS(y~cov + (1|ID),~1,~0,data = data,X = X,X_ID = 'ID',h2_step = 0.01,mc.cores=1,algorithm = 'Fast',
                  centerX = T,relmat = list(ID = list(K=K,p=2*p)),h2_start = c(ID = 0.43),
                  proximal_markers = proximal_markers,diagonalize = F,method='ML')
# qq_p_plot(g4$results$p_value_REML)
qq_p_plot(list(g3=g3$results$p_value_ML,g4=g4$results$p_value_ML))
head(g3$results)
head(g4$results)

plot(-log10(g4$results$p_value_ML))
plot(-log10(g1$results$p_value_RML))


# Run GxEMMAnet
library(glmnet)

h2 = 0.5
prop_X = 0.6

p = 1000
p_g = 10 # include 10 coefficients
p_gxe = 0 # include 0 gxe coefficients
beta = c(rep(1,p_g),rep(0,p-p_g))
beta_gxe = c(rep(0,p_g),rep(1,p_gxe),rep(0,p-p_g-p_gxe))  

X = cbind(1,matrix(sample(c(0,1),n*(1000-1),replace=T),nrow = n))
rownames(X) = data$ID
g = X %*% beta + data$cov*(X %*% beta_gxe);g=g/sd(g)
e = rnorm(n);e=e/sd(e)
y = data$cov * beta_cov + sqrt(prop_X) * c(g) + sqrt(1-prop_X) * e + sqrt(h2/(1-h2))* Z %*% rnorm(nG) 
y = y/sd(y)
data$y = y

# LASSO with random effect
i = sample(1:nrow(X))
gLASSO = GridLMMnet(y~cov + (1|Group),data,X,X_ID = 'ID',alpha = 1,h2_divisions = 100,diagonalize = T)

# a = GridLMM_ML(y~cov + (1|Group),data,X=X[,which(gLASSO$beta[-c(1:2),80] != 0)],ML=T,REML=T,initial_step = 0.01);a$results

# LASSO without random effect for comparison
gLASSO_0 = glmnet(cbind(1,data$cov,X),y,alpha = 1,penalty.factor = c(0,0,rep(1,p)),standardize = FALSE)
# gLASSO_0$beta[1,] = gLASSO_0$a0
# gLASSO_0$a0 = 0

par(mfrow=c(2,1))
cols = c(1,1,rep(2,p_g),rep(3,p_gxe),rep(3,p-p_g-p_gxe))
plot(gLASSO,col=cols[-1],'lambda')
plot(gLASSO_0,col=cols[-1],'lambda')

gLcv = cv.glmnet(cbind(1,data$cov,X),y,alpha = 1,penalty.factor = c(0,0,rep(1,p)),standardize = FALSE,keep = T)
gLASSOcv = GridLMMnet(y~cov + (1|Group),data,X,alpha = 1,h2_divisions = 8*4,diagonalize = F,foldid = gLcv$foldid,mc.cores = 8, save_V_folder = 'V_folder')
plot(gLASSOcv);points(log(gLcv$lambda),gLcv$cvm)
plot(gLcv)

gLASSO_setup = GridLMMnet_setup(y~cov + (1|Group),data,X,alpha = 1,diagonalize = F,foldid = gLcv$foldid,mc.cores = 8, save_V_folder = 'V_folder')
gLASSO_setup$h2s_matrix = setup_Grid(names(gLASSO_setup$V_setup$RE_setup),h2_step = 1/(8*4))
cl = start_cluster(8,'FORK')
gLASSO_setup$V_setup = calculate_Grid(gLASSO_setup$V_setup,gLASSO_setup$h2s_matrix)
stop_cluster(cl)

lambda = get_lambda_sequence(gLASSO_setup,alpha = 1)
cl = start_cluster(8,'FORK')
results = run_GridLMMnet(gLASSO_setup,alpha = 1,lambda=lambda)
gLASSOcv2 = collect_results_GridLMMnet(gLASSO_setup,results,lambda)
stop_cluster(cl)
plot(gLASSOcv)
plot(gLASSOcv2)

library(sommer)
library(lme4)

# simulation with 1 RE, n=100
n = 200
nG = 50
data = data.frame(Group = factor(rep(1:nG,each = n/nG)),cov = sample(c(-1,1),n,replace=T)) # cov is a covariate for all tests
p = 3
Z = model.matrix(~0+Group,data)
K = Z %*% t(Z)
Z2 = data$cov * Z

X = matrix(sample(c(0,1),n*p,replace=T),nrow = n)
data$ID = 1:nrow(data)
data = data.frame(data,X)

h2 = 0.7
prop_X = 0.

beta = c(1,1,0)  # first SNP is real

g = X %*% beta;g=g/sd(g)
ge = sqrt(h2)* Z %*% rnorm(nG) 
e = sqrt(1-h2)*rnorm(n)#,0,1+(ge-min(ge))/sd(ge))
y = sqrt(prop_X) * c(g) + sqrt(1-prop_X) * (ge+e)
y = y/sd(y)
data$y = y
# data$y[order(ge)[1:50]] = min(data$y) - 6
quantNorm = function(x){qnorm(rank(x,ties.method = "average")/(length(x)+1))}
# data$y = quantNorm(data$y)

lme1 = lmer(y~X1+X2+X3 + (1|Group),data,REML=F)
vcov = as.data.frame(VarCorr(lme1))$vcov
vcov[1]/sum(vcov)
m1 = mmer2(y~X1+X2+X3,random = ~Group,data = data)
pin(m1,h2~V1/(V1+V2))
gxe1 = GridLMM_posterior(y~X1+X2+X3 + (1|Group),data,thresh_nonzero = 50,mc.cores=1)
g2 = GridLMM_ML(y~X1+X2+X3 + (1|Group),data,tolerance = 1e-5)
g2$results
plot(gxe1$h2s_results$posterior~gxe1$h2s_results$Group,type='l');abline(v=vcov[1]/sum(vcov));abline(v=sum(gxe1$h2s_results$posterior*gxe1$h2s_results$Group),col=2);abline(v=h2,col=3)
abline(v=c(g2$results$Group.ML,g2$results$Group.REML),col=3:4)
sum(gxe1$h2s_results$posterior*gxe1$h2s_results$Group)
