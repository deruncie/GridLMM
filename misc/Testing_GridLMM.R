library(GridLMM)

set.seed(1)
# simulation with 1 RE, n=100
n = 200
nG = 50
data = data.frame(Group = factor(rep(1:nG,each = n/nG)),cov = sample(c(-1,1),n,replace=T)) # cov is a covariate for all tests
p = 1000
Z = model.matrix(~0+Group,data)
K = Z %*% t(Z)

Xg = matrix(sample(c(0,1),nG*p,replace=T),nrow = nG)
Xg[,4:6] = Xg[,1:3]
rownames(Xg) = 1:nG
X = Z %*% Xg
# X = matrix(sample(c(0,1),n*p,replace=T),nrow = n)
X = sweep(X,2,colMeans(X),'-')
data$ID = factor(1:nrow(data))
rownames(X) = data$ID
# K = X %*% t(X) / p
rownames(K) = colnames(K) = data$ID

h2 = 0.5
prop_X = 0.5

beta = c(1,1,0,rep(0,p-3))  # first SNP is real
beta_cov = 0
beta_gxe = c(1,0,1,rep(0,p-3))

g = X %*% beta + data$cov*(X %*% beta_gxe);g=g/sd(g)
e = sqrt(h2)* Z %*% rnorm(nG) + sqrt(1-h2)*rnorm(n);e=e/sd(e)
y = data$cov * beta_cov + sqrt(prop_X) * c(g) + sqrt(1-prop_X) * e
y = y/sd(y)
data$y = y





# Run GxEMMA_GWAS

prop_X= 0.5
beta = c(1,1,0,rep(0,p-3))  # first SNP is real
beta_cov = 0
beta_gxe = 0*c(1,0,1,rep(0,p-3))

g = X %*% beta + data$cov*(X %*% beta_gxe);g=g/sd(g)
e = sqrt(h2)* Z %*% rnorm(nG) + sqrt(1-h2)*rnorm(n);e=e/sd(e)
y = data$cov * beta_cov + sqrt(prop_X) * c(g) + sqrt(1-prop_X) * e
y = y/sd(y)
data$y = y

g0 = GridLMM_ML(y~cov+(1|Group)+(0+cov|Group),data,relmat = list(Group = tcrossprod(scale_SNPs(Xg[,1:50],scaleX = F))/50))
g0$results
g1 = GridLMM_GWAS(y~cov + (1|Group) + (0+cov|Group),~1,~0,data = data,X = Xg[,1:50],X_ID = 'Group',h2_divisions = 10, save_V_list = 'V_folder',mc.cores=1)
proximal_matrix = diag(1,50)
diag(proximal_matrix[-1,]) = 1
diag(proximal_matrix[-c(1:2),]) = 1
diag(proximal_matrix[-c(1:3),]) = 1
proximal_matrix[1:4,1:4]=1
g2 = GridLMM_GWAS(y~cov + (1|Group) + (0+cov|Group),~1,~0,data = data,X = Xg[,1:50],X_ID = 'Group',h2_divisions = 10,mc.cores=1,proximal_matrix = proximal_matrix)
g3 = GridLMM_GWAS(y~cov + (1|Group) + (0+cov|Group),~1,~0,data = data,X = Xg[,1:50],X_ID = 'Group',proximal_matrix = proximal_matrix,h2_divisions = 10)
g3b = GridLMM_GWAS(y~cov + (1|Group) + (0+cov|Group),~1,~0,data = data,X = Xg[,1:50],X_ID = 'Group',h2_divisions = 10,
                   RE_setup = g3$setup$RE_setup, V_list = g3$setup$V_list, downdate_Xs = g3$setup$downdate_Xs)
head(g1$results)
head(g2$results)
head(g3$results)
head(g3b$results)

g4 = GridLMM_GWAS_fast(y~cov + (1|Group) + (0+cov|Group),~1,~0,data = data,X = Xg[,1:50],X_ID = 'Group',max_step = 100,h2_step = 0.01,mc.cores=1)
g4b = GridLMM_GWAS_fast(y~cov + (1|Group) + (0+cov|Group),~1,~0,data = data,X = Xg[,1:50],X_ID = 'Group',max_step = 100,h2_step = 0.01,mc.cores=1,proximal_matrix = proximal_matrix)
head(g4$results)
head(g4b$results)


# Run GxEMMAnet
library(glmnet)

h2 = 0.7
prop_X = 0.6

p_g = 10 # include 10 coefficients
p_gxe = 0 # include 0 gxe coefficients
beta = c(rep(1,p_g),rep(0,p-p_g))
beta_gxe = c(rep(0,p_g),rep(1,p_gxe),rep(0,p-p_g-p_gxe))  

# X = cbind(1,matrix(sample(c(0,1),n*(1000-1),replace=T),nrow = n))
g = X %*% beta + data$cov*(X %*% beta_gxe);g=g/sd(g)
e = rnorm(n);e=e/sd(e)
y = data$cov * beta_cov + sqrt(prop_X) * c(g) + sqrt(1-prop_X) * e + sqrt(h2/(1-h2))* Z %*% rnorm(nG) 
y = y/sd(y)
data$y = y

# LASSO with random effect
i = sample(1:nrow(X))
gLASSO = GridLMMnet(y~cov + (1|Group),data,X,alpha = 1,h2_divisions = 100,diagonalize = T)

# a = GridLMM_ML(y~cov + (1|Group),data,X=X[,which(gLASSO$beta[-c(1:2),80] != 0)],ML=T,REML=T,initial_step = 0.01);a$results

# LASSO without random effect for comparison
gLASSO_0 = glmnet(cbind(1,data$cov,X),y,alpha = 1,penalty.factor = c(0,0,rep(1,p)),standardize = FALSE)
# gLASSO_0$beta[1,] = gLASSO_0$a0
# gLASSO_0$a0 = 0

par(mfrow=c(2,1))
cols = c(1,1,rep(2,p_g),rep(3,p_gxe),rep(3,p-p_g-p_gxe))
plot(gLASSO,col=cols[-1],'lambda')
plot(gLASSO_0,col=cols[-1],'lambda')

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
