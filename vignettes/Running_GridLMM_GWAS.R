## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(GridLMM)

## ------------------------------------------------------------------------
n_geno = 300
n_chr = 5
n_markers_per_chr = 1000

# layout of Plots
Plots = matrix(1:n_geno,nc = 10) # rectangular array of plots
Plot_Row = row(Plots)
Plot_Col = col(Plots)

# create genotypes and randomize in field
data = data.frame(Geno = paste0('Geno',1:n_geno), Plot = sample(Plots))
data$Row = Plot_Row[data$Plot]
data$Col = Plot_Col[data$Plot]

# create marker map and random genotypes
map = expand.grid(pos = 1:n_markers_per_chr,Chr = 1:n_chr)[,2:1]
map = data.frame(snp = paste0('snp_',1:nrow(map)),map)
X = matrix(sample(c(0,2),n_geno*nrow(map),replace=T),ncol = nrow(map),dimnames = list(data$Geno,map$snp))

# simulate 1 effect per chromosome
beta = rep(0,nrow(map))
beta[(1:n_chr-1)*n_markers_per_chr + sample(1:n_markers_per_chr,n_chr)] = 1

# add background variation: modeled as small effects across all the markers
beta_background = rnorm(length(beta),0,0.1)

# create simulated data
background_h2 = 0.6 # relative contribution of the background
data$micro_env = rnorm(nrow(data)); data$micro_env = data$micro_env / sd(data$micro_env)
data$background = X %*% beta_background; data$background = data$background / sd(data$background)
data$SNPs = X %*% beta;data$SNPs = data$SNPs/sd(data$SNPs)
data$y = data$SNPs + sqrt(background_h2) * data$background + sqrt(1-background_h2) * data$micro_env
data$y = data$y/sd(data$y)

## ------------------------------------------------------------------------
X_centered = sweep(X,2,colMeans(X),'-') # center marker genotypes
K = tcrossprod(X_centered) / ncol(X_centered)
rownames(K) = colnames(K) = rownames(X)

## ------------------------------------------------------------------------
library(KRLS)
field = data[,c('Row','Col')]
dists = as.matrix(dist(field))
h = median(dists)
K_plot = gausskernel(field,h^2/2); diag(K_plot)=1 # 
rownames(K_plot) = colnames(K_plot) = data$Plot

## ------------------------------------------------------------------------
h2_spatial = 0.9
s2_spatial = 1
chol_K_plot = chol(h2_spatial*K_plot + diag(1-h2_spatial,nrow(data)))
data$spatial_error = t(chol_K_plot) %*% rnorm(nrow(data),0,sqrt(s2_spatial))
field_error = 0*Plots
field_error[data$Plot] = data$spatial_error
image(field_error)
data$y = data$y + data$spatial_error #

## ------------------------------------------------------------------------
gwas = GridLMM_GWAS_fast(
                        formula = y~1 + (1|Geno) + (1|Plot), # the same error model is used for each marker. It is specified similarly to lmer
                        test_formula = ~1, # this is the model for each marker. ~1 means an intercept for the marker. ~1 + cov means an intercept plus a slope on `cov` for each marker
                        reduced_formula = ~0, # This is the null model for each test. ~0 means just the error model. ~1 means the null includes the intercept of the marker, but not anything additional
                        data = data, # The dataframe to look for terms from the 3 models
                        weights = NULL, # optional observation-specific weights
                        X = X, # The matrix of markers. Note: This can be of dimension n_g x p, where n_g is the number of unique genotypes.
                        X_ID = 'Geno', # The column of the data that identifies each genotype. Each level of data$Geno should match a rowname of X
                        h2_start = NULL, # An optional vector of h2s to use as starting values for each random effect in the model. If NULL, will be calculated from the error model using relmatLmer
                        h2_step = 0.01, # step size per random effect for testing alternate values of h2
                        max_steps = 100, # maximum number of steps of size h2_step to take from h2_start
                        X_map = map, # Optional. The marker positions.
                        relmat = list(Plot = K_plot), # A list of Kernel matrices for the random effects. If X_ID (here Geno) is not included in this list, then it is calculated as tcrossprod(Xc)/ncol(Xc) where Xc is the centered (and optionally scaled) X. If any random effects are described in `error_model` but not provided here, the Kernel is assumed to be the identity matrix
                        centerX = TRUE, # Should the markers be centered when calculating the GRM (only will be done if needed for calculating the GRM),
                        scaleX = FALSE, # Should the markers be scaled to have constant variance when calculating the GRM?
                        fillNAX = FALSE, # Should missing marker data be filled in with the mean allele frequency?
                        method = 'REML', # REML = Wald test, ML = LRT, BF = calculate Bayes factors
                        mc.cores = 1, # How many cores should be used for parallel processing. Unless X is large, tends to actually be faster with mc.cores = 1
                        verbose = FALSE # Should progress be printed to the screen?
)

## ------------------------------------------------------------------------
names(gwas)

## ------------------------------------------------------------------------
gwas2 = run_GWAS_fast(X,gwas$setup,h2_step = 0.01,max_steps = 100,method = 'ML',verbose = FALSE)

## ------------------------------------------------------------------------
head(gwas$results)

## ----message=FALSE-------------------------------------------------------
library(qqman)
# manhattan plot with true SNPs highlighted
manhattan(gwas$results,'Chr','pos','p_value_REML','snp',highlight = colnames(X)[beta != 0])

# qqplot 
qq(gwas$results$p_value_REML)

## ------------------------------------------------------------------------
gwas1 = GridLMM_GWAS_fast(
                        formula = y~1 + (1|Geno), # the same error model is used for each marker. It is specified similarly to lmer
                        test_formula = ~1, # this is the model for each marker. ~1 means an intercept for the marker. ~1 + cov means an intercept plus a slope on `cov` for each marker
                        reduced_formula = ~0, # This is the null model for each test. ~0 means just the error model. ~1 means the null includes the intercept of the marker, but not anything additional
                        data = data, # The dataframe to look for terms from the 3 models
                        weights = NULL, # optional observation-specific weights
                        X = X, # The matrix of markers. Note: This can be of dimension n_g x p, where n_g is the number of unique genotypes.
                        X_ID = 'Geno', # The column of the data that identifies each genotype. Each level of data$Geno should match a rowname of X
                        h2_start = NULL, # An optional vector of h2s to use as starting values for each random effect in the model. If NULL, will be calculated from the error model using relmatLmer
                        h2_step = 0.01, # step size per random effect for testing alternate values of h2
                        max_steps = 100, # maximum number of steps of size h2_step to take from h2_start
                        X_map = map, # Optional. The marker positions.
                        relmat = list(Plot = K_plot), # A list of Kernel matrices for the random effects. If X_ID (here Geno) is not included in this list, then it is calculated as tcrossprod(Xc)/ncol(Xc) where Xc is the centered (and optionally scaled) X. If any random effects are described in `error_model` but not provided here, the Kernel is assumed to be the identity matrix
                        centerX = TRUE, # Should the markers be centered when calculating the GRM (only will be done if needed for calculating the GRM),
                        scaleX = FALSE, # Should the markers be scaled to have constant variance when calculating the GRM?
                        fillNAX = FALSE, # Should missing marker data be filled in with the mean allele frequency?
                        method = 'REML', # Should the best model be selected by REML (if False, will be selected by ML)
                        mc.cores = 1, # How many cores should be used for parallel processing. Unless X is large, tends to actually be faster with mc.cores = 1
                        verbose = FALSE # Should progress be printed to the screen?
)

## ------------------------------------------------------------------------
plot(-log10(gwas$results$p_value_REML),-log10(gwas1$results$p_value_REML),col = (beta != 0)+1);abline(0,1)

u = sort(-log10(runif(nrow(gwas$results))))
plot(sort(-log10(gwas$results$p_value_REML))~u,pch=21,cex=.5);abline(0,1)
points(sort(-log10(gwas1$results$p_value_REML))~u,pch=21,cex=.5,bg=2,col=NULL)

## ------------------------------------------------------------------------
# add a random environmental vector to data
data$env = rnorm(nrow(data))

gxe_gwas = GridLMM_GWAS_fast(
                        formula = y~1 + env + (1 + env|Geno), # the same error model is used for each marker. It is specified similarly to lmer
                        test_formula = ~1 + env, # this is the model for each marker. ~1 means an intercept for the marker. ~1 + cov means an intercept plus a slope on `cov` for each marker
                        reduced_formula = ~1, # This is the null model for each test. ~0 means just the error model. ~1 means the null includes the intercept of the marker, but not anything additional
                        data = data, # The dataframe to look for terms from the 3 models
                        weights = NULL, # optional observation-specific weights
                        X = X, # The matrix of markers. Note: This can be of dimension n_g x p, where n_g is the number of unique genotypes.
                        X_ID = 'Geno', # The column of the data that identifies each genotype. Each level of data$Geno should match a rowname of X
                        h2_start = NULL, # An optional vector of h2s to use as starting values for each random effect in the model. If NULL, will be calculated from the error model using relmatLmer
                        h2_step = 0.01, # step size per random effect for testing alternate values of h2
                        max_steps = 100, # maximum number of steps of size h2_step to take from h2_start
                        X_map = map, # Optional. The marker positions.
                        relmat = NULL, # A list of Kernel matrices for the random effects. If X_ID (here Geno) is not included in this list, then it is calculated as tcrossprod(Xc)/ncol(Xc) where Xc is the centered (and optionally scaled) X. If any random effects are described in `error_model` but not provided here, the Kernel is assumed to be the identity matrix
                        centerX = TRUE, # Should the markers be centered when calculating the GRM (only will be done if needed for calculating the GRM),
                        scaleX = FALSE, # Should the markers be scaled to have constant variance when calculating the GRM?
                        fillNAX = FALSE, # Should missing marker data be filled in with the mean allele frequency?
                        method = 'REML', # Should the best model be selected by REML (if False, will be selected by ML)
                        mc.cores = 1, # How many cores should be used for parallel processing. Unless X is large, tends to actually be faster with mc.cores = 1
                        verbose = FALSE # Should progress be printed to the screen?
)

## ------------------------------------------------------------------------
u = sort(-log10(runif(nrow(gxe_gwas$results))))
plot(sort(-log10(gxe_gwas$results$p_value_REML.1))~u,pch=21,cex=.5,main = 'SNP main effect');abline(0,1)
plot(sort(-log10(gxe_gwas$results$p_value_REML.2))~u,pch=21,cex=.5,main = 'SNP gxe effect');abline(0,1)

## ------------------------------------------------------------------------
library(snpStats)

write.plink(file.base = 'test_data',snps = as(X,'SnpMatrix'),id = data$Geno,phenotype = data$y,
            chromosome = map$Chr,position = map$pos)

## ------------------------------------------------------------------------
K_G = tcrossprod(X)/ncol(X)

## ------------------------------------------------------------------------
null_model = GridLMM_ML(formula = y~1 + (1|Geno) + (1|Plot),data = data,relmat = list(Geno = K_G, Plot = K_plot),REML = T,save_V_folder = 'V_folder',tolerance = 1e-3)

## ------------------------------------------------------------------------
null_model$results[,c('Geno.REML','Plot.REML')]

## ------------------------------------------------------------------------
V_setup = null_model$setup

## ------------------------------------------------------------------------
markers = read.delim('test_data.bim',header = F,row.names = NULL)
nrow(markers)

## ------------------------------------------------------------------------
chunks = split(1:nrow(markers),gl(10,nrow(markers)/10))

## ------------------------------------------------------------------------
library(doParallel)
library(foreach)
full_results = foreach(chunk = chunks,i = 1:length(chunks),.combine = 'rbind') %do% {
  print(sprintf('Running chunk %d',i))
  X_chunk = read.plink('test_data',select.snps = chunk)
  X_chunk = as(X_chunk$genotypes,'numeric') # convert to numeric matrix
  # Run the GWAS for the chunk
  results_chunk = GridLMM_GWAS_fast(formula = y~1 + (1|Geno) + (1|Plot),
                                    test_formula = ~1,
                                    reduced_formula = ~0,
                                    data = data,
                                    X = X_chunk,
                                    X_ID = 'Geno',
                                    relmat = list(Geno = K_G, Plot = K_plot),
                                    V_setup = V_setup,  # This tells the function to reuse existing matrices
                                    h2_start = null_model$results[,c('Geno.REML','Plot.REML')],
                                    method = 'REML',
                                    verbose = F)
  results_chunk$results # return just the results
}

## ------------------------------------------------------------------------
plot(-log10(full_results$p_value_REML),-log10(gwas$results$p_value_REML));abline(0,1)

