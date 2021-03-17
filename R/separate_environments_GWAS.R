library(GridLMM)
library(data.table)
library(foreach)
library(doParallel)
library(tidyverse)

accessions = fread('../SeeDs/MexicoAccession.txt',data.table=F,h=F)
# 
# # make LOO K's
# # drop markers with >10% missing data and MAF < 0.05
# # impute missing values with mean
prep_genos = function(mat) {
  individuals = tibble(Name = mat[,1],ID = sapply(mat[,1],function(x) strsplit(x,'.',fixed=T)[[1]][1]))
  individuals = individuals[!duplicated(individuals$ID),]
  individuals = individuals[match(accessions$V1,individuals$ID),]
  mat = mat[match(individuals$Name,mat[,1]),]
  rownames(mat) = mat[,1]
  mat = mat[,-1]
  mat = as.matrix(mat)

  mat = mat[,apply(mat,2,function(x) mean(!is.na(x)) >= 0.1 & (.5-abs(.5-mean(x,na.rm=T)))>0.05)]
  means = colMeans(mat,na.rm=T)
  mat[is.na(mat)] = matrix(means,nrow = nrow(mat),ncol = ncol(mat),byrow=T)[is.na(mat)]
  mat = 2*mat
  mat = sweep(mat,2,colMeans(mat),'-')
  mat
}
# Ks = list()
# for(chr in 1:10) {
#   print(chr)
#   mat<-fread(cmd = paste('unzip -p ./Genotypes/Ch',chr,'Merged.hmp.txt.zip',sep=""),skip=1,sep='\t',header=TRUE,data.table=F)
#   mat = prep_genos(mat)
#   Ks[[chr]] = list(K = tcrossprod(mat)/ncol(mat),m = ncol(mat))
# }
# 
# # names = do.call(cbind,lapply(Ks,function(x) rownames(x$K)))
# LOO_ks = list()
# try(dir.create('LOO_Ks'))
# for(chr in 1:10) {
#   K = 0
#   m = 0
#   for(chr2 in 1:10) {
#     if(chr == chr2) next
#     Ki = Ks[[chr2]]
#     K = K + Ki$K*Ki$m
#     m = m + Ki$m
#   }
#   K = K/m
#   LOO_ks[[chr]] = K
#   write.csv(K,file = sprintf('LOO_Ks/K_chr_%02d.csv',chr))
# }

# #load in the high quality unimputed genos for K matrix
# load('./Genotypes/HighFilteredReduced.Rimage')
# genoMDS<-as.matrix(genoMDS)
# genoMDS[which(is.na(genoMDS))]<-0.5
# genoMDS = genoMDS[match(accessions$V1,rownames(genoMDS)),]
# K_centered=sweep(genoMDS,2,colMeans(genoMDS,na.rm = TRUE),'-')
# K = tcrossprod(K_centered) / ncol(K_centered)
# rownames(K) = colnames(K) = rownames(genoMDS)
# K = K/mean(diag(K))
# dim(K)



phenotypes<-c('Days to anthesis','Plant height','Anthesis silking interval','Field weight','Bare cob weight','Grain weight per hectare')
#environments<-c("wet6MonthGrowAvg","pre6MonthGrowAvg","cld6MonthGrowAvg","Isothermality")
environments<-c('altitude','meanTemp','annualPrecipitation')

phenot = phenotypes[1]
env = environments[1]

chr = as.numeric(commandArgs(t=T)[1])
if(is.na(chr)) chr = 10
# chr = 10
results = c()
try(dir.create('Separate_GWAS_results'))
# for(chr in 1:10) {
  mat<-fread(cmd = paste('unzip -p ./Genotypes/Ch',chr,'Merged.hmp.txt.zip',sep=""),skip=1,sep='\t',header=TRUE,data.table=F)
  
  mat = prep_genos(mat)
  
  K = fread(sprintf('LOO_Ks/K_chr_%02d.csv',chr),data.table=F)
  rownames(K) = colnames(K) = sapply(K[,1],function(x) strsplit(x,'.',fixed=T)[[1]][1])
  K = as.matrix(K[,-1])
  rownames(mat) = rownames(K)
  
  results = c()
for(phenot in phenotypes) {
  for(env in environments) {
    try({
      load(paste('./GWASResiduals/',phenot,'_',env,'_GWASResiduals.Rimage',sep=""))
      data = jimbo
      data$env = data$trialElv
      
      data$observation_ID = paste(data$SampleID,data$tester,data$year,data$locale,sep=':')
      data$trial = paste(data$locale,data$year,sep=':')
      data = subset(data,SampleID %in% rownames(K))
      
      trial. = unique(data$trial)[1]
      for(trial. in unique(data$trial)) {
        data_trial = subset(data,trial == trial.)
        
        res = GridLMM_GWAS(value ~ tester + (1|SampleID),data = data_trial,
                           test_formula=~1,reduced_formula=~0,
                           X_ID = 'SampleID',X = mat,
                           relmat = list(SampleID = K),method='REML',verbose=T)
        res = res$results
        res$beta = res[,max(which(grepl('beta',colnames(res))))]
        res$SE = abs(res$beta/sqrt(res$F.1))
        res$MAF = .5-abs(.5-colMeans(sweep(mat,2,apply(mat,2,min),'-'))/2)
        results = bind_rows(results,data.frame(phenotype = phenot,environment = env,chr=chr,trail = trial.,env = data_trial$env[1],
                                               X_ID = res$X_ID,beta = res$beta,SE = res$SE,MAF = res$MAF))
      }
    }
  }
}
fwrite(results,file = sprintf('Separate_GWAS_results/Chr_%02d_results.csv',chr))