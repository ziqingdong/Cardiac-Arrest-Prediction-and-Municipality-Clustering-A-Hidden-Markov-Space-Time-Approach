# preliminaries
rm(list = ls())
library(coda)
library(splines)
library(LaplacesDemon)
library(pscl)
library(tidyverse)
source("MCMC_HMM.R")
# source("compile_MCMC.R")     # for Mac
source("compile_MCMC_win.R")   # for Windows
source("MLE_HMM.R")

# load data
setwd("..")
setwd("Data")
data = readRDS("ohca_data_rit.rds") 
C = as.matrix(readRDS("adjacency_matrix.rds"))
setwd("..")
setwd("R")

# remove municipalities with population 0
data = data %>%
  filter(Year < 2017)
data0 = data
data = data[data0$Pop>0,]
C = C[data0$Pop[1:117]>0,data0$Pop[1:117]>0]
datae = data[data$Year<2016,]

# organize responses and offset 
n = 115; TT = length(unique(datae$Year))+1; TTe = TT-1
municipalities = data$Municipality[1:n]
years = sort(unique(data$Year)); yearse = years[-TT]
Y = matrix(data$ohca,n,TT); dimnames(Y) = list(municipalities,years)
Ye = Y[,-TT]
M = matrix(data$Pop,n,TT); dimnames(M) = list(municipalities,years)
Me = M[,-TT]

# organize covariates affecting the latent states
Ze = array(0,c(7,n,TTe))
for(i in 1:n){
  tmp1 = datae$proportion.f[datae[,1]==municipalities[i]]
  #  tmp2 = is.na(tmp1)|is.nan(tmp1); tmp1[tmp2] = 0  # in case of missing data
  tmp3 = datae$big.city[datae[,1]==municipalities[i]]
  Ze[,i,] = rbind(tmp1,
                  # tmp2,
                  datae$proportion.2[datae[,1]==municipalities[i]],
                  datae$proportion.3[datae[,1]==municipalities[i]],
                  datae$proportion.4[datae[,1]==municipalities[i]],
                  datae$proportion.5[datae[,1]==municipalities[i]],
                  datae$proportion.6[datae[,1]==municipalities[i]],
                  datae$big.city[datae[,1]==municipalities[i]])
}
namesZ = c("proportion.f","proportion.2","proportion.3","proportion.4",
           "proportion.5","proportion.6","big.city")
dimnames(Ze) = list(namesZ,municipalities,yearse)

# organize covariates affecting directly the responses
Xe = NULL
covariates_to_move = c("proportion.f","proportion.2","proportion.3","proportion.4",
                       "proportion.5","proportion.6","big.city")  
Xe = array(0, c(length(covariates_to_move), n, TTe))
dimnames(Xe) = list(covariates_to_move, municipalities, yearse)
for (j in 1:length(covariates_to_move)) {
  covariate = covariates_to_move[j]
  cov_index = match(covariate, dimnames(Ze)[[1]])  
  Xe[j, , ] = Ze[cov_index, , ]
}
covariate_indices = match(covariates_to_move, dimnames(Ze)[[1]])
Ze = Ze[-covariate_indices, , , drop = FALSE]
namesZ = dimnames(Ze)[[1]]

# transform adjacency matrix 
cmax = max(rowSums(C))  
cN = matrix(0,n,cmax)
for(i in 1:n){
  ind = which(C[i,]==1)
  if(length(ind)>0) cN[i,1:length(ind)] = ind
}

# covariates for predictions two times ahead
years_new = c(2016,2017)
Xn = NULL
Xn = array(Xe[,,TTe], c(dim(Xe)[1], n, 2))
dimnames(Xn) = list(dimnames(Xe)[[1]], municipalities, years_new)
Zn = array(Ze[,,TTe],c(dim(Ze)[1],n,2))
dimnames(Zn) = list(namesZ,municipalities,years_new)
Mn = matrix(Me[,TTe],n,2)
dimnames(Mn) = list(municipalities,years_new)
timesn = 2016:2017


#---- naive prediction ----
# previous value
mYh0a = Y[,TT-1]
err0a = Y[,TT]-mYh0a
chi20a = sum(err0a^2/(mYh0a+0.5))
# exponential average of order 3
tmp = rowMeans(Y[,1:3])
for(t in 4:TTe) tmp = tmp/2+Y[,t]/2
mYh0b = tmp
err0b = Y[,TT]-mYh0b
chi20b = sum(err0b^2/(mYh0b+0.5))

#---- models without spatial component ----
nit = 100000 # number of iterations

# with 1 latent state
set.seed(0)
taube = 0.1 # sd for proposal on beta 
est01mle = MLE_HMM(Y=Ye,X=Xe,M=Me,k=1)
est01 = MCMC_HMM(Y=Ye,X=Xe,Z=Ze,M=Me,k=1,nit=nit,taube=taube,same=TRUE,Xn=Xn,Zn=Zn,Mn=Mn,
                        Be=est01mle$Be,timesn=timesn) # from MLE
mYh01 = est01$mYh[,1]
err01 = Y[,TT]-mYh01
chi201 = sum(err01^2/(mYh01+0.5))

# with 2 latent states
set.seed(0)
est02mle = MLE_HMM(Y=Ye,X=Xe,M=Me,k=2)
for(it in 1:5){
  print(it)
  tmp = MLE_HMM(Y=Ye,X=Xe,M=Me,k=2,rand.start = TRUE)
  if(tmp$lk>est02mle$lk) est02mle = tmp
}
taube = 0.07; tauga = 3; taude = 1.4
est02 = MCMC_HMM(Y=Ye,X=Xe,Z=Ze,M=Me,k=2,nit=nit,taube=taube,tauga=tauga,taude=taude,
                        same=TRUE,Xn=Xn,Zn=Zn,Mn=Mn,Be=est02mle$Be,U=est02mle$Ul,timesn=timesn) # from MLE
mYh02 = est02$mYh[,1]
err02 = Y[,TT]-mYh02
chi202 = sum(err02^2/(mYh02+0.5))

# with 3 latent states
set.seed(0)
est03mle = MLE_HMM(Y=Ye,X=Xe,M=Me,k=3)
for(it in 1:10){
  print(it)
  tmp = MLE_HMM(Y=Ye,X=Xe,M=Me,k=3,rand.start = TRUE)
  if(tmp$lk>est03mle$lk) est03mle = tmp
}
set.seed(0)
taube = 0.07; tauga = 1.3; taude = 2
est03 = MCMC_HMM(Y=Ye,X=Xe,Z=Ze,M=Me,k=3,nit=nit,taube=taube,tauga=tauga,taude=taude,
                        same=TRUE,Xn=Xn,Zn=Zn,Mn=Mn,Be=est03mle$Be,U=est03mle$Ul,timesn=timesn) # from MLE
mYh03 = est03$mYh[,1]
err03 = Y[,TT]-mYh03
chi203 = sum(err03^2/(mYh03+0.5))

# with 4 latent states
set.seed(0)
est04mle = MLE_HMM(Y=Ye,X=Xe,M=Me,k=4)
for(it in 1:20){
  print(it)
  tmp = MLE_HMM(Y=Ye,X=Xe,M=Me,k=4,rand.start = TRUE)
  if(tmp$lk>est04mle$lk) est04mle = tmp
}
set.seed(0)
taube = 0.07; tauga = 20; taude = 12 
est04 = MCMC_HMM(Y=Ye,X=Xe,Z=Ze,M=Me,k=4,nit=nit,taube=taube,tauga=tauga,taude=taude,
                        same=TRUE,Xn=Xn,Zn=Zn,Mn=Mn,Be=est04mle$Be,U=est04mle$Ul,timesn=timesn) # from MLE
mYh04 = est04$mYh[,1]
err04 = Y[,TT]-mYh04
chi204 = sum(err04^2/(mYh04+0.5))

#---- models with spatial component ----
# with 1 latent state (with spatial component)
est1 = est01
err1 = err01
chi21 = chi201

# with 2 latent state (with spatial component)
set.seed(0)
taube = 0.08; tauga = 2.1; taude = 1.5 # sd for proposal on beta, gamma, delta
est2 = MCMC_HMM(Y=Ye,X=Xe,Z=Ze,M=Me,cN=cN,k=2,nit=nit,taube=taube,tauga=tauga,taude=taude,
                       same=TRUE,Xn=Xn,Zn=Zn,Mn=Mn,Be=est02mle$Be,U=est02mle$Ul,
                       timesn=timesn)   
mYh2 = est2$mYh[,1]
err2 = Y[,TT]-mYh2
chi22 = sum(err2^2/(mYh2+0.5))

# with 3 latent states (with spatial component)
set.seed(0)
taube = 0.07; tauga = .9; taude = 2 # sd for proposal on beta, gamma, delta   
est3 = MCMC_HMM(Y=Ye,X=Xe,Z=Ze,M=Me,cN=cN,k=3,nit=nit,taube=taube,tauga=tauga,taude=taude,
                       same=TRUE,Xn=Xn,Zn=Zn,Mn=Mn,Be=est03mle$Be,
                       U=est03mle$Ul,timesn=timesn)  
mYh3 = est3$mYh[,1]
err3 = Y[,TT]-mYh3
chi23 = sum(err3^2/(mYh3+0.5))

# with 4 latent states
set.seed(0)
taube = 0.07; tauga = 20; taude = 10 # sd for proposal on beta, gamma, delta  
est4 = MCMC_HMM(Y=Ye,X=Xe,Z=Ze,M=Me,cN=cN,k=4,nit=nit,taube=taube,tauga=tauga,taude=taude,
                       same=TRUE,Xn=Xn,Zn=Zn,Mn=Mn,Be=est04mle$Be,U=est04mle$Ul,timesn=timesn)   # with 4 latent state (with spatial component)
mYh4 = est4$mYh[,1]
err4 = Y[,TT]-mYh4
chi24 = sum(err4^2/(mYh4+0.5))

#---- model selection ----
eldpw01 = est01$eldpwv[length(est01$eldpwv)]; seldpw01 = est01$seldpwv[length(est01$seldpwv)]
eldpw02 = est02$eldpwv[length(est02$eldpwv)]; seldpw02 = est02$seldpwv[length(est02$seldpwv)]
eldpw03 = est03$eldpwv[length(est03$eldpwv)]; seldpw03 = est03$seldpwv[length(est03$seldpwv)]
eldpw04 = est04$eldpwv[length(est04$eldpwv)]; seldpw04 = est04$seldpwv[length(est04$seldpwv)]
eldpw1 = est1$eldpwv[length(est1$eldpwv)]; seldpw1 = est1$seldpwv[length(est1$seldpwv)]
eldpw2 = est2$eldpwv[length(est2$eldpwv)]; seldpw2 = est2$seldpwv[length(est2$seldpwv)]
eldpw3 = est3$eldpwv[length(est3$eldpwv)]; seldpw3 = est3$seldpwv[length(est3$seldpwv)]
eldpw4 = est4$eldpwv[length(est4$eldpwv)]; seldpw4 = est4$seldpwv[length(est4$seldpwv)]
Tab = data.frame(spatial=c(0,0,0,0,1,1,1,1),k=c(1:4,1:4),
                 eldpw=c(eldpw01,eldpw02,eldpw03,eldpw04,eldpw1,eldpw2,eldpw3,eldpw4),
                 seldpw=c(seldpw01,seldpw02,seldpw03,seldpw04,seldpw1,seldpw2,seldpw3,seldpw4))
print(Tab) 
