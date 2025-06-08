MCMC_HMM <- function(Y,X=NULL,Z=NULL,M,cN=NULL,k,nit=100000,taueta=0.01,taube = 0.01,tauga = 0.25,
                            taude = 0.25,regany = 50,eta=NULL,Be=NULL,Ga=NULL,De=NULL,U=NULL,
                            same=FALSE,Xn=NULL,Zn=NULL,Mn=NULL,zip=FALSE,rand.start=FALSE,output=FALSE,
                            label.switch=TRUE,timesn=NULL){

# Augmented MCMC algorithm for MCMC Bayesian estimation of the hidden Markov spatial-temporal model
#
# INPUT:
# Y        = matrix of responses (n x TT)
# X        = covariate at response level (ncov1 x n x TT)
# Z        = covariate at latent level (np2 x n x TT)
# M        = matrix of offsets (n x TT)
# k        = number of latent states
# cN       = matrix of neighborhood for each unit (n x cmax) - if missing an LM model is fitted
# taueta    = sd of the proposal for eta (for ZIP model)
# taube    = sd of the proposal for beta
# tauga    = sd of the proposal for gamma
# taude    = sd of the proposal for delta
# regany   = register and iteration any
# eta      = initial value of eta parameter (1)
# Be       = initial value of beta parameters (k x 1+np1), np1 = ncov1+1
# Ga       = initial value of gamma parameters (k-1 x (1+np2+k))
# De       = initial value of delta parameters (k-1 x (1+np2+k) x k)
# U        = initial value of beta parameters (n x TT)
# same     = TRUE for same regression parameters across states
# Xn       = covariate at response level for prediction (np1 x n x TTn)
# Zn       = covariate at latent level for prediction (np2 x n x TTh)
# Mn       = matrix of offsets for predictions (n x TT1)
# rand.start = for random starting values
# output   = to have the full output with estimates at any registered iteration
# label.switch = to require label switching
# timesn   = labels of the new times
#
# OUTPUT:
# Eta      = parameters affecting the conditional zero response probability (nit x 2) (for ZIP model)
# Be0      = starting value of beta parameters
# mBe      = estimate of beta parameters
# seBe     = standard error for beta parameters
# Ga0      = starting value of gamma parameters
# mGa      = estimate of gamma parameters
# seGa     = standard error for gamma parameters
# De0      = starting value of delta parameters
# mDe      = estimate of delta parameters
# seDe     = standard error for delta parameters
# U0       = starting value of U
# MU       = MAP assignement of U
# mpiv     = estimate of the average initial probabilities
# sepiv    = standard error for the average initial probabilities
# mPi      = estimate of the average transition probabilities
# sePi     = standard error for the average transition probabilities
# mYh      = prediction of Y
# seYh     = standard error for prediction of Y
# lklprv0  = vector of loglikelihood+logpraio (nit) for the only component involving beta
# MUh      = prediction of U
# lklprv   = vector of loglikelihood+logpraio (nit)
# eldpwv   = vector of eldpw values
# seldpwv  = vector of eldpw s.e.
# acceta   = acceptances of eta (for ZIP model)
# accbe    = acceptances of Be
# accU     = acceptances of U
# accga    = acceptances of Ga
# accde    = acceptances of De
# ESSBe    = effective sample size for beta
# ESSGa    = effective sample size for gamma
# ESSDe    = effective sample size for delts
# *optional if output=TRUE
# BE       = parameters affecting the conditional response probabilities (k x (1+np1) x nit)
# UU       = matrix of latent state (n x TT x nit)
# GA       = parameters affecting the initial probabilities (k-1 x (1+np2+k) x nit)
# DE       = parameters affecting the transition probabilities (k-1 x (np2+k) x k x nit)
# UUn      = predicted latent states (n x TT1 x nit1)
# YYn      = predicted responses (n x TT1 x nit1)
# ELDPW    = single components of eldpw (n x TT x nit1)
# LPR      = posterior expected value for any observation and interation (nT x nit1)
  
t1 = proc.time()

#---- hyperparameters ----
  sieta = sibe = siga = side = sqrt(1000)

#---- preliminaries ----
  sites = dimnames(Y)[[1]]
  times = dimnames(Y)[[2]]
  Xnames = c("(intercept)",dimnames(X)[[1]])
  Znames = c("(intercept)",dimnames(Z)[[1]])
  if(!is.null(cN)) Znames = c(Znames,paste("state_nei",2:k,sep="."))
  n = as.integer(nrow(Y)); TT = as.integer(ncol(Y))
  Y = matrix(as.double(Y),n,TT)
  k = as.integer(k)
  MM = array(0,c(k,k-1,k))
  for(u in 1:k) MM[,,u] = diag(k)[,-u]
  if(is.null(X)) X = array(0,c(0,n,TT))
  X1 = X; X1[is.na(X1)] = 0
  M1 = matrix(as.integer(pmax(1,M)),n,TT)
  if(!is.null(cN)){
    cmax = ncol(cN)
    cN = matrix(as.integer(cN),n,cmax)
    nnv = as.integer(rowSums(cN>0))
    C = matrix(as.integer(0),n,n)
    for(i in 1:n) if(nnv[i]>0) C[i,cN[1:nnv[i]]] = 1
    nnv = as.integer(nnv)
  }
  ncov1 = dim(X)[1]
  if(is.null(Z)) Z = array(0,c(0,n,TT))
  ncov2 = dim(Z)[1]
  np1 = as.integer(1+ncov1)
  if(is.null(cN)) np2 = as.integer(1+ncov2) else np2 = as.integer(ncov2+k)
  mmax = as.integer(max(M))
  nT = as.integer(n*TT)
  Y2 = c(Y)
  lYf = lgamma(Y2+1)
  M2 = c(M1)
  lM2 = log(M2)
  X2 = t(matrix(X1,ncov1,nT))
  if(!is.null(Mn)){
    Mn = as.matrix(Mn)
    TT1 = ncol(Mn)
    lMn = log(Mn)
  }

#---- starting values ----
  if(is.null(U)){
    if(rand.start){
      U = matrix(sample(1:k,nT,rep=TRUE),n,TT)
    }else{
      Int = matrix(1,TT)%x%diag(n)
      if(zip){
        if(ncov1==0){
          est = zeroinfl(Y2~Int-1|M2,offset = lM2)
        }else{
          est = zeroinfl(Y2~cbind(Int,X2)-1|M2,offset = lM2)
        }
        out = kmeans(est$coefficients$count[1:n],k,nstart = 1000)
      }else{
        if(ncov1==0){
          est = glm(Y2~Int-1,offset = lM2,family=poisson())
        }else{
          est = glm(Y2~cbind(Int,X2)-1,offset = lM2,family=poisson())
        }
        out = kmeans(est$coefficients[1:n],k,nstart = 1000)
      }
      U = matrix(out$cluster,n,TT)
    }
  }
  U = matrix(as.integer(U),n,TT)
  if(is.null(Be)){
    Int = matrix(0,nT,k)
    Int[cbind(1:nT,c(U))] = 1
    if(ncov1==0){
      est = glm(Y2~Int-1,offset = lM2,family=poisson())
    }else{
      est = glm(Y2~cbind(Int,X2)-1,offset = lM2,family=poisson())
    }
    if(zip){
      if(ncov1==0){
        est = zeroinfl(Y2~Int-1|M2,offset = lM2)
      }else{
        est = zeroinfl(Y2~cbind(Int,X2)-1|M2,offset = lM2)
      }
      eta = est$coefficients$zero
      Be1 = est$coefficients$count[1:k]; Be2 = est$coefficients$count[-(1:k)]
    }else{
      if(ncov1==0){
        est = glm(Y2~Int-1,offset = lM2,family=poisson())
      }else{
        est = glm(Y2~cbind(Int,X2)-1,offset = lM2,family=poisson())
      }
      Be1 = est$coefficients[1:k]; Be2 = est$coefficients[-(1:k)]
    }
    Be = cbind(Be1,rep(1,k)%o%Be2)
  }else{
    if(is.vector(Be)) Be = t(Be)
    Be1 = Be[,1]; Be2 = colMeans(Be[,-1,drop=FALSE])
  }
  if(is.null(Ga)){
    if(rand.start) Ga = matrix(rnorm(np2*(k-1),0,tauga),k-1,np2) else Ga = matrix(0,k-1,np2)
  }
  if(is.null(De)){
    if(rand.start) De = array(rnorm(np2*(k-1)*k,0,taude),c(k-1,np2,k)) else De = array(0,c(k-1,np2,k))
  }

# label switching
  if(k>1 & label.switch){
    pit = 0
    ind = order(Be1)
    if(!all(ind==1:k)){
      Ga0 = Ga
      cat("label switch of starting values\n")
      ind = order(Be1)
      Be1 = Be1[ind]
      Be = Be[ind,,drop=FALSE]
      U0 = U
      for(u in 1:k) U[U0==ind[u]] = u
      if(!is.null(cN)){
        Tmp = cbind(0,Ga[,(ncov2+2):np2])[,ind,drop=FALSE]
        Ga[,1] = Ga[,1]+Tmp[,1]
        Ga[,(ncov2+2):np2] = Tmp[,-1,drop=FALSE]-Tmp[,1]
        for(up in 1:k){
          Tmp = cbind(0,De[,(ncov2+2):np2,up])[,ind,drop=FALSE]
          De[,1,up] = De[,1,up]+Tmp[,1]
          De[,(ncov2+2):np2,up] = Tmp[,-1,drop=FALSE]-Tmp[,1]
        }
      }
      MGa = MM[,,1]%*%Ga
      MGa = MGa[ind,,drop=FALSE]
      MGa = MGa-rep(1,k)%o%MGa[1,]
      Ga = MGa[-1,,drop=FALSE]
      MDe = MDes = array(0,c(k,np2,k))
      for(up in 1:k) MDe[,,up] = matrix(MM[,,up],k,k-1)%*%matrix(De[,,up],k-1,np2)
      MDe = MDe[ind,,ind,drop=FALSE]
      for(up in 1:k){
        MDe[,,up] = MDe[,,up]-rep(1,k)%o%MDe[up,,up]
        De[,,up] = MDe[-up,,up]
      }
    }
  }
  dimnames(Be) = list(state=1:k,covariate=Xnames)
  dimnames(Ga) = dimnames(De) = NULL
  Be0 = Be; Ga0 = Ga; De0 = De; U0 = U

#---- MCMC algorithm ----
  acceta = accbe = accU = accga = accde = 0
  nit1 = nit/regany
  Eta = matrix(0,nit1,2)
  BE = array(0,c(k,np1,nit1))
  dimnames(BE) = list(state=1:k,covariate=Xnames,iteraton=(1:(nit/regany))*regany)
  UU = array(0,c(n,TT,nit1))
  dimnames(UU) = list(site=sites,time=times,iteraton=(1:(nit/regany))*regany)
  GA = array(0,c(k-1,np2,nit1))
  if(k>1) dimnames(GA) = list(state=2:k,covariate=Znames,iteraton=(1:(nit/regany))*regany)
  DE = array(0,c(k-1,np2,k,nit1))
  if(k>1) dimnames(DE) = list(state=2:k,covariate=Znames,prev_state=1:k,iteration=(1:(nit/regany))*regany)
  if(k>1){
    Mpiv = matrix(0,nit1,k)
    dimnames(Mpiv) = list(iteration=(1:(nit/regany))*regany,state=1:k)
    MPi = array(0,c(k,k,nit1))
    dimnames(MPi) = list(prev_state=1:k,state=1:k,iteration=(1:(nit/regany))*regany)
  }
  if(!is.null(Mn)){
    UUn = array(0,c(n,TT1,nit1))
    dimnames(UUn) = list(site=sites,time=timesn,iteraton=(1:(nit/regany))*regany)
    YYn = array(0,c(n,TT1,nit1))
    dimnames(YYn) = list(site=sites,time=timesn,iteraton=(1:(nit/regany))*regany)
  }
  lklprv0 = lklprv = rep(0,nit1)
  names(lklprv0) = names(lklprv) = (1:(nit/regany))*regany
  eldpwv = seldpwv = rep(0,nit1)
  names(eldpwv) = names(seldpwv) = (1:(nit/regany))*regany
  provavv = rep(0,nit1)
  ELDPW = array(0,c(n,TT,nit1))
  dimnames(ELDPW) = list(site=sites,time=times,iteration=(1:(nit/regany))*regany)
  LPR = matrix(0,nT,nit1)
  dimnames(LPR) = list(site.time=NULL,iteration=(1:(nit/regany))*regany)
  eldpwaic = seldpwaic = 0
  t0 = proc.time()[3]
  # print(table(U))
  if(zip){
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
    cat("     it     |      k      |    acceta   |    accbe    |    accU     |    accga    |    accde    |    lklpr    |  eldpwaic   |  seldpwaic  |     time    |\n")
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
  }else{
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
    cat("     it     |      k      |    accbe    |    accU     |    accga    |    accde    |    lklpr    |  eldpwaic   |  seldpwaic  |     time    |\n")
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
  }
  for(it in (-nit/4):nit){

# Update eta parameter
    if(zip){
      etas = eta+rnorm(2,0,taueta)
      out = .Fortran("lk_be4_STZ",eta,etas,Be,Be,U,lM2,X2,Y2,n,TT,k,np1,nT,lkv=rep(0,k),
                     lksv=rep(0,k),lYf)
      lk = sum(out$lkv); lks = sum(out$lksv)
      lpr = sum(dnorm(eta,0,sieta,log=TRUE))
      lprs = sum(dnorm(etas,0,sieta,log=TRUE))
      al = min(1,exp(lprs+lks-lpr-lk))
      if(runif(1)<al){
        eta = etas
        acceta = acceta+1
        lk = lks; lpr = lprs
      }
      lklpr = lpr
    }else{
      lklpr = 0
    }

# Update beta parameters
    if(same){
      Be1s = Be1+rnorm(k,0,taube)
      Be2s = Be2+rnorm(np1-1,0,taube)
      Bes = cbind(Be1s,rep(1,k)%o%Be2s)
      if(zip){
        out = .Fortran("lk_be4_STZ",eta,eta,Be,Bes,U,lM2,X2,Y2,n,TT,k,np1,nT,lkv=rep(0,k),
                       lksv=rep(0,k),lYf)
      }else{
        out = .Fortran("lk_be3_STP",Be,Bes,U,lM2,X2,Y2,n,TT,k,np1,nT,lkv=rep(0,k),
                       lksv=rep(0,k),lYf)
      }
      lk = sum(out$lkv); lks = sum(out$lksv)
      lpr = sum(dnorm(Be1,0,sibe,log=TRUE))+ sum(dnorm(Be2,0,sibe,log=TRUE))
      lprs = sum(dnorm(Be1s,0,sibe,log=TRUE))+ sum(dnorm(Be2s,0,sibe,log=TRUE))
      al = min(1,exp(lprs+lks-lpr-lk))
      if(runif(1)<al){
        Be1 = Be1s; Be2 = Be2s; Be = Bes
        accbe = accbe+1
        lk = lks; lpr = lprs
      }
      lklpr = lklpr+lpr
      if(k==1) lklpr = lklpr+lk
    }else{
      Bes = Be+rnorm(k*np1,0,taube)
      if(zip){
        out = .Fortran("lk_be3_STZ",eta,eta,Be,Bes,U,lM2,X2,Y2,n,TT,k,np1,nT,lkv=rep(0,k),
                       lksv=rep(0,k),lYf)
      }else{
        out = .Fortran("lk_be3_STP",Be,Bes,U,lM2,X2,Y2,n,TT,k,np1,nT,lkv=rep(0,k),lksv=rep(0,k))
      }
      for(u in 1:k){
        beu = Be[u,]; beus = Bes[u,]
        lk = out$lkv[u]; lks = out$lksv[u]
        lpr = sum(dnorm(beu,0,sibe,log=TRUE))
        lprs = sum(dnorm(beus,0,sibe,log=TRUE))
        al = min(1,exp(lprs+lks-lpr-lk))
        if(runif(1)<al){
          Be[u,] = beus
          accbe = accbe+1/k
          lk = lks; lpr = lprs
        }
        lklpr = lklpr+lk+lpr
      }
    }
    if(is.nan(lklpr)){
      print(1)
      browser()
    }
    lklpr0 = lklpr

# update gamma parameters
    if(k>1) for(u in 1:(k-1)){
      gau = Ga[u,]
      gaus = gau+rnorm(np2,0,tauga)
      Gas = Ga; Gas[u,] = gaus
      MGa = MM[,,1]%*%Ga
      MGas = MM[,,1]%*%Gas
      if(is.null(cN)){
        out = .Fortran("lk_ga2_TP",MGa,k,np2,MGas,n,U,TT,Z,lk=0,lks=0)
      }else{
        out = .Fortran("lk_ga2",MGa,k,np2,MGas,cN,cmax,n,U,TT,Z,nnv,lk=0,lks=0)
      }
      lk = out$lk; lks = out$lks
      lpr = sum(dnorm(gau,0,siga,log=TRUE))
      lprs = sum(dnorm(gaus,0,siga,log=TRUE))
      al = min(1,exp(lprs+lks-lpr-lk))
      if(!is.nan(al)) if(runif(1)<al){
        Ga = Gas
        accga = accga+1/(k-1)
        lk = lks; lpr = lprs
      }
      # lklpr = lklpr+lk+lpr
      lklpr = lklpr+lpr
    }
    if(is.nan(lklpr) | is.infinite(lklpr)){
      print(3)
      browser()
    }

# update delta parameters
    if(k>1) for(up in 1:k){
      for(u in 1:(k-1)){
        deupu = De[u,,up]
        deupus = deupu+rnorm(np2,0,taude)
        Des = De; Des[u,,up] = deupus
        MDe = MDes = array(0,c(k,np2,k))
        for(u1 in 1:k){
          MDe[,,u1] = matrix(MM[,,u1],k,k-1)%*%matrix(De[,,u1],k-1,np2)
          MDes[,,u1] = matrix(MM[,,u1],k,k-1)%*%matrix(Des[,,u1],k-1,np2)
        }
        if(is.null(cN)){
          out = .Fortran("lk_de2_TP",MDe,k,np2,MDes,n,U,TT,Z,up,lk=0,lks=0)
        }else{
          out = .Fortran("lk_de2",MDe,k,np2,MDes,C,n,U,TT,Z,nnv,up,lk=0,lks=0)
        }
        lk = out$lk; lks = out$lks
        lpr = sum(dnorm(deupu,0,side,log=TRUE))
        lprs = sum(dnorm(deupus,0,side,log=TRUE))
        al = min(1,exp(lprs+lks-lpr-lk))
        if(!is.nan(al)) if(runif(1)<al){
          De = Des
          accde = accde+1/(k*(k-1))
          lk = lks; lpr = lprs
        }
        # lklpr = lklpr+lk+lpr
        lklpr = lklpr+lpr
        if(is.nan(lklpr) | is.infinite(lklpr)){
          print(4)
          browser()
        }
      }
    }

# update latent variables
    if(k>1){
      MGa = MM[,,1]%*%Ga
      MDe = array(0,c(k,np2,k))
      for(u in 1:k) MDe[,,u] = matrix(MM[,,u],k,k-1)%*%matrix(De[,,u],k-1,np2)
      U0 = U
      if(is.null(cN)){
        if(zip){
          out = .Fortran("update_u6_TZ",eta,Be,k,np1,MGa,np2,MDe,U=U,lM2,X2,nT,Z,Y2,n,TT,
                         accU=accU,LR=array(0,c(k,n,TT)),lk=0,lYf)
        }else{
          out = .Fortran("update_u5_TP",Be,k,np1,MGa,np2,MDe,U=U,lM2,X2,nT,Z,Y2,n,TT,
                         accU=accU,LR=array(0,c(k,n,TT)),lk=0,lYf)
        }
      }else{
        if(zip){
          out = .Fortran("update_u6_STZ",eta,Be,k,np1,MGa,np2,MDe,cN,cmax,U=U,lM2,X2,nT,Z,Y2,n,TT,
                         nnv,accU=accU,LR=array(0,c(k,n,TT)),lk=0,lYf)
        }else{
          out = .Fortran("update_u5_STP",Be,k,np1,MGa,np2,MDe,cN,cmax,U=U,lM2,X2,nT,Z,Y2,n,TT,
                         nnv,accU=accU,LR=array(0,c(k,n,TT)),lk=0,lYf)
        }
      }
      lklpr = lklpr+out$lk
    }
    if(is.nan(lklpr)){
      print(2)
      browser()
    }

# label switching
    if(k>1 & label.switch){
      ind = order(Be1)
      if(!all(ind==1:k)){
          if(pit%%15==0){
          if(pit>0) cat("\n")
          cat("label switch at iteration",it)
        }else{
          cat("",it)
        }
        pit = pit+1
        ind = order(Be1)
        Be1 = Be1[ind]
        Be = Be[ind,,drop=FALSE]
        dimnames(Be) = list(state=1:k,covariate=Xnames)
        U0 = U
        for(u in 1:k) U[U0==ind[u]] = u
        if(!is.null(cN)){
          Tmp = cbind(0,Ga[,(ncov2+2):np2])[,ind,drop=FALSE]
          Ga[,1] = Ga[,1]+Tmp[,1]
          Ga[,(ncov2+2):np2] = Tmp[,-1,drop=FALSE]-Tmp[,1]
          for(up in 1:k){
            Tmp = cbind(0,De[,(ncov2+2):np2,up])[,ind,drop=FALSE]
            De[,1,up] = De[,1,up]+Tmp[,1]
            De[,(ncov2+2):np2,up] = Tmp[,-1,drop=FALSE]-Tmp[,1]
          }
        }
        MGa = MM[,,1]%*%Ga
        MGa = MGa[ind,,drop=FALSE]
        MGa = MGa-rep(1,k)%o%MGa[1,]
        Ga = MGa[-1,,drop=FALSE]
        MDe = MDes = array(0,c(k,np2,k))
        for(up in 1:k) MDe[,,up] = matrix(MM[,,up],k,k-1)%*%matrix(De[,,up],k-1,np2)
        MDe = MDe[ind,,ind,drop=FALSE]
        for(up in 1:k){
          MDe[,,up] = MDe[,,up]-rep(1,k)%o%MDe[up,,up]
          De[,,up] = MDe[-up,,up]
        }
      }
    }

# average initial and transition probabilities
    if(k>1 & it>0 & it%%regany==0){
      MGa = MM[,,1]%*%Ga
      if(is.null(cN)){
        mpiv = .Fortran("mpiv_TP",MGa,k,np2,n,TT,Z,mpiv=rep(0,k))$mpiv
      }else{
        mpiv = .Fortran("mpiv2",MGa,k,np2,cN,cmax,n,U,TT,Z,nnv,mpiv=rep(0,k))$mpiv
      }
      MDe = array(0,c(k,np2,k))
      for(u1 in 1:k) MDe[,,u1] = matrix(MM[,,u1],k,k-1)%*%matrix(De[,,u1],k-1,np2)
      mPi = matrix(0,k,k)
      if(is.null(cN)){
        for(up in 1:k) mPi[up,] = .Fortran("mPi_TP",MDe,k,np2,n,TT,Z,up,mpiv=rep(0,k))$mpiv
      }else{
        for(up in 1:k) mPi[up,] = .Fortran("mPi2",MDe,k,np2,C,n,U,TT,Z,nnv,up,mpiv=rep(0,k))$mpiv
      }
    }

# prediction
    if(!is.null(Mn) & it>0 & it%%regany==0){
      if(k==1){
        Un = matrix(1,n,TT1)
      }else{
        MDe = MDes = array(0,c(k,np2,k))
        for(u1 in 1:k) MDe[,,u1] = matrix(MM[,,u1],k,k-1)%*%matrix(De[,,u1],k-1,np2)
        Un = matrix(U[,TT],n,TT1)
        for(t in 1:TT1){
          if(t==1) u0 = U[,TT] else u0 = Un[,t-1]
          if(is.null(cN)){
            out = .Fortran("predict_u_TP",k,np2,MDe,u0,Zn[,,t],n,u=u0)
          }else{
            out = .Fortran("predict_u_STP",k,np2,MDe,cN,cmax,u0,Zn[,,t],n,nnv,u=u0)
          }
          Un[,t] = out$u
        }
      }
      Yn = matrix(0,n,TT1)
      for(t in 1:TT1){
        Tmp = Be[Un[,t],]
        if(is.null(Xn)){
          la = exp(lMn[,t]+Tmp)
        }else{
          if(ncov1==1)   #FB changed here
            la = exp(lMn[,t]+rowSums(cbind(1,Xn[,,1])*Tmp))
          else
            la = exp(lMn[,t]+rowSums(cbind(1,t(Xn[,,1]))*Tmp))
        }
        if(zip){
          eta0 = eta[1]+Mn[,t]/1000*eta[2]
          p0 = 1/(1+exp(-eta0))
          Yn[,t] = (runif(n)>p0)*rpois(n,la)
        }else{
          Yn[,t] = rpois(n,la)
        }
      }
    }

# record iteration
    if(it>0 & it%%regany==0){
      it1 = it/regany
      if(zip) Eta[it1,] = eta
      BE[,,it1] = Be
      UU[,,it1] = U
      GA[,,it1] = Ga
      DE[,,,it1] = De
      if(k>1){
        Mpiv[it1,] = mpiv
        MPi[,,it1] = mPi
      }
      if(!is.null(Mn)){
        UUn[,,it1] = Un
        YYn[,,it1] = Yn
      }
      lklprv0[it1] = lklpr0
      lklprv[it1] = lklpr
      llav = rep(0,nT)
      for(t in 1:TT){
        ind = (t-1)*n+(1:n)
        llav[ind] = lM2[ind]+rowSums(cbind(1,X2[ind,])*Be[U[,t],])
      }
      lav = exp(llav)
      if(zip){
        eta0 = eta[1]+M2/1000*eta[2]
        p0 = 1/(1+exp(-eta0)); p1 = 1/(1+exp(eta0))
        LPR[,it1] = log(p0*(Y2<0.5)+p1*dpois(Y2,lav))
      }else{
        LPR[,it1] = dpois(Y2,lav,log=TRUE)
      }
      Tmp = LPR[,1:it1,drop=FALSE]
      tmpv = apply(Tmp,1,max)
      tmp1v = log(rowMeans(exp(Tmp-tmpv)))+tmpv
      tmp2v = rowMeans((Tmp-rowMeans(Tmp))^2)
      provavv[it1] = sum(tmp1v)
      ELDPW[,,it1] = tmp1v-tmp2v
      eldpwaic = sum(ELDPW[,,it1],na.rm=TRUE)
      seldpwaic = sqrt(nT*var(as.vector(ELDPW[,,it1]),na.rm=TRUE))
      eldpwv[it1] = eldpwaic
      seldpwv[it1] = seldpwaic
      if(it%%(10*regany)==0){
        tmp1 = lklprv[1:it1]-max(lklprv[1:it1])
        tmp2 = eldpwv[1:it1]-max(eldpwv[1:it1])
        tmp3 = provavv[1:it1]-max(provavv[1:it1])
        plot(tmp1,ylim = c(min(c(tmp1,tmp2)),max(c(tmp1,tmp2))),type="l"); lines(c(1,it1),c(tmp1[it1],tmp1[it1]),col="grey")
        lines(tmp2,type="l",col=2); lines(c(1,it1),c(tmp2[it1],tmp2[it1]),col="grey")
        lines(tmp3,type="l",col=3); lines(c(1,it1),c(tmp3[it1],tmp3[it1]),col="grey")
      }
    }

# display acceptance rate
    if(it%%1000==0){
      if(k>1 & label.switch){
        if(pit>0) cat("\n")
        pit = 0
      }
      tt = proc.time()[3]-t0
      it2 = it+nit/4
      if(zip){
        cat(sprintf("%11g",c(it,k,acceta/it2,accbe/it2,accU/it2,accga/it2,accde/it2,lklpr,eldpwaic,seldpwaic,tt/it2*1000)),
            "\n",sep=" | ")
        
      }else{
        cat(sprintf("%11g",c(it,k,accbe/it2,accU/it2,accga/it2,accde/it2,lklpr,eldpwaic,seldpwaic,tt/it2*1000)),
            "\n",sep=" | ")
      }
      # print(table(U))
    }
  }
  if(zip){
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
  }else{
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
  }

# final output
  mBe = apply(BE,c(1,2),mean)
  seBe = apply(BE,c(1,2),sd)
  out = list(Be0=Be0,mBe=mBe,seBe=seBe)
  if(k>1){
    mGa = apply(GA,c(1,2),mean)
    seGa = apply(GA,c(1,2),sd)
    mDe = apply(DE,c(1,2,3),mean)
    seDe = apply(DE,c(1,2,3),sd)
    Mode <- function(x){
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    MU = apply(UU,1:2,Mode)
    mpiv = colMeans(Mpiv)
    sepiv = apply(Mpiv,2,sd)
    mPi = apply(MPi,1:2,mean)
    sePi = apply(MPi,1:2,sd)
    out = c(out,list(Ga0=Ga0,mGa=mGa,seGa=seGa,De0=De0,mDe=mDe,seDe=seDe,U0=U0,MU=MU,mpiv=mpiv,
                     sepiv=sepiv,mPi=mPi,sePi=sePi))
  }
  if(!is.null(Mn)){
    mYh = apply(YYn,1:2,mean)
    seYh = apply(YYn,1:2,sd)
    out = c(out,list(mYh=mYh,seYh=seYh))
    if(k>1){
      MUh = apply(UUn,c(1,2),Mode)
      out$MUh = MUh
    }
  }
  out = c(out,list(lklprv0=lklprv0,lklprv=lklprv,eldpwv=eldpwv,seldpwv=seldpwv,accbe=accbe/it2))
  if(k>1) out = c(out,list(accU=accU/it2,accga=accga/it2,accde=accde/it2))
  ESSBe = apply(BE,1:2,ESS)
  out$ESSBe = ESSBe
  if(k>1){
    ESSGa = apply(GA,1:2,ESS)
    ESSDe = apply(DE,1:3,ESS)
    out = c(out,list(ESSGa=ESSGa,ESSDe=ESSDe))
  }
  if(output){
    out$BE = BE
    if(k>1){
      out = c(out,list(GA=GA,DE=DE,UU=UU,Mpiv,MPi))
    }
    if(!is.null(Mn)){
      if(k>1) out$UUn = UUn
      out$YYn = YYn
    }
    out = c(out,list(ELDPW=ELDPW,LPR=LPR))
  }
  if(zip) out = c(out,list(Eta=Eta,acceta=acceta))
  return(out)
  
}