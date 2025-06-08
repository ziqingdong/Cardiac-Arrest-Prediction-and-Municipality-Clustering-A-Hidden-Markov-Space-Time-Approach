MLE_HMM <- function(Y,X=NULL,M,k,eta=NULL,Be=NULL,piv=NULL,Pi=NULL,
                           zip=FALSE,rand.start=FALSE,sta=FALSE){

# EM algorithm for MLE estimation of the hidden Markov spatial-temporal model
#
# INPUT:
# Y          = matrix of responses (n x TT)
# X          = covariate at response level (ncov x n x TT)
# M          = matrix of offsets (n x TT)
# k          = number of latent states
# eta        = initial value of eta parameter (1)
# Be         = initial value of beta parameters (k x 1+ncov)
# piv        = initial value of initial probability vector (k)
# Pi         = initial value of transition matrix (k x k)
# rand.start = random start
# sta        = initial distribution equal to the stationary distribution
#
# OUTPUT:
# eta      = parameters affecting the conditional zero response probability (nit x 2) (for ZIP model)
# Be0      = starting value of beta parameters without states
# Be       = parameters affecting the conditional response probabilities (k x (1+ncov))
# Pc       = matrix of posterior probabilities (k x n x TT)
# pi       = initial probability vector (k)
# Pi       = transition matrix (k x k)
# lk       = final log-likelihood
# aic.     = value of AIC
# bic.     = value of BIC
# Ul       = local decoding (n x TT)

#---- preliminaries ----
  n = as.integer(nrow(Y)); TT = as.integer(ncol(Y))
  Y = matrix(as.double(Y),n,TT)
  k = as.integer(k)
  MM = array(0,c(k,k-1,k))
  for(u in 1:k) MM[,,u] = diag(k)[,-u]
  if(!is.null(X)){
    X1 = X; X1[is.na(X1)] = 0
  }
  M1 = matrix(as.integer(pmax(1,M)),n,TT)
  if(!is.null(cN)){
    cmax = ncol(cN)
    cN = matrix(as.integer(cN),n,cmax)
    nnv = as.integer(rowSums(cN>0))
    C = matrix(as.integer(0),n,n)
    for(i in 1:n) if(nnv[i]>0) C[i,cN[1:nnv[i]]] = 1
    nnv = as.integer(nnv)
  }
  if(is.null(X)) ncov = 0 else ncov = dim(X)[1]
  mmax = as.integer(max(M))
  nT = as.integer(n*TT)
  Y2 = as.integer(c(Y))
  if(zip) lYf = lgamma(Y2+1)
  lM = log(M)
  M2 = c(M1)
  lM2 = log(M2)
  if(!is.null(X)) X2 = t(matrix(X1,ncov,nT))
  if(!is.null(Mn)){
    Mn = as.matrix(Mn)
    TT1 = ncol(Mn)
    lMn = log(Mn)
  }
  YY2 = rep(Y2,each=k)
  if(is.null(X)){
    XX2 = matrix(1,nT,1)%x%diag(k)
  }else{
    XX2 = cbind(matrix(1,nT,1)%x%diag(k),X2%x%matrix(1,k,1))
  }
  if(zip) MM2 = rep(M2,each=k)
  LM2 = rep(lM2,each=k)
  if(sta & k>1){
    ind = k*(0:(k-1))+(1:k)
    I = diag(k^2); I = I[,-ind]
  }

#---- initial values ----
  if(is.null(Be)){
    if(k==1){
      if(zip){
        est = zeroinfl(Y2~X2|M2,offset = lM2)
        Be = est$coefficients$count
        eta = est$coefficients$zero
      }else{
        if(is.null(X)){
          est = glm(Y2~1,offset = lM2,family=poisson())
        }else{
          est = glm(Y2~X2,offset = lM2,family=poisson())
        }
        Be = est$coefficients
      }
      Be1 = Be[1]
      if(!is.null(X)) Be2 = Be[-1]
      piv = 1; Pi = 1
    }else{
      if(rand.start){
        if(zip){
          est = zeroinfl(Y2~X2|M2,offset = lM2)
          V1 = vcov(est)[1:(ncov+1),1:(ncov+1)]
          Be = matrix(rmvnorm(k,est$coefficients$count,k*V1),k,ncov+1)
          V2 = vcov(est)[(ncov+2):(ncov+3),(ncov+2):(ncov+3)]
          eta = rmvnorm(k,est$coefficients$zero,k*V2)
        }else{
          if(is.null(X)){
            est = glm(Y2~1,offset = lM2,family=poisson())
          }else{
            est = glm(Y2~X2,offset = lM2,family=poisson())
          }
          Be = matrix(rmvnorm(k,est$coefficients,k*vcov(est)),k,ncov+1)
        }
        Be1 = Be[,1]
        if(!is.null(X)){
          Be2 = colMeans(as.matrix(Be[,-1]))
          Be = cbind(int=Be1,cov=rep(1,k)%o%Be2)
        }
        if(is.null(piv)){
          if(is.null(Pi)){
            Pi = matrix(runif(k^2),k); Pi = (1/rowSums(Pi))*Pi
          }
          if(sta){
            Tmp = Pi
            for(it1 in 1:1000) Tmp = Tmp%*%Pi
            piv = colMeans(Tmp)
          }else{
            piv = runif(k); piv = piv/sum(piv)
          }
        }
      }else{
        Int = matrix(1,TT)%x%diag(n)
        if(zip){
          if(ncov==0){
            est = zeroinfl(Y2~Int-1|M2,offset = lM2)
          }else{
            est = zeroinfl(Y2~cbind(Int,X2)-1|M2,offset = lM2)
          }
          out = kmeans(est$coefficients$count[1:n],k,nstart = 1000)
          eta = est$coefficients$zero
          Be2 = est$coefficients$count[-(1:n)]
        }else{
          if(ncov==0){
            est = glm(Y2~Int-1,offset = lM2,family=poisson())
          }else{
            est = glm(Y2~cbind(Int,X2)-1,offset = lM2,family=poisson())
          }
          out = kmeans(est$coefficients[1:n],k,nstart = 1000)
          Be2 = est$coefficients[-(1:n)]
        }
        if(any(is.na(Be2))) Be2[is.na(Be2)] = 0
        ord = order(out$centers)
        Be1 = out$centers[ord]
        Be = cbind(int=Be1,cov=rep(1,k)%o%Be2)
        clus = rep(0,n)
        for(u in 1:k) clus[out$cluster==ord[u]] = u
        U = matrix(as.integer(clus),n,TT)
        if(is.null(piv)) piv = rep(1/k,k)
        if(is.null(Pi)) Pi = (1+diag(k)*k)/(2*k)
      }
    }
  }

#---- EM algorithm ----
# compute log-likelihood
  Pc = A = array(0,c(k,n,TT))
  for(i in 1:n) for(t in 1:TT){
    if(is.null(X)){
      la = exp(lM[i,t] + Be1)
    }else{
      la = exp(lM[i,t] + Be1 + c(X[,i,t]%*%Be2))
    }
    if(zip){
      eta0 = eta[1]+M[i,t]*eta[2]
      p0 = 1/(1+exp(-eta0)); p1 = 1/(1+exp(eta0))
      Pc[,i,t] = p0*(Y[i,t]==0)+p1*dpois(Y[i,t],la)
    }else{
      Pc[,i,t] = dpois(Y[i,t],la)
    }
  }
  Pc = pmax(Pc,10^-10)
  if(k==1){
    pv = apply(Pc[1,,],1,prod)
  }else{
    Pit = t(Pi)
    for(i in 1:n){
      A[,i,1] = piv*Pc[,i,1]
      for(t in 2:TT) A[,i,t] = (Pit%*%A[,i,t-1])*Pc[,i,t]
    }
    pv = colSums(A[,,TT])
  }
  lk = sum(log(pv))
  cat("------------|-------------|-------------|-------------|\n")
  cat("     it     |      k      |      lk     |     diff    |\n")
  cat("------------|-------------|-------------|-------------|\n")
  cat(sprintf("%11g",c(0,k,lk)),"\n",sep=" | ")
  if(k>1){
    if(sta){
      th = NULL
      for(u in 1:k) th = c(th,log(Pi[u,-u]/Pi[u,u]))
    }
    lko = lk; it = 0
    while(((lk-lko)/abs(lko)>10^-10 & it<1000) | it==0){
      lko = lk; it = it+1
# E-step
      B = Pp1 = array(0,c(k,n,TT))
      Pp2 = matrix(0,k,k)
      for(i in 1:n){
        B[,i,TT] = 1
        Pp1[,i,TT] = A[,i,TT]/pv[i]
        Pp2 = Pp2+(A[,i,TT-1]%o%(B[,i,TT]*Pc[,i,TT]))*Pi/pv[i]
        for(t in (TT-1):1){
          B[,i,t] = Pi%*%(Pc[,i,t+1]*B[,i,t+1])
          Pp1[,i,t] = (A[,i,t]*B[,i,t])/pv[i]
          if(t>1) Pp2 = Pp2+(A[,i,t-1]%o%(B[,i,t]*Pc[,i,t]))*Pi/pv[i]
        }
      }
# M-step
      if(sta){
        th = optim(th,lksta,method="Nelder-Mead",control=list(fnscale=-1,reltol=10^-8),
                   pp1=rowSums(Pp1[,,1]),Pp2=Pp2,I=I)$par
        Pi = matrix(exp(I%*%th),k,k,byrow = TRUE)
        Pi = (1/rowSums(Pi))*Pi
        Tmp = Pi
        for(it1 in 1:1000) Tmp = Tmp%*%Pi
        piv = colMeans(Tmp)
      }else{
        piv = rowSums(Pp1[,,1])/sum(Pp1[,,1])
        Pi = (1/rowSums(Pp2))*Pp2
      }
      if(zip){
        YY2 = as.integer(YY2)
        est = zeroinfl(YY2~XX2-1|MM2,offset = LM2,weights = c(Pp1))
        eta = est$coefficients$zero
        Be1 = est$coefficients$count[1:k]
        Be2 = est$coefficients$count[-(1:k)]
      }else{
        est = glm(YY2~XX2-1,offset = LM2,family=poisson(),weights = c(Pp1))
        Be1 = est$coefficients[1:k]
        Be2 = est$coefficients[-(1:k)]
      }
      Be = cbind(Be1,rep(1,k)%o%Be2)
# compute log-likelihood
      Pc = A = array(0,c(k,n,TT))
      for(i in 1:n) for(t in 1:TT){
        if(is.null(X)){
          la = exp(lM[i,t] + Be1)
        }else{
          la = exp(lM[i,t] + Be1 + c(X[,i,t]%*%Be2))
        }
        if(zip){
          eta0 = eta[1]+M[i,t]*eta[2]
          p0 = 1/(1+exp(-eta0)); p1 = 1/(1+exp(eta0))
          Pc[,i,t] = p0*(Y[i,t]==0)+p1*dpois(Y[i,t],la)
        }else{
          Pc[,i,t] = dpois(Y[i,t],la)
        }
      }
      Pit = t(Pi)
      for(i in 1:n){
        A[,i,1] = piv*Pc[,i,1]
        for(t in 2:TT) A[,i,t] = (Pit%*%A[,i,t-1])*Pc[,i,t]
      }
      pv = colSums(A[,,TT])
      lk = sum(log(pv))
      if(it%%10==0) cat(sprintf("%11g",c(it,k,lk,lk-lko)),"\n",sep=" | ")
    }
    if(it%%10>0) cat(sprintf("%11g",c(it,k,lk,lk-lko)),"\n",sep=" | ")
  }
  cat("------------|-------------|-------------|-------------|\n")
  
#---- output ----
  np = k^2-1+k+ncov+2*zip
  if(sta) np = np-(k-1)
  aic = -2*lk+2*np
  bic = -2*lk+log(n)*np
  if(k==1){
    Pp1 = array(1,c(1,n,TT))
    Ul = matrix(1,n,TT)
  }else{
    Ul = apply(Pp1,c(2,3),which.max)
  }
  if(zip){
    out = list(eta=eta,Be=Be,piv=piv,Pi=Pi,Pp1=Pp1,lk=lk,np=np,aic=aic,bic=bic,Ul=Ul)
  }else{
    out = list(Be=Be,piv=piv,Pi=Pi,Pp1=Pp1,lk=lk,np=np,aic=aic,bic=bic,Ul=Ul)
  }
  return(out)

}