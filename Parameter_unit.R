#install.packages("np")
require("np") # for the function b.star which is used in the function block.boot
#########################################   INPUT   #########################################
# data     : Error correction data or univariate time series data matrix(row : time )
# lag      : lag of the data
# kn       : minimum number of observations in each regime
# ngird    : number of grid for iterating when finding gamma hat
# r1       : number of simulation
# r2       : number of bootstrapping
# b or bl  : number of block for block bootstrapping
# mod.b    : method of determining block number b 
#           - if "fix", user can choose any positive integer
#           - if "opt", the function automatically choose the optimal number of block from the fucntion b.star
# mod.boot : type of test ("bnd" for block bstp test / "adf" for adf test)

# g1,g2    : threshold parameters(gamma)
# Z        : error correction term matrix(here, it is lagged data matrix y)
# rhs      : deltay_t-1 term data matrix(design matrix)
# lhs      : deltay_t term data matrix

block.boot.detail(data, lag, kn, ngrid, r2, bl, mod.b, mod.boot)





### RBB Sup-W p-value function
# OUTPUT : Detailed results of RBB supwald or adf test according to the options in the function
#          (results of rbbbnd.detail or rbbadf.detail function)
block.boot.detail<-function(data, lag, kn, ngrid, r2, bl, mod.b, mod.boot){ #mod.b : fix or opt / mod.boot : bnd or adf
  if(mod.b=="fix" & mod.boot=="bnd"){
    b<-bl
    print(rbbbnd.detail(data, lag, kn, ngrid, b, r2))
  }
  
  if(mod.b=="fix" & mod.boot=="adf"){
    b<-bl
    print(rbbadf.detail(data, lag ,b, r2))
  }
  
  if(mod.b=="opt" & mod.boot=="bnd"){
    b<-ceiling(b.star(data, kn)[,1])
    print(rbbbnd.detail(data, lag, kn, ngrid, b, r2))
  }
  
  if(mod.b=="opt" & mod.boot=="adf"){
    b<-ceiling(b.star(data, kn)[,1])
    print(rbbadf.detail(data, lag, b, r2))
  }
}


# Individual Residual Block Bootstrap function for bnd, adf case
# OUTPUT - Sup-wald stat., estimation of gamma & coefficients, p-value of RBB sup-wald test,
#          minimum number of observation in each regime, Number of iterating grid and Number of block 
#        - ADF stat., estimation of coefficients, p-value of RBB_ADF test, Number of block
rbbbnd.detail<-function(data, lag, kn, ngrid, b, r2){
  #Construct basic parameters and data matrix
  y <- as.matrix(data,ncol=1)
  n <- nrow(y)
  
  res.bnd <- bndadf.det(y, lag, kn, ngrid)
  uhat<-res.bnd[[1]]           #n-1 vector
  wald<-as.numeric(res.bnd[[2]])
  alp<-as.numeric(res.bnd[[3]])
  est.co<-res.bnd[[4]]
  gamma.hat<-as.numeric(res.bnd[[5]])
  
  #Centring for blocking procedure
  cdy<-matrix(0,n-b,1)
  for(i in 1:(n-b)){
    cdy[i,]<-mean(uhat[i:(i+b-1)]) #uhat : vector form /cdy: n-b x 1
  }
  u.tilda<- uhat - mean(cdy)       #n-p-1 length vector
  
  #Block bootstrapping and calculate p-value
  w.boot <- matrix(0,r2,1)
  k<-floor((n-1)/b)
  l<-k*b+1
  util.mat<-matrix(0,n-b,b)
  for(iu in 1:(n-b)){
    util.mat[iu,]<-u.tilda[iu:(iu+b-1)]
  }
  
  list.boot<-lapply(1:r2,function(xx){
    sam.num<-sample(1:(n-b),k,replace = T)
    u.til.boot<-c()
    for(is in sam.num){
      u.til.boot<-c(u.til.boot,util.mat[is,])   #k*b vector
    }
    
    y.boot <<- cumsum(c(y[1],u.til.boot)) #k*b+1(=l) vector
    
    wald.boot<-bndadf.det(y.boot, lag, kn, ngrid)[[2]]
    return(wald.boot)
  })
  
  for(i in 1:r2){
    w.boot[i,1]<-list.boot[[i]]
  }
  pval.b<-1-sum(w.boot<=wald)/r2
  
  #generate results
  gamma.res<-c( paste0(round(gamma.hat[1],4),"(",gamma.hat[3],")"),paste0(round(gamma.hat[2],4),"(",gamma.hat[4],")") )
  gamma.mat<-data.frame();gamma.mat[1,1]<-gamma.res[1];gamma.mat[2,1]<-gamma.res[2]
  rownames(gamma.mat)<-c("gamma1","gamma2");colnames(gamma.mat)<-"value(qunatile)"
  if(p>0){
    rownames(est.co)<-c("Coefficient of ECT1","Coefficient of ECT2","Constants",
                          paste( rep("Coefficients of lag",p) ,rep(1:p) ) )
  }else{
    rownames(est.co)<-c("Coefficient of ECT1","Coefficient of ECT2","Constants")
  }
  colnames(est.co)<-"variable"
  conseq <- list("Sup-Wald stat."=wald,
                 "Gamma hat and its quantile in ECT values"=gamma.mat,
                 "Estimation of Coefficients"=est.co,
                 "P-value of RBB Sup-wald stat."=pval.b,
                 "Minimum number of observation in each regime"=kn,
                 "Number of iterating grid"=ngrid,
                 "Number of block"=b)
  
  return(conseq)
}
rbbadf.detail<-function(data, lag ,b, r2){
  #Construct model data matrix and main parameters
  y<-as.matrix(data,ncol=1)
  n<-nrow(y)
  p<-lag
  
  #ADF results
  list.res<-adf.det(y,p)
  dyht<-list.res[[1]]
  t<-as.numeric(list.res[[2]])
  para.hat<-list.res[[3]]
  
  #Block bootstrapping and calculate p-value
  t.boot <- matrix(0,r2,1)
  k<-floor((n-1)/b)
  l<-k*b+1
  util.mat<-matrix(0,n-b,b)
  for(iu in 1:(n-b)){
    util.mat[iu,]<-dyht[iu:(iu+b-1)]
  }
  
  list.boot<-lapply(1:r2,function(xx){
    sam.num<-sample(1:(n-b),k,replace = T)
    u.til.boot<-c()
    for(is in sam.num){
      u.til.boot<-c(u.til.boot,util.mat[is,])
    }
    
    y.boot <- y[1]
    for(j in 1:(l-1)){
      y.boot<-c(y.boot,sum(u.til.boot[1:j]))
    }
    return(adf.det(y.boot, p)[[2]])
  })
  
  for(i in 1:r2){
    t.boot[i,1]<-list.boot[[i]]
  }
  pval<-sum(t.boot<=t)/r2
  conseq <- list("ADF stat."=t,
                 "Estimation of Coefficients"=para.hat,
                 "p-value of RBB_ADF test"=pval,
                 "Number of block"=b)
  return(conseq)
}



# Sup-Wald statistics
# OUTPUT : list of the followings
#        - uhat(residual) : n-1 length vector
#        - sup-wald statistics : scalar
#        - estimation of ECT terms : vector
#        - estimation of coefficients : vector
#        - estimation of gamma and its quantile : vector
bndadf.det<-function(y, lag, m, ngrid){        # model : delta(t)=mu+delta lags
  #Construct model data matrix
  y <- as.matrix(y,ncol=1)
  n <<- nrow(y)
  p <<-lag
  dy_lhs<-embed(diff(y),(p+1))[,1]    #n-p-1 x 1
  dy_lag<-embed(diff(y),(p+1))[,-1]   #n-p-1 x p
  dy_rhs<-cbind(1,dy_lag)             #n-p-1 x (p+1)
  
  if(p>0){
    z<- y[-c(1:p,n),]                   #n-p-1 x 1
  }else{
    z<- y[-n,]
  }  
  
  #Setting grid for finding gamma hat
  ys <- max(abs(y))
  sprt <- seq(from=-ys,by=2*ys/ngrid,length.out=ngrid+1)
  m1<-length(sprt)
  W<-matrix(0,m1,m1)
  
  for(i in 1:m1){
    g1<-sprt[i]
    for(j in i:m1){
      g2<-sprt[j]
      result <- try(find.w2.det(g1, g2, z, dy_rhs, dy_lhs, m), silent = TRUE)
      if (inherits(result, "try-error")) {
        W[i, j] <-NA
      }else{
        W[i, j] <- result$Wald
      }
    }
  }
  sup <- max(W,na.rm = T)
  argsup<- which(W==sup,arr.ind = T)
  gmL<-sprt[argsup[1]]
  gmU<-sprt[argsup[2]]
  res<-find.w2.det(gmL,gmU,z,dy_rhs,dy_lhs,m)
  ap<-res$alpha_hat
  para.hat<-res$para.hat
  
  quan1<-round(sum(ifelse(z <= gmL, 1, 0))/length(z),3) 
  quan2<-round(sum(ifelse(z <= gmU, 1, 0))/length(z),3)
  #derive RBB residual u.hat
  ylag<-y[-nrow(y),]
  u <- diff(y)-ap[1]*ylag*ifelse(ylag<=gmL,1,0)-ap[2]*ylag*ifelse(ylag>gmU,1,0)
  
  return(list("RBB_residual"=u,"sup-Wald"=sup,"alpha hat"=ap,"Estimation of coefficients"=para.hat,"gamma hat"=c(gmL, gmU, quan1, quan2)))
} 

# ADF test statistics
# OUTPUT : list of the followings
#        - adf model residual : n-1 length vector
#        - ADF statistic      : scalar
#        - estimation of coefficients : vector
adf.det<-function(data, lag){                # model : y(t)=mu+y(t-1)+delta_y of lags
  #Construct model data matrix
  y<-as.matrix(data,ncol=1)
  p<-lag
  n<-nrow(y)
  y_lhs<-embed(y,(p+2))[,1]           #n-p-1 x 1
  y_rhs<-embed(y,(p+2))[,2]           #n-p-1 x 1
  dy_lag<-embed(diff(y),(p+1))[,-1]   #n-p-1 x p
  x <- cbind(1,y_rhs,dy_lag)          #n-p-1 x (p+2)
  
  #Estimation from LS method
  ix <- solve(crossprod(x))         #(p+2)x(p+2)
  bt <- ix%*%t(x)%*%y_lhs           #betahat : (p+2) x 1
  eht <- y_lhs-x%*%bt               #residual of the model above : (n-p-1) x 1
  sig <- crossprod(eht)/nrow(eht)   #est. of variance of residual : 1x1
  
  t <- (bt[2]-1)/sqrt(sig*ix[2,2])  #bt[2] : y(t-1)'s coefficient estimation : 1x1
  uht <- y[-1,]-bt[2]*y[-n,]        #delta(t)-bt[2]*delta(t-1) : n-1 x 1
  uht <- uht - mean(uht)            #u centering
  
  return(list("ADF residual" = uht,"ADF statistic"=t,"para.hat"=bt))
}

# Find Wald statistics
# OUTPUT : list of estimation of coefficients and wald statistics 
find.w2.det<-function(g1, g2, Z, rhs, lhs, kn){
  if(sum(ifelse(Z <= g1, 1, 0))>kn & sum(ifelse(Z > g2, 1, 0))>kn ){
    Z.gam<-cbind(ifelse(Z <= g1, 1, 0) * Z,ifelse(Z > g2, 1, 0) * Z) #n-p-1 x 2
    
    X.tot<-cbind(Z.gam, rhs)                                         #n-p-1 x 3+p
    bhat<-solve(crossprod(X.tot))%*%t(X.tot)%*%lhs                   #3+p x 3+p
    
    Mx<-diag(1,n-p-1)-rhs%*%solve(crossprod(rhs))%*%t(rhs)           #n-p-1 x n-p-1
    My <- Mx%*%lhs                                                   #n-p-1 x 1
    Mz<-Mx%*%Z.gam                                                   #n-p-1 x 2
    invMz<-solve(crossprod(Mz))                                      #2x2
    ahat<-invMz%*%t(Mz)%*%My                                         #2x1
    res<-My-Mz%*%ahat                                                #n-p-1 x 1
    sig <- crossprod(res)/nrow(res)                                  #1x1
    wij <- t(ahat)%*%t(Mz)%*%Mz%*%ahat/sig                           #1x1
  }
  return(list(alpha_hat=ahat,Wald=wij,para.hat=bhat))
}


