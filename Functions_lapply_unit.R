#install.packages("np")
require("np") # for the function b.star which is used in the function block.boot
#########################################   INPUT   #########################################
# data     : Error correction data or univariate time series data matrix(row : time )
# lag      : lag of the data
# kn       : minimum number of observations in each regime
# ngird    : number of grid for iterating when finding gamma hat
# r1       : number of simulation
# r2       : number of bootstrapping
# sam      : sample size
# b or bl  : number of block for block bootstrapping
# mod.b    : method of determining block number b 
#           - if "fix", user can choose any positive integer
#           - if "opt", the function automatically choose the optimal number of block from the fucntion b.star
# mod.boot : type of test ("bnd" for block bstp test / "adf" for adf test)

# g1,g2    : threshold parameters(gamma)
# Z        : error correction term matrix(here, it is lagged data matrix y)
# rhs      : deltay_t-1 term data matrix(design matrix)
# lhs      : deltay_t term data matrix

# gm       : gamma value(symmetric gamma : (-theta,theta) )
# ap       : scalar value alpha(coefficient of ECM term)
# mu       : level of the model

# rho      : AR coefficient value(scalar)
# tht      : MA coefficient value(scalar)

# mod.type : "size" or "power" (specify the type of test)


###############################################################################################
###############################   Data simulation procedure    ################################

#  simulation code for all the cases
#  OUTPUT : Rejection ratio of two tests according to the input information 
simulation2<-function(r1, r2, sam, gm0, rho0, tht0, b, lag, mu, mod.type){ #mod.type =="size" or "power"

 if(mod.type=="size"){
   res5<-res10<-as.data.frame(matrix(0,2,5)) 
   rownames(res5)<-rownames(res10)<-c("supW","ADF")
   colnames(res5)<-colnames(res10)<-c("(0,0)","(-0.5,0)","(0.5,0)","(0,-0.5)","(0,0.5)")
   
   list.res<-list()
   list.res<-list( res5,  res10)
   p5<-paste("@@Result of NULL : size(sign.lvl : 0.05) ( n =", sam, "/ alpha = 0 )")
   p10<-paste("@@Result of NULL : size(sign.lvl : 0.1) ( n =", sam, "/ alpha = 0 )")
   names(list.res)<-c(p5,p10)
   
   n<-sam
   kn<-10
   ngrid<-round(n/10)+20
   
   for(am in 1:length(rho0)){
     rho<-rho0[am];  tht<-tht0[am]
     sim.res<-simu.null(r1, r2, n, rho, tht, b, lag, mu, kn, ngrid)
     for(i in 1:2){
       list.res[[2]][i,am]<-sim.res[[i]][1]
       list.res[[1]][i,am]<- sim.res[[i]][2]
     }
   }

   res5.a<-rbind(c("(0,0)","(-0.5,0)","(0.5,0)","(0,-0.5)","(0,0.5)"),list.res[[1]])
   res10.a<-rbind(c("(0,0)","(-0.5,0)","(0.5,0)","(0,-0.5)","(0,0.5)"),list.res[[2]]) 
   rownames(res5.a)<-rownames(res10.a)<-c("Rho&Tht","supW","ADF")
   
   cat("<< Rejection frequency >>\n", file="Test of size_unit result.csv")
   cat(p5, file="Test of size_unit result.csv",append = T,"\n")
   write.table(as.matrix(res5.a),"Test of size_unit result.csv",append = T,sep = ",",col.names = F)
   cat(p10, file="Test of size_unit result.csv",append = T,"\n")
   write.table(as.matrix(res10.a),"Test of size_unit result.csv",append = T,sep = ",",col.names = F)
   return(list.res)
  }else{
    
    n<-sam
    ap<-(-0.1)
    kn<-10
    ngrid<-round(n/10)+20
    
    res5.4<-res5.8<-res10.4<-res10.8<-as.data.frame(matrix(0,2,5)) 
    rownames(res5.4)<-rownames(res5.8)<-rownames(res10.4)<-rownames(res10.8)<-c("supW","ADF")
    colnames(res5.4)<-colnames(res5.8)<-colnames(res10.4)<-colnames(res10.8)<-c("(0,0)","(-0.5,0)","(0.5,0)","(0,-0.5)","(0,0.5)")
    
    list.res2<-list()
    list.res2<-list( res5.4, res10.4, res5.8,  res10.8 )
    p5.4<-paste("@@Result of Alternative : power(sign.lvl : 0.05) ( alpha =", ap, "/ n =", sam, "/ gamma =", gm0[1],")")
    p10.4<-paste("@@Result of Alternative : power(sign.lvl : 0.1) ( alpha =", ap, "/ n =", sam, "/ gamma =", gm0[1],")")
    p5.8<-paste("@@Result of Alternative : power(sign.lvl : 0.05) ( alpha =", ap, "/ n =", sam, "/ gamma =", gm0[2],")")
    p10.8<-paste("@@Result of Alternative : power(sign.lvl : 0.1) ( alpha =", ap, "/ n =", sam, "/ gamma =", gm0[2],")")
    names(list.res2)<-c(p5.4,p10.4,p5.8,p10.8)
    
    for(ig in 1:length(gm0)){
      gm<-gm0[ig]
      for(am in 1:length(rho0)){
        rho<-rho0[am];  tht<-tht0[am]
        sim.res<-simu.alt(r1, r2, n, gm, ap, rho, tht, b, lag, mu, kn, ngrid)
        for(i in 1:2){
          list.res2[[(ig-1)*2+1]][i,am]<-sim.res[[i]][2]
          list.res2[[(ig-1)*2+2]][i,am]<-sim.res[[i]][1]
        }
      }
    }
    
    res5.4a<-rbind(c("(0,0)","(-0.5,0)","(0.5,0)","(0,-0.5)","(0,0.5)"),list.res2[[1]])
    res10.4a<-rbind(c("(0,0)","(-0.5,0)","(0.5,0)","(0,-0.5)","(0,0.5)"),list.res2[[2]])
    res5.8a<-rbind(c("(0,0)","(-0.5,0)","(0.5,0)","(0,-0.5)","(0,0.5)"),list.res2[[3]])
    res10.8a<-rbind(c("(0,0)","(-0.5,0)","(0.5,0)","(0,-0.5)","(0,0.5)"),list.res2[[4]])
    
    rownames(res5.4a)<-rownames(res10.4a)<-rownames(res5.8a)<-rownames(res10.8a)<-c("Rho&Tht","supW","ADF")
    
    cat("<< Rejection frequency >>\n", file="Test of power_unit result.csv")
    cat(p5.4, file="Test of power_unit result.csv",append = T,"\n")
    write.table(as.matrix(res5.4a),"Test of power_unit result.csv",append = T,sep = ",",col.names = F)
    cat(p10.4, file="Test of power_unit result.csv",append = T,"\n")
    write.table(as.matrix(res10.4a),"Test of power_unit result.csv",append = T,sep = ",",col.names = F)
    cat(p5.8, file="Test of power_unit result.csv",append = T,"\n")
    write.table(as.matrix(res5.8a),"Test of power_unit result.csv",append = T,sep = ",",col.names = F)
    cat(p10.8, file="Test of power_unit result.csv",append = T,"\n")
    write.table(as.matrix(res10.8a),"Test of power_unit result.csv",append = T,sep = ",",col.names = F)
    
    return(list.res2)
  }
}

#  simulation code for each case
#  OUTPUT : Rejection ratio of RBB sup-wald/ADF tests under NULL
simu.null<-function(r1, r2, n, rho, tht, b, lag, mu, kn, ngrid){  
  wn<-adf1<-matrix(0,r1,1)
  lapply(1:r1, function(i){
    xs<-as.matrix(dgp_lin(n, rho, tht, mu))
    wn[i,]  <<- rbb.bnd(xs, lag, kn, ngrid, b, r2)
    adf1[i,]<<- rbb.adf(xs, lag ,b, r2)
  })

  rej.wn10 <- sum(wn<0.1)/r1
  rej.adf10<- sum(adf1<0.1)/r1 
  rej.wn5 <- sum(wn<0.05)/r1
  rej.adf5<- sum(adf1<0.05)/r1 
  
  list.res<-list( Boot.sup = c(rej.wn10,rej.wn5), ADF = c(rej.adf10,rej.adf5) )
  return(list.res)
}

#  OUTPUT : Rejection ratio of RBB sup-wald/ADF tests under Alternative
simu.alt<-function(r1, r2, n, gm, ap, rho, tht, b, lag, mu, kn, ngrid){
  wn<-adf1<-matrix(0,r1,1)
  lapply(1:r1, function(i){
    xs<-as.matrix(dgp_thr(n, rho, tht, gm, ap, mu))
    wn[i,]  <<- rbb.bnd(xs, lag, kn, ngrid, b, r2)
    adf1[i,]<<- rbb.adf(xs, lag ,b, r2)
  })

  rej.wn10 <- sum(wn<0.1)/r1
  rej.adf10<- sum(adf1<0.1)/r1 
  rej.wn5 <- sum(wn<0.05)/r1
  rej.adf5<- sum(adf1<0.05)/r1 
  
  list.res<-list( Boot.sup = c(rej.wn10,rej.wn5), ADF = c(rej.adf10,rej.adf5) )
  return(list.res)
}  




###############################################################################################
#############################   Data construction for simulation   ############################

# 1. dxt=mu+alpha%*%ZU ZL + u / u~ARMA(1,1)
# OUTPUT : Band-Threshold based random generated data matrix
dgp_thr<-function(n, rho, tht, gm, ap, mu){
  e <- matrix(rnorm(n+5,0,1),ncol=1)
  u <- x <- matrix(0,n+5,1)
  u[1,] <- x[1,] <- e[1,]
  
  for(j in 2:(n+5)){
    u[j,] <- rho*u[(j-1),] + e[j,] + tht*e[(j-1),]
    dx<- mu + ap*x[(j-1),]*ifelse(abs(x[(j-1),])> gm,1,0) + u[j,]
    x[j,] <- x[(j-1),]+dx 
  }
  x <- x[-(1:5),]
  return(x)
}
#dgp_thr(100,0,0.5,4,-0.1,0)

# 2. dxt=mu+u / u~ARMA(1,1)
# OUTPUT : ARMA (linear) based random generated data matrix
dgp_lin<-function(n, rho, tht, mu){
  e <- matrix(rnorm(n+5,0,1),ncol=1)
  u <- x <- matrix(0,n+5,1)
  u[1,] <- x[1,] <- e[1,]
  
  for(j in 2:(n+5)){
    u[j,] <- rho*u[(j-1),] + e[j,] + tht*e[(j-1),]
    dx<- mu + u[j,]
    x[j,] <- x[(j-1),]+dx 
  }
  x <- x[-(1:5),]
  return(x)
}
#dgp_lin(100,0,0.5,0)




###############################################################################################
###############################   Derive test_statistics   ####################################

### 1. RBB Sup-W p-value function
block.boot<-function(data, lag, kn, ngrid, r2, bl, mod.b, mod.boot){ #mod.b : fix or opt / mod.boot : bnd or adf
  if(mod.b=="fix" & mod.boot=="bnd"){
    b<-bl
    return(rbb.bnd(data, lag, kn, ngrid, b, r2))
  }
  
  if(mod.b=="fix" & mod.boot=="adf"){
    b<-bl
    return(rbb.adf(data, lag ,b, r2))
  }
  
  if(mod.b=="opt" & mod.boot=="bnd"){
    b<-ceiling(b.star(data, kn)[,1])
    return(rbb.bnd(data, lag, kn, ngrid, b, r2))
  }
  
  if(mod.b=="opt" & mod.boot=="adf"){
    b<-ceiling(b.star(data, kn)[,1])
    return(rbb.adf(data, lag, b, r2))
  }
}


# 1) Individual Residual Block Bootstrap function for bnd, adf case
# OUTPUT : p-value of RBB sup-wald test / RBB ADF test
rbb.bnd<-function(data, lag, kn, ngrid, b, r2){
  #Construct basic parameters and data matrix
  y <- as.matrix(data,ncol=1)
  n <- nrow(y)
  
  res.bnd <- bndadf(y, lag, kn, ngrid)
  uhat<-res.bnd[[1]]           #n-1 vector
  wald<-res.bnd[[2]]
  
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
    
    wald.boot<-bndadf(y.boot, lag, kn, ngrid)[[2]]
    return(wald.boot)
  })
  
  for(i in 1:r2){
    w.boot[i,1]<-list.boot[[i]]
  }
  return(1-sum(w.boot<=wald)/r2)
}
rbb.adf<-function(data, lag ,b, r2){
  #Construct model data matrix and main parameters
  y<-as.matrix(data,ncol=1)
  n<-nrow(y)
  
  #ADF results
  res.adf <- adf(y, lag)
  dyht<-res.adf[[1]]           #n-1 vector
  t<-as.numeric(res.adf[[2]])
  
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
    
    y.boot <<- cumsum(c(y[1],u.til.boot)) #k*b+1(=l) vector
    t.b <- as.numeric(adf(y.boot, lag)[[2]])
    return(t.b)
  })
  
  for(i in 1:r2){
    t.boot[i,1]<-list.boot[[i]]
  }
  return(sum(t.boot<=t)/r2)
}


# 2) Derive Sup-Wald statistics
# OUTPUT : list of the followings
#        - uhat(residual) : n-1 length vector
#        - sup-wald statistics : scalar
bndadf<-function(y, lag, m, ngrid){        # model : delta(t)=mu+delta lags
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
      result <- try(find.w2(g1, g2, z, dy_rhs, dy_lhs, m), silent = TRUE)
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
  res<-find.w2(gmL,gmU,z,dy_rhs,dy_lhs,m)
  ap<-res$alpha_hat
  
  #derive RBB residual u.hat
  ylag<-y[-nrow(y),]
  u <- diff(y)-ap[1]*ylag*ifelse(ylag<=gmL,1,0)-ap[2]*ylag*ifelse(ylag>gmU,1,0)
  
  return(list("RBB_residual"=u,"sup-Wald"=sup))
} 

# 3) Find Wald statistics
# OUTPUT : list of estimation of coefficients and wald statistics 
find.w2<-function(g1, g2, Z, rhs, lhs, kn){
  if(sum(ifelse(Z <= g1, 1, 0))>kn & sum(ifelse(Z > g2, 1, 0))>kn ){
    Z.gam<-cbind(ifelse(Z <= g1, 1, 0) * Z,ifelse(Z > g2, 1, 0) * Z) #n-p-1 x2
    Mx<-diag(1,n-p-1)-rhs%*%solve(crossprod(rhs))%*%t(rhs)           #n-p-1 x n-p-1
    My <- Mx%*%lhs                                                   #n-p-1 x 1
    Mz<-Mx%*%Z.gam                                                   #n-p-1 x 2
    invMz<-solve(crossprod(Mz))                                      #2x2
    ahat<-invMz%*%t(Mz)%*%My                                         #2x1
    res<-My-Mz%*%ahat                                                #n-p-1 x 1
    sig <- crossprod(res)/nrow(res)                                  #1x1
    wij <- t(ahat)%*%t(Mz)%*%Mz%*%ahat/sig                           #1x1
  }
  return(list(alpha_hat=ahat,Wald=wij))
}


### 2. ADF test statistics
# OUTPUT : list of the followings
#        - adf model residual : n-1 length vector
#        - ADF statistic : scalar
adf<-function(data, lag){                # model : y(t)=mu+y(t-1)+delta_y of lags
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
  
  return(list("ADF residual" = uht,"ADF statistic"=t))
}






