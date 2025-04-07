
#####generate model#######################################################################
#inverse 
matpower <- function(a,b){
  a <- round((a+t(a))/2,7)
  tem <- eigen(a)
  return(tem$vectors%*%diag((tem$values)^b)%*%t(tem$vectors))
}

#Square root
Sqrt <- function(Omega) {
  EgSig0 <- eigen(Omega)
  EgVec <-EgSig0$vectors
  Gamma <- EgVec %*% diag(sqrt(EgSig0$values)) %*% t(EgVec)
  return(Gamma)
}

#variable selection ,high dimension p>n
Sample_COV <- function(X){
  n <- nrow(X)
  Q <- diag(n)-1/n*rep(1, n)%*%t(rep(1, n))
  out <- t(X)%*%Q%*%X/(n-1)
  return(out)
}

Gendata <- function(n, p, d, m, mu, struc,rho, eps, Xd, Sc){
  
  #covariance structure
  if (struc == "AR1") {
    #Autoregressive structure
    Sigma <-  diag(1,p)
    for (i in 1:p){
      for(j in 1:p){
        Sigma[i,j] <- rho^(abs(i-j))
      }
    }
  }
  if (struc == "CS") {
    # #Compound symmetry structure
    Sigma <- diag(p)
    for(i in 1:p){
      for(j in 1:p){
        if(i!=j)   Sigma[i,j] <- rho
      }
    }
  }
  
  # # #generate covarites
  mean <- rep(0,p)
  dat <- mvrnorm(n,mean,Sigma)
  if (Xd == "norm"){
    X <- dat
  }
  if (Xd == "t3"){
    mean <- rep(0,p)
    X <- rmt(n, mean,  Sigma, 3)
  }
  if (Xd == "t5"){
    mean <- rep(0,p)
    X <- rmt(n, mean, Sigma, 5)
  }
  if (Xd == "mixed3normCORt"){
    pp <- floor(p/3)
    meanpp <- rep(0,pp)
    Sigmapp <- diag(1,pp)
    for (i in 1:pp){
      for(j in 1:pp){
        Sigmapp[i,j] <- rho^(abs(i-j))
      }
    }
    X1 <- mvrnorm(n,meanpp,Sigmapp)   
    
    X2 <- mvrnorm(n,meanpp,diag(1,pp))  
    
    
    mean3 <- rep(0,p-2*pp)
    Sigma3 <- diag(1,p-2*pp)     #diag
    dat_temp3 <- mvrnorm(n,mean3,Sigma3)
    dat3 <- dat_temp3%*%diag(1/sqrt(diag(Sigma3)))
    pv3 <- 2*pnorm(abs(dat3), lower.tail = FALSE)
    #t5
    k <- 5
    X3 <- qt(pv3, df=k)/sqrt(k/(k-2))
    
    X <- cbind(X1,X2,X3)
    mean <- apply(X,2,mean)
    Sigma <- Sample_COV(X)   # equal to cov(X)
    
  }
  
  #gernerate error term
  if (eps =="norm")  {ep <- rnorm(n,0,1)/sqrt(1)}
  if (eps =="chisq") {ep <- (rchisq(n,1,0)-1)/sqrt(2)}
  if (eps =="t3")    {ep <- rt(n,3,0)/sqrt(3)}
 
  #model: Y=G(X,beta,g)+ep
  #p<n
  if (Sc == 1){
    
    #ind <- sample(1:p,m,replace = F)   
    ind <- 1:m                          
    ss <- rep(0,p)
    ss[ind]<- 1
    ss[1:floor(m/2)] <- 1
    
    Y <- (mu)*((X%*%(ss))) + 3*ep        #linear model  ,p <n 
    
  }
  if (Sc == 2){
    
    ss <- c(rep(1,m),rep(0,p-m))   
    
    mt <- floor(m/2)
    beta_1 <- as.vector(c(rep(1,mt),rep(0,p-mt)))
    beta_2 <- as.vector(c(rep(0,mt),rep(1,m-mt),rep(0,p-m)))
    
    #low dimensional
    Y <- mu*((abs(X%*%beta_1+3))+exp((X%*%beta_2))) + ep 

  }
  if (Sc == 3){
    
    ss <- c(rep(1,m),rep(0,p-m))
    
    mt <- floor(m/3)
    beta_1 <- 1*as.vector(c(rep(1,mt),rep(0,p-(mt))))
    beta_2 <- 1*as.vector(c(rep(0,mt),rep(1,mt),rep(0,p-2*mt)))
    beta_3 <- 1*as.vector(c(rep(0,2*mt),rep(1,m-2*mt),rep(0,p-m)))
    
    Y <- mu*((X%*%beta_1)+(((X%*%beta_2)+3)^2)+exp(1*(X%*%beta_3))) + ep 
    
  }
  #p>n
  if (Sc == 11){
    
    ind <- 1:m                         
    ss <- rep(0,p)
    ss[ind]<- 1

    Y <- (mu)*1*exp(5+((X%*%(ss)))) + ep        #nonlinear model ,p>n       
  }
  if (Sc == 22){
    
    ss <- c(rep(1,m),rep(0,p-m))   

    mt <- floor(m/2)
    beta_1 <- 1*as.vector(c(rep(1,mt),rep(0,p-mt)))
    beta_2 <- 1*as.vector(c(rep(0,mt),rep(1,m-mt),rep(0,p-m)))
    
    Y <- mu*(2*((X%*%beta_1)+0))+3*exp(1*(X%*%beta_2)) +1*ep
  }
  if (Sc == 33){
    
    ss <- c(rep(1,m),rep(0,p-m))
    mt <- floor(m/3)
    beta_1 <- 1*as.vector(c(rep(1,mt),rep(0,p-(mt))))
    beta_2 <- 1*as.vector(c(rep(0,mt),rep(1,mt),rep(0,p-2*mt)))
    beta_3 <- 1*as.vector(c(rep(0,2*mt),rep(1,m-2*mt),rep(0,p-m)))

    Y <- mu*((X%*%beta_1)+abs((X%*%beta_2)+5)+exp(1*(X%*%beta_3))) + ep 

  }

  X <- t((t(X)-apply(X,2,mean))/ apply(X, 2, sd))

  return(list(X = X, Y = Y, ss = ss, ep = ep, Sigma = Sigma, mean = mean)) 
}


#####estimate the dimension reduction matrix B and kernel matrix M#########################
discretize <- function(y,h){
  n=length(y);m=floor(n/h)
  y=y+.00001*mean(y)*rnorm(n)
  yord = y[order(y)]
  divpt=numeric();for(i in 1:(h-1)) divpt = c(divpt,yord[i*m+1])
  y1=rep(0,n);y1[y<divpt[1]]=1;y1[y>=divpt[h-1]]=h
  for(i in 2:(h-1)) y1[(y>=divpt[i-1])&(y<divpt[i])]=i
  return(y1)}
SIR <- function(d,x,y){
  p <- ncol(x)
  n <- nrow(x)
  h <- 6
  ytype <-"continuous"
  N <- var(x)
  signrt <- matpower(N,-1/2)
  xc <- t(t(x)-apply(x,2,mean))
  xst <- xc%*%signrt
  if(ytype=="continuous") ydis <- discretize(y,h)
  if(ytype=="categorical") ydis <- y
  yless <- ydis
  ylabel <- numeric()
  for(i in 1:n) {if(var(yless)!=0) {
    ylabel <- c(ylabel,yless[1])
    yless <- yless[yless!=yless[1]]}}
  ylabel <- c(ylabel,yless[1])
  prob <- numeric()
  exy <- numeric()
  for(i in 1:h) prob <- c(prob,length(ydis[ydis==ylabel[i]])/n) 
  for(i in 1:h) exy <- rbind(exy,apply(xst[ydis==ylabel[i],],2,mean)) 
  sirmat <- t(exy)%*%diag(prob)%*%exy
  #M <- sirmat
  sig <- matpower(N,1/2)
  M <- sig%*%sirmat%*%sig
  B <- (eigen(sirmat)$vectors[,1:d])  #Li, 2018
  B <- as.matrix(B)
  
  return(list(B = B, M = sirmat, N = N, mat = sirmat))
}
SAVE <- function(d,x,y){
  p <- ncol(x)
  n <- nrow(x)
  h <- 6
  ytype <- "continuous"
  N <- var(x)
  signrt <- matpower(N,-1/2)
  xc <- t(t(x)-apply(x,2,mean))
  xst <- xc%*%signrt
  if(ytype=="continuous") ydis <- discretize(y,h)
  if(ytype=="categorical") ydis <- y
  yless <- ydis
  ylabel <- numeric()
  for(i in 1:n) {if(var(yless)!=0) {
    ylabel <- c(ylabel,yless[1])
    yless <- yless[yless!=yless[1]]}
  }
  
  ylabel <- c(ylabel,yless[1])
  prob <- numeric()
  for(i in 1:h) prob <- c(prob,length(ydis[ydis==ylabel[i]])/n)
  vxy  <-  array(0,c(p,p,h))
  for(i in 1:h) vxy[,,i]<-var(xst[ydis==ylabel[i],])
  savemat <- 0
  for(i in 1:h){
    savemat <- savemat+prob[i]*(vxy[,,i]-diag(p))%*%(vxy[,,i]-diag(p))
  }
  
  sig <- matpower(N,1/2)
  M <- sig%*%savemat%*%sig
  B <- (eigen(savemat)$vectors[,1:d])  #Li, 2018,??S(Y|Z)
  B <- as.matrix(B)
  
  return(list(B = B, M = savemat, N = N, mat = savemat))
}

##different transformation function for response and centralization response to meet E(f(Y))=0##
discretize_trans <- function(y,H,trans){
  
  n <- length(y)
  
  yless <- dr.slices(y,H)$slice.indicator 
  indexlist <- vector("list",H)  #vector("list",length)
  for(i in 1:H){
    indexlist[[i]] <- which(yless==i)
  }
  
  if (trans=="indicator"){
    sirdis <- rep(0,n)
    for(i in 1:H){sirdis[indexlist[[i]]] <- i }
    ytemp <- sirdis-mean(sirdis)
    ydis <- rep(0,n)
    for (i in 1:H){ydis[indexlist[[i]]] <- ytemp[indexlist[[i]]]}
  }   #slice function in SIR
  if (trans=="CN") {
    #ydis <- y
    ytemp <- y - mean(y)
    ydis <- rep(0,n)
    for (i in 1:H){ydis[indexlist[[i]]] <- ytemp[indexlist[[i]]]}
  }             #cook and Ni (2006)
  if (trans=="poly"){
    #ydis <- y^2 #2-moment
    #ytemp <- y^2 - mean(y^2)
    ytemp <- y^i - mean(y^i)
    ydis <- rep(0,n)
    for (i in 1:H){ydis[indexlist[[i]]] <- ytemp[indexlist[[i]]]}
  }                                      #Yin and cook (2002),k-th moment 
  
  ylabel <- numeric()
  ytemp <- yless
  for(i in 1:n){
    if(var(ytemp)!=0){
      ylabel <- c(ylabel,ytemp[1])
      ytemp <- ytemp[ytemp!=ytemp[1]]
    }
  }
  ylabel <- c(ylabel,ytemp[1])
  
  return(list(yless = yless, ylabel = ylabel, ydis = ydis,indexlist = indexlist))
}

#######lassso, 1996#################################
LASSO <- function(X, Y, method){
  n <- nrow(X)
  p <- ncol(X)
  
  if (method == "CV"){
    fit_1 <- cv.glmnet(X, Y, family="gaussian")  #10 fold CV
    lambda_select <- fit_1$lambda.min
    #beta_min <- coef(fit_1, s = "lambda.min")
    #beta_lse <- coef(fit_1, s = "lambda.1se")
    
  }
  
  if(method == "AIC"){
    fit_1 <- glmnet(X, Y, family = "gaussian")    #Default is standardize=TRUE.
    k <- fit_1$df
    AIC <- deviance(fit_1)+2*k
    i_min <- which.min(AIC)
    lambda_select <- fit_1$lambda[i_min]
  }
  
  if(method == "BIC"){
    fit_1 <- glmnet(X, Y, family = "gaussian")    #Default is standardize=TRUE.
    k <- fit_1$df
    BIC <- deviance(fit_1)+2*k*log(n)
    i_min <- which.min(BIC)
    lambda_select <- fit_1$lambda[i_min]
  } 
  
  fit <- glmnet(X, Y, family = "gaussian", lambda = lambda_select)
  beta_lasso <- as.numeric(fit$beta[,1])  #coef(fit)
  #beta_lasso <- (coef(fit)[-1])
  
  signal_est <- as.vector(which(beta_lasso!=0))
  signal_num <- length(signal_est)
  
  return(list(beta_lasso = beta_lasso, lambda = lambda_select,
              signal_est = signal_est, signal_num = signal_num))
  
}

####multivariate response linear model estimation by response transformation#####################################
Func_slice <- function(X, Y, H, trans, debiased){
  tuning <- "CV"
  
  n <- nrow(X)
  p <- ncol(X)

  yless <- dr.slices(Y,H)$slice.indicator #slice from small to large
  
  indexlist <- vector("list",H)  #vector("list",length)
  for(i in 1:H){
    indexlist[[i]] <- which(yless==i)
  }
  ylabel <- numeric()
  ytemp <- yless
  for(i in 1:n){
    if(var(ytemp)!=0){
      ylabel <- c(ylabel,ytemp[1])
      ytemp <- ytemp[ytemp!=ytemp[1]]
    }
  }
  ylabel <- c(ylabel,ytemp[1])
  
  B <- matrix(0,p,H)
  for (i in 1:H){
    
    if (trans == "indicator"){ydis <- rep(1,length(Y))}   #SIR,1991
    if (trans == "CIRE"){ydis <- Y}                 #cook and Ni (2006)
    if (trans == "poly"){ydis <- Y^i}  #Y^2             #Yin and cook (2002),k-th moment 
    
    ydis[which(yless!=i)] <- 0
    ydis <- ydis - mean(ydis)
    
    YY <- ydis[indexlist[[i]]]
    XX <- X[indexlist[[i]],]
    
    if (debiased == "OLS"){
      B[,i] <- as.numeric(lm(ydis~X-1)$coefficients)
      sdB <- sqrt(diag(solve( t(X)%*%X)*var(ydis)))
      B[,i] <- B[,i]/sdB
      
    }
    if (debiased == "lasso"){
      fit <- LASSO(X,ydis,method = tuning)
      beta_lasso <- fit$beta_lasso
      B[,i] <- beta_lasso
    }
    if (debiased == "debiased_lasso"){
      
      fit <- LASSO(X,ydis,method = tuning)
      beta_lasso <- fit$beta_lasso
  
      #debiased lasso method 1: theta_lasso + 1/nMt(X)(Y-Xtheta_lasso)
      M <- debiasingMatrix(X, is_wide = TRUE, length(ydis), rows = 1:p) 
      B[,i] <- beta_lasso + (1/(length(ydis)))*M%*%t(X)%*%(ydis-X%*%beta_lasso)
    }
  }
  
  Bp <- eigen(B%*%t(B)/H)$vectors[,1:H]     #SDR subspace
  
  return(list(B = B, Bp = Bp, ylabel = ylabel, ydis = ydis,indexlist = indexlist))
  
}
#two splits, statistics 
MyFunc <-function(dat, H, trans, sp, debiased){
  X <- dat$X
  Y <- dat$Y
  n <- nrow(X)
  p <- ncol(X)
   
  #sample splitting
  n_1 <- floor(n/2)
  n_2 <- n-n_1
  Ind_1 <- sample(1:n,n_1,replace = F)
  Ind_2 <- (1:n)[-c(Ind_1)]
  
  dat_1 <- list(X = X[Ind_1,], Y = Y[Ind_1])
  dat_2 <- list(X = X[Ind_2,], Y = Y[Ind_2])
  

  #process in dat1
  X_1 <- dat_1$X
  Y_1 <- dat_1$Y
  fit_1 <- Func_slice(X_1,Y_1,H,trans,debiased)
  B_1 <- fit_1$B

  # #process in dat2
  X_2 <- dat_2$X
  Y_2 <- dat_2$Y
  fit_2 <- Func_slice(X_2,Y_2,H,trans,debiased)
  B_2 <- fit_2$B

  
  #ranking statistics for variable selection
  if(sp=="full")   {
    #process in full data
    fit_full <- Func_slice(X,Y,H,trans,debiased)
    B_full <- fit_full$B
    Wj <- sapply(1:p,function(j){B_full[j,]%*%B_full[j,]})
    Td <- sapply(1:H,function(x){B_full[,x]%*%B_full[,x]})
  }
  if (sp=="allone"){
    fit_full <- Func_slice(X,Y,H,trans,debiased)
    B_full <- fit_full$B
    Wj <- sapply(1:p,function(j){sum(abs(B_full[j,]))})
    Td <- sapply(1:H,function(x){sum(abs(B_full[,x]))})
  }
  if(sp=="split2") {
    Wj <- sapply(1:p,function(j){B_1[j,]%*%B_2[j,]})
    Td <- sapply(1:H,function(x){B_1[,x]%*%B_2[,x]})
  }
  if(sp=="half")   {
    Wj <- sapply(1:p,function(j){B_1[j,]%*%B_1[j,]})
    Td <- sapply(1:H,function(x){B_1[,x]%*%B_1[,x]})
  }
  
  return(list(Wj = Wj,Td = Td))
}

#variable selection for p < n/(2H)
SDS<- function(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H,trans,sp,debiased,alpha,threshold,seedNo){
  set.seed(123+1000*seedNo)
  
  tic <- proc.time()
  dat <- Gendata(n, p, d, m, mu, struc,rho, eps, Xd, Sc)   #simulation newdata
  tru1 <- which(dat$ss==1)
  
  X <- dat$X
  Y <- dat$Y
  
  temp <- MyFunc(dat,H,trans,sp,debiased)              
  Wj <-temp$Wj
  
  t <- sort(abs(Wj))   #threshold is positive
  
  
  ###FDR 
  
  if (threshold == "-"){
    Ta <- sapply(t,function(x){(0+length(Wj[Wj<=(-x)]))/max(1,length(Wj[Wj>=x]))})
  }  
  if (threshold == "+"){
    Ta <- sapply(t,function(x){(1+length(Wj[Wj<=(-x)]))/max(1,length(Wj[Wj>=x]))}) 
  }  #SDA+, size conservative. The results does not different much with SDA if p sufficient large 
  
  best_L <- min(t[which(Ta<=alpha)])
  
  det1 <- which(Wj>=best_L)                                #the estimated signals
  FDR <- length(setdiff(det1,tru1))/max(1,length(det1))    #under H1, false discovery
  AP <- length(intersect(det1,tru1))/max(1,length(tru1))   #under H1, true discovery
  NFD <- length(setdiff(det1,tru1))                        #number of false discovery
  ND  <- length(det1)                                      #number of discovery
  NN <- length(which(Wj<0))
  
  if (all(tru1%in%det1)==TRUE){Pa <- 1}                   #selection consistency 
  if (all(tru1%in%det1)!=TRUE){Pa <- 0}
  
  # mFDR <- NFD/max(1,ND)            #marginal FDR
  # mAP <- (ND-NFD)/max(1,m)         $marginal AP
  
  tic2 <- proc.time()-tic
  time <- tic2[3]
  
  return( c(FDR,AP,NFD,ND,NN,Pa,time))
  
}
Process_SDS <- function(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H,trans,sp,debiased,alpha,threshold,BB){
  
  source("SDR_selection_function.R",local=TRUE)
  
  set.seed(1234567)
  
  
  result <- foreach(seedNo=1:BB,.combine='rbind',.packages= packs)%dopar%
    SDS(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H,trans,sp,debiased,alpha,threshold,seedNo)
  
  result_mean <- apply(result,2,mean)  
  result_sd <- apply(result,2,sd)   
  
  return(list(result = result,result_mean = result_mean, result_sd = result_sd))
  
} 

Func_slice_HD <- function(X,Y,H,trans,debiased){
  n <- nrow(X)
  p <- ncol(X)

  yless <- dr.slices(Y,H)$slice.indicator 
  
  indexlist <- vector("list",H)  #vector
  for(i in 1:H){
    indexlist[[i]] <- which(yless==i)
  }
  ylabel <- numeric()
  ytemp <- yless
  for(i in 1:n){
    if(var(ytemp)!=0){
      ylabel <- c(ylabel,ytemp[1])
      ytemp <- ytemp[ytemp!=ytemp[1]]
    }
  }
  ylabel <- c(ylabel,ytemp[1])
  
  
  
  if (trans == "indicator"){ydis <- rep(1,length(Y))}   #SIR,1991
  if (trans == "CIRE"){ydis <- Y}                                      #cook and Ni (2006)
  if (trans == "poly"){ydis <- Y^2}                                   #Yin and cook (2002),k-th moment 
  
  B <- matrix(0,p,H)
  B_index <- matrix(0,p,H)
  for (i in 1:H){
    ydis[which(yless!=i)] <- 0
    ydis <- ydis - mean(ydis)
    
    YY <- ydis[indexlist[[i]]]
    XX <- X[indexlist[[i]],]
    
    if (debiased == "OLS"){
      #B[,i] <- solve(var(X))%*%(t(X)%*%ydis)/n
      beta <- as.numeric(lm(ydis~X-1)$coefficients)
      sv <- as.vector(which(beta!=0))
      B[,i] <- beta
      B_index[sv,i] <- 1 
      
    }
    if (debiased == "lasso"){
      
      #lasso and debiased lasso the estimate of beta has zero coefficients
      fit <- glmnet(X, ydis, family = "gaussian",alpha=1)    #Default is standardize=TRUE.
      k <- fit$df
      AIC <- deviance(fit)+2*k
      #BIC <- deviance(fit)+log(length(ydis))*k
      i_min <- which.min(AIC)
      lambda_select <- fit$lambda[i_min]/(n_1)
      
      fit_AIC <- glmnet(X, ydis, family = "gaussian", lambda = lambda_select)
      beta_lasso <- as.numeric(fit_AIC$beta[,1])
      sv <- as.vector(which(beta_lasso!=0))     #location of singal
      
      #upper bound p/3
      k <- ceiling(p/3)   #####
      wv <- which(fit$df==max(fit$df[fit$df<k]))[1]
      sv_1 <- which(fit$beta_lasso[, wv]!=0)
      w1_1 <- fit$beta[, wv]
      
      if (k<length(sv)){
        sv <- sv_1
        beta_lasso <- w1_1
      }
      
      B[,i] <- beta_lasso
      B_index[sv,i] <- 1
      
    }
    if (debiased == "debiased_lasso"){
      
      #lasso and debiased lasso the estimate of beta has zero coefficients
      fit <- glmnet(X, ydis, family = "gaussian",alpha=1)    #Default is standardize=TRUE.
      k <- fit$df
      AIC <- deviance(fit)+2*k
      BIC <- deviance(fit)+log(length(ydis))*k
      i_min <- which.min(AIC)
      lambda_select <- fit$lambda[i_min]
      
      fit_AIC <- glmnet(X, ydis, family = "gaussian", lambda = lambda_select,standardize=FALSE)
      beta_lasso <- as.numeric(fit_AIC$beta[,1])
      
      #debiased lasso method 1: theta_lasso + 1/nMt(X)(Y-Xtheta_lasso)
      M <- debiasingMatrix(X, is_wide = TRUE, length(ydis), rows = 1:p)
      beta_debiasedlasso <- beta_lasso + (1/(length(ydis)))*M%*%t(X)%*%(ydis-X%*%beta_lasso)
      B[,i] <- beta_debiasedlasso
      
      sv<- as.vector(which(beta_debiasedlasso!=0))
      B_index[sv,i] <- 1
      
    }
  }
  
  Bp <- eigen(B%*%t(B)/H)$vectors[,1:H]     #SDR subspace
  
  return(list(B = B, Bp = Bp,B_index = B_index, ylabel = ylabel, ydis = ydis,indexlist = indexlist))
  
}

MyFunc_HD <-function(dat,H,trans,sp){
  tuning <- "CV"
  
  X <- dat$X
  Y <- dat$Y
  n <- nrow(X)
  p <- ncol(X)
  
  #sample splitting
  n_1 <- 1*floor(n/2)
  n_2 <- n-n_1
  Ind_1 <- sample(1:n,n_1,replace = F)
  Ind_2 <- (1:n)[-c(Ind_1)]
  dat_1 <- list(X = X[Ind_1,], Y = Y[Ind_1])
  dat_2 <- list(X = X[Ind_2,], Y = Y[Ind_2])
  
  #process in dat1
  X_1 <- dat_1$X
  Y_1 <- dat_1$Y
  
  yless_1 <- dr.slices(Y_1,H)$slice.indicator 
  
  B_1 <- matrix(0,p,H)
  B_1_index <- matrix(0,p,H)
  indexlist_1 <- vector("list",H)  
  for (i in 1:H){
    
    if (trans == "indicator"){ydis_1 <- rep(1,length(Y_1))}   #SIR,1991
    if (trans == "CIRE"){ydis_1 <- Y_1}      #cook and Ni (2006)
    if (trans == "poly"){ydis_1 <- Y_1^2}   #Yin and cook (2002),k-th moment 
    
    ydis_1[which(yless_1!=i)] <- 0
    ydis_1 <- ydis_1 - mean(ydis_1)
    
    indexlist_1[[i]] <- which(yless_1==i)
    
    YY_1 <- ydis_1[indexlist_1[[i]]]
    XX_1 <- X_1[indexlist_1[[i]],]
    
    fit <- LASSO(X_1,ydis_1,method= tuning)
    beta_lasso <- fit$beta_lasso
    sv <- fit$signal_est     #location of singal
    
    p2 <- floor((n_2/H))
    if (length(sv)>p2){
      sv <- (order(abs(beta_lasso),decreasing = T)[1:p2])
      beta_lasso[-sv] <- 0
    }
    
    B_1[,i] <- as.numeric(beta_lasso)
    B_1_index[sv,i] <- 1
    
  }
  
  #process in dat2 with OLS refined by dat1
  X_2 <- dat_2$X
  Y_2 <- dat_2$Y
  
  yless_2 <- dr.slices(Y_2,H)$slice.indicator 
  
  B_2 <- matrix(0,p,H)
  indexlist_2 <- vector("list",H)  
  for (i in 1:H){
    
    if (trans == "indicator"){ydis_2 <- rep(1,length(Y_2))}   #SIR,1991
    if (trans == "CIRE"){ydis_2 <- Y_2}                                      #cook and Ni (2006)
    if (trans == "poly"){ydis_2 <- Y_2^2}                                   #Yin and cook (2002),k-th moment 
    
    ydis_2[which(yless_2!=i)] <- 0
    ydis_2 <- ydis_2 - mean(ydis_2)
    
    indexlist_2[[i]] <- which(yless_2==i)
    
    YY_2 <- ydis_2[indexlist_2[[i]]]
    XX_2 <- X_2[indexlist_2[[i]],]
    
    temp <- which(B_1_index[,i]==1)
    if (length(temp)>0&&length(temp)<=p2){
      beta <- as.numeric(lm(ydis_2~X_2[,temp]-1)$coefficients)
      B_2[temp,i] <- beta
      
      sigma <- rep(1, p)
      DIAG_1 <- diag(solve( t(X_1[,temp])%*%X_1[,temp])*var(ydis_1))
      DIAG_2 <- diag(solve( t(X_2[,temp])%*%X_2[,temp])*var(ydis_2))
      
      #sigma[temp] <- sqrt(DIAG_2*DIAG_1)
      sigma[temp] <-  DIAG_2
      B_2[,i] <- B_2[,i]/sigma
      
    }else{
      B_2[,i] <- rep(0,p)
    }
  } 
  
  #ranking statistics for variable selection
  if(sp=="full")   {
    
    yless <- dr.slices(Y,H)$slice.indicator 
    
    B_full <- matrix(0,p,H)
    B_index <- matrix(0,p,H)
    indexlist <- vector("list",H)  
    for (i in 1:H){
      
      if (trans == "indicator"){ydis <- rep(1,length(Y))}   #SIR,1991
      if (trans == "CIRE"){ydis <- Y}      #cook and Ni (2006)
      if (trans == "poly"){ydis <- Y^2}   #Yin and cook (2002),k-th moment 
      
      ydis[which(yless!=i)] <- 0
      ydis <- ydis - mean(ydis)
      
      indexlist[[i]] <- which(yless==i)
      
      YY <- ydis[indexlist[[i]]]
      XX <- X[indexlist[[i]],]
      
      fit <- LASSO(X,ydis,method= tuning)
      beta_lasso <- fit$beta_lasso
      sv <- fit$signal_est     #location of singal
      
      p2 <- floor((n_2/H))
      if (length(sv)>p2){
        sv <- (order(abs(beta_lasso),decreasing = T)[1:p2])
        beta_lasso[-sv] <- 0
      }
      
      B_full[,i] <- as.numeric(beta_lasso)
      B_index[sv,i] <- 1
    }
    
    
    Wj <- sapply(1:p,function(j){B_full[j,]%*%B_full[j,]})
    Td <- sapply(1:H,function(x){B_full[,x]%*%B_full[,x]})
  }
  if(sp=="split2") {
    Wj <- sapply(1:p,function(j){B_1[j,]%*%B_2[j,]})
    Td <- sapply(1:H,function(x){B_1[,x]%*%B_2[,x]})
  }
  if(sp=="half")   {
    Wj <- sapply(1:p,function(j){B_1[j,]%*%B_1[j,]})
    Td <- sapply(1:H,function(x){B_1[,x]%*%B_1[,x]})
  }
  
  return(list(Wj = Wj,Td = Td))
}
SDS_HD<- function(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H,trans,sp,alpha,threshold,seedNo){
  set.seed(123+1000*seedNo)
  
  tic <- proc.time()
  dat <- Gendata(n, p, d, m, mu, struc,rho, eps, Xd, Sc)   
  tru1 <- which(dat$ss==1)
  
  temp <- MyFunc_HD(dat,H,trans,sp)              
  Wj <-temp$Wj
  
  t <- sort(abs(Wj))   
  
  
  ###FDR 
  
  if (threshold == "-"){
    Ta <- sapply(t,function(x){(0+length(Wj[Wj<=(-x)]))/max(1,length(Wj[Wj>=x]))})
  }  
  if (threshold == "+"){
    Ta <- sapply(t,function(x){(1+length(Wj[Wj<=(-x)]))/max(1,length(Wj[Wj>=x]))}) 
  }  #SDA+
  
  best_L <- min(t[which(Ta<=alpha)])
  
  det1 <- which(Wj>=best_L)  
  FDR <- length(setdiff(det1,tru1))/max(1,length(det1))    
  AP <- length(intersect(det1,tru1))/max(1,length(tru1))  
  NFD <- length(setdiff(det1,tru1))                        
  ND  <- length(det1)                                      
  NN <- length(which(Wj<0))
  
  if (all(tru1%in%det1)==TRUE){Pa <- 1}   #selection consistency 
  if (all(tru1%in%det1)!=TRUE){Pa <- 0}
  
  # mFDR <- NFD/max(1,ND)            #marginal FDR
  # mAP <- (ND-NFD)/max(1,m)         $marginal AP
  
  tic2 <- proc.time()-tic
  time <- tic2[3]
  
  return( c(FDR,AP,NFD,ND,NN,Pa,time))
}
Process_SDS_HD <- function(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H,trans,sp,alpha,threshold,BB){
  
  source("SDR_selection_function.R",local=TRUE)
  
  set.seed(1234567)
  
  
  result <- foreach(seedNo=1:BB,.combine='rbind',.packages= packs)%dopar%
    SDS_HD(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H,trans,sp,alpha,threshold,seedNo)
  
  result_mean <- apply(result,2,mean) 
  result_sd <- apply(result,2,sd)   
  
  return(list(result = result,result_mean = result_mean, result_sd = result_sd))
  
} 

################compared methods#######################################################
#Candes et al (2018), model-X knockoff
MX_knockoff <- function(n,p,d,m,mu,struc,rho,eps,Xd,Sc,alpha,threshold,seedNo){
  set.seed(123+1000*seedNo)
  
  tic <- proc.time()
  dat <- Gendata(n, p, d, m, mu, struc,rho, eps,Xd,Sc)  
  tru1 <- which(dat$ss==1)
  
  X <- dat$X
  Y <- dat$Y
  
  Xtilde <- create.second_order(X, method='equi')  #model-X knockoff,estimated Sigma
 
  #different antisymmetric function
  Wj <- stat.lasso_lambdasmax(X,Xtilde,Y)     #lasso signed maximum,LSM, linear models
  ##Wj <-  stat.lasso_coefdiff(X,Xtilde,Y)     #lasso coefficient difference(LCD), power ??LSM??
  
  ##Wj <- stat.glmnet_lambdasmax(X,Xtilde,Y)  #select lambda by CV with glmnet,the same result with lasso
  #Wj <- stat.random_forest(X,Xtilde,Y)
  t <- sort(abs(Wj))
  
  ###FDR 
  if (threshold == "-"){
    Ta <- sapply(t,function(x){(0+length(Wj[Wj<=(-x)]))/max(1,length(Wj[Wj>=x]))}) #SDA,FDR??Щ???Ʋ?ס?ʸĽ?ΪSDA+
  }
  if (threshold == "+"){
    Ta <- sapply(t,function(x){(1+length(Wj[Wj<=(-x)]))/max(1,length(Wj[Wj>=x]))}) #SDA+,sizeƫ????,??p?ܴ?ʱ??Ч????SDA?????
  }
  
  best_L <- min(t[which(Ta<=alpha)])
  det1 <- which(Wj>=best_L)  
  FDR <- length(setdiff(det1,tru1))/max(1,length(det1))   
  AP <- length(intersect(det1,tru1))/max(1,length(tru1))  
  NFD <- length(setdiff(det1,tru1))                        
  ND  <- length(det1)                                      
  NN <- length(which(Wj<0))
  
  if (all(tru1%in%det1)==TRUE){Pa <- 1}   #selection consistency 
  if (all(tru1%in%det1)!=TRUE){Pa <- 0}
  
  # mFDR <- NFD/max(1,ND)            #marginal FDR
  # mAP <- (ND-NFD)/max(1,m)         $marginal AP
  
  tic2 <- proc.time()-tic
  time <- tic2[3]
  
  return( c(FDR,AP,NFD,ND,NN,Pa,time))
  
}
Process_MX_knockoff <- function(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H,alpha,threshold,BB){
  
  source("SDR_selection_function.R",local=TRUE)
  
  set.seed(1234567)
  
  result <- foreach(seedNo=1:BB,.combine='rbind',.packages= packs)%dopar%
    MX_knockoff(n,p,d,m,mu,struc,rho,eps,Xd,Sc,alpha,threshold,seedNo)
  
  result_mean <- apply(result,2,mean)  
  result_sd <- apply(result,2,sd)   
  
  return(list(result = result,result_mean = result_mean, result_sd = result_sd))
  
} 

#cook(2004)
#marginal cooridinate hypotheses test or Yu (2016), marginal SIR,p<n
MCH_SIR <- function(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H,alpha,seedNo){
  set.seed(123+1000*seedNo)
  
  tic <- proc.time()
  
  dat <- Gendata(n, p, d, m, mu, struc,rho, eps,Xd,Sc)
  tru1 <- which(dat$ss==1)
  
  X <- dat$X
  Y <- dat$Y
  fit_sir <- dr(Y~., method = "sir",nslices = H,numdir = H, data = as.data.frame(cbind(X))) #numdirĬ????4
  #dr.coordinate.test(fit_sir,~-X[1:m])
  MCH_sir <- drop1(fit_sir,update = FALSE)

  pvalue <- MCH_sir$P.value
  
  orderp <- order(pvalue)
  sortp <- sort(pvalue)/(1:p)
  k <- max(which(sortp<=(alpha/p)))    #Bonferroni
  detBH <- NULL
  if(!is.infinite(k)){detBH = orderp[1:k]}             
  
  FDR <- length(setdiff(detBH,tru1))/max(1,length(detBH))     
  AP <- length(intersect(detBH,tru1))/max(1,length(tru1))     
  NFD <- length(setdiff(detBH,tru1))                           
  ND <- length(detBH)                                         
  NN <- length(which(pvalue< alpha))                         
  
  if (all(tru1%in%detBH)==TRUE){Pa <- 1}   #selection consistency 
  if (all(tru1%in%detBH)!=TRUE){Pa <- 0}
  
  tic2 <- proc.time() -tic
  time <- tic2[3]
  
  return( c(FDR,AP,NFD,ND,NN,Pa,time))             
  
}
Process_MCH_SIR <- function(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H,alpha,BB){
  
  source("SDR_selection_function.R",local=TRUE)
  
  set.seed(1234567)
  
  result <- foreach(seedNo=1:BB,.combine='rbind',.packages= packs,.errorhandling="remove")%dopar%
    MCH_SIR(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H,alpha,seedNo)
  
  result_mean <- apply(result,2,mean) 
  result_sd <- apply(result,2,sd)   
  
  return(list(result = result,result_mean = result_mean, result_sd = result_sd))
  
} 

#marginal SIR with screening under high dimension p>n
MCH_SIR_HD <- function(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H,alpha,seedNo){
  set.seed(123+1000*seedNo)
  
  tic <- proc.time()
  
  tuning <- "CV"
  
  dat <- Gendata(n, p, d, m, mu, struc,rho, eps,Xd,Sc)
  tru1 <- which(dat$ss==1)
  
  X <- dat$X
  Y <- dat$Y
  
  #sample splitting
  n_1 <- 1*floor(n/2)
  n_2 <- n-n_1
  Ind_1 <- sample(1:n,n_1,replace = F)
  Ind_2 <- (1:n)[-c(Ind_1)]
  dat_1 <- list(X = X[Ind_1,], Y = Y[Ind_1])
  dat_2 <- list(X = X[Ind_2,], Y = Y[Ind_2])
  
  #process in dat1
  X_1 <- dat_1$X
  Y_1 <- dat_1$Y
  
  fit <- LASSO(X_1,Y_1,method= tuning)
  beta_lasso <- fit$beta_lasso
  sv <- fit$signal_est     #location of signal
  p2 <- fit$signal_num
  m2 <- floor(n_2/H)
  if (p2>m2){
    sv <- (order(abs(beta_lasso),decreasing = T)[1:m2])
    beta_lasso[-sv] <- 0
  }
  
  B_1 <- as.numeric(beta_lasso)
  sv <- sort(sv)
  
  #process in dat2 with selected covariates
  X_2 <- dat_2$X
  Y_2 <- dat_2$Y
  dat2 <- data.frame(X_2 = X_2[,sv],Y_2 = Y_2)
  
  fit_sir <- dr(Y_2~., method = "sir",nslices = H,numdir = H, data = dat2) #numdir default=4
  MCH_sir <- drop1(fit_sir,update = FALSE)  #update=F equals to drop1.dr
  pvalue <- MCH_sir$P.value
  
  orderp <- order(pvalue)
  t2 <- min(p2,m2)
  sortp <- sort(pvalue)/(1:t2)
  k <- max(which(sortp<=(alpha/t2)))    #Bonferroni
  detBH <- NULL
  if(!is.infinite(k)){detBH = orderp[1:k]}               
  
  FDR <- length(setdiff(detBH,tru1))/max(1,length(detBH))      
  AP <- length(intersect(detBH,tru1))/max(1,length(tru1))      
  NFD <- length(setdiff(detBH,tru1))                            
  ND <- length(detBH)                                          
  NN <- length(which(pvalue< alpha))                          
  
  if (all(tru1%in%detBH)==TRUE){Pa <- 1}   #selection consistency 
  if (all(tru1%in%detBH)!=TRUE){Pa <- 0}
  
  tic2 <- proc.time() -tic
  time <- tic2[3]
  
  return( c(FDR,AP,NFD,ND,NN,Pa,time))             
  
}
Process_MCH_SIR_HD <- function(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H,alpha,BB){
  
  source("SDR_selection_function.R",local=TRUE)
  
  set.seed(1234567)
  
  result <- foreach(seedNo=1:BB,.combine='rbind',.packages= packs,.errorhandling="remove")%dopar%
    MCH_SIR_HD(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H,alpha,seedNo)
  
  result_mean <- apply(result,2,mean) 
  result_sd <- apply(result,2,sd)   
  
  return(list(result = result,result_mean = result_mean, result_sd = result_sd))
  
} 

#Yu,Dong and Shao(2016) marginal SIR(marginal independence SIR, let Sigma=Ip)
#marginal independence SIR, I-MSIR in Yu,Dong and Shao(2016)
IM_SIR <- function(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H,c,alpha,seedNo){
  set.seed(123+1000*seedNo)
  
  tic <- proc.time()
  
  dat <- Gendata(n, p, d, m, mu, struc,rho, eps,Xd,Sc)
  tru1 <- which(dat$ss==1)
  
  X <- dat$X
  Y <- dat$Y
  
  sir <- function(d,x,y,H){
    p <- ncol(x)
    n <- nrow(x)
    
    N <- var(x)
    xst <- x
    ydis <- dr.slices(y,H)$slice.indicator
    
    yless <- ydis
    ylabel <- numeric()
    for(i in 1:n) {if(var(yless)!=0) {
      ylabel <- c(ylabel,yless[1])
      yless <- yless[yless!=yless[1]]}}
    ylabel <- c(ylabel,yless[1])
    prob <- numeric()
    exy <- numeric()
    for(i in 1:H) prob <- c(prob,length(ydis[ydis==ylabel[i]])/n) #the proportion for each slice
    for(i in 1:H) exy <- rbind(exy,apply(xst[ydis==ylabel[i],],2,mean)) #the z mean for each slice
    sirmat <- t(exy)%*%diag(prob)%*%exy
    M <- sirmat
    sig <- matpower(N,1/2)
    B <- (eigen(sirmat)$vectors[,1:d])  
    B <- as.matrix(B)
    
    return(list(B = B, M = sirmat, N = N, mat = sirmat))
  }
  M <- sir(d,X,Y,H)$M
  Wj <- diag(M)       #marginal independence SIR, I-MSIR
  
  
  best_L <- floor(c*n/log(n))
  det1 <- order(Wj,decreasing = TRUE)[1:best_L]  #select the top L large as signals

  
  FDR <-length(setdiff(det1,tru1))/max(1,length(det1))    #under H1, false discovery
  AP <- length(intersect(det1,tru1))/max(1,length(tru1))   #under H1, true discovery
  NFD <- length(setdiff(det1,tru1))                        #number of false discovery
  ND  <- length(det1)                                      #number of discovery
  NN <- length(which(Wj<0))
  
  if (all(tru1%in%det1)==TRUE){Pa <- 1}   #selection consistency 
  if (all(tru1%in%det1)!=TRUE){Pa <- 0}
  
  tic2 <- proc.time() -tic
  time <- tic2[3]
  return( c(FDR,AP,NFD,ND,NN,Pa,time))             
  
}
Process_IM_SIR <- function(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H,c,alpha,BB){
  
  source("SDR_selection_function.R",local=TRUE)
  
  set.seed(1234567)
  
  result <- foreach(seedNo=1:BB,.combine='rbind',.packages= packs)%dopar%
    IM_SIR(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H,c,alpha,seedNo)
  
  result_mean <- apply(result,2,mean)  
  result_sd <- apply(result,2,sd) 
  
  return(list(result = result,result_mean = result_mean, result_sd = result_sd))
  
} 
