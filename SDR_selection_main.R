##############################################################################
###Code for "Model-free controlled variable selection via data splitting"#####
###Yixin Han, Xu Guo, and Changliang Zou (2024) SCIENTIA SINICA Mathematica###
###https://arxiv.org/abs/2210.12382###########################################
###https://www.sciengine.com/SSM/doi/10.1360/SCM-2022-0066####################
##############################################################################

rm(list=ls())
#n         total sample size 
#n1        subdata1 sample size
#n2        subdata2 sample size
#p         variable dimension 
#m         important variable number
#d         SDR structural dimension
#BB        simulation replication times
#rho       covariance matrix parameter
#eps       random noise term
#alpha     significant level
#mu        signal amplitude
#H         slice number
#trans     response transformation function
#debiased  debiased lasso 
#Xd        covariate distribution
#Sc        Model example, Scenarios 


library(MASS)
library(mnormt)
library(dr)
library(MAVE)
library(ggplot2)
library(latex2exp)
library(knockoff)
library(selectiveInference)
library(hdi)
library(CovTools)
library(ranger)
library(ncvreg)


######################################################################################
#global fixed parameters 
BB <- 500
sp <- "split2"
alpha <- 0.2
mu <- 1
trans <- "indicator"

#packages
packs <- c("dr","MASS","knockoff","selectiveInference",
           "mnormt","ranger","ncvreg")

setwd("C:/Users/STAT/hyx/SDR selection/code")
source("SDR_selection_function.R",local=TRUE)

#d <- 1; m <- 10    #  Scenario I
#d <- 2; m <- 10     #  Scenario II
#d <- 3; m <- 10     #  Scenario III

#######slice number vs response transformation in low dimension#############################################################################################################
#result
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
cl <- makeCluster(35)
registerDoParallel(cl)

BB <- 500
Xd <- "norm"
#Sc <- 1; d <-1
#Sc <- 2; d <-2
Sc <- 3; d <-3 
m <- 10
p <- 20
n <- 500        #n/2 > Hp
rho <- 0.5         #for simplicity
struc <- "AR1"
eps <- "norm"
mu <- 1
sp <- "split2"
debiased <- "OLS"

H <-c(3,4,5,6,7,8,10)     


result_temp <-  matrix(0,length(H),7)
colnames(result_temp) <- c("FDR","AP","NFD","ND","NN","Pa","time")
rownames(result_temp) <- H
result_H <- list(result_indicator = result_temp,
                 result_CIRE = result_temp,
                 result_poly = result_temp)
result_H_sd <- result_H

for(k in 1:length(H)){
  tic=proc.time()
  
  temp_indicator <- Process_SDS(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H[k],trans="indicator",sp,debiased,alpha,threshold="-",BB)
  result_H$result_indicator[k,] <- temp_indicator$result_mean
  result_H_sd$result_indicator[k,] <- temp_indicator$result_sd
  
  temp_CIRE <- Process_SDS(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H[k],trans="CIRE",sp,debiased,alpha,threshold="-",BB)
  result_H$result_CIRE[k,] <- temp_CIRE$result_mean
  result_H_sd$result_CIRE[k,] <- temp_CIRE$result_sd
  
  temp_poly <- Process_SDS(n,p,d,m,mu,struc,rho,eps,Xd,Sc,H[k],trans="poly",sp,debiased,alpha,threshold="-",BB)
  result_H$result_poly[k,] <- temp_poly$result_mean
  result_H_sd$result_poly[k,] <- temp_poly$result_sd
  
  tic2=proc.time()-tic
  
  print(paste("H",H[k],"-p",p,"-m",m,"-time:",tic2[3],sep=""))
  
}
stopCluster(cl)

result_H
write.csv(result_H,file='result_H.csv')
save(list=c("result_H"),file='result_H.Rdata')


dat_result_indicator <- data.frame(power = c(result_H$result_indicator[,1],result_H$result_indicator[,2]),
                                   H = rep(H,2),
                                   Type = factor(rep(c("FDR","AP"),each=length(H)),levels=c("FDR","AP")),
                                   method = rep("SDS",2*length(H)),
                                   transform = rep("indicator",2*length(H)))
dat_result_CIRE <- data.frame(power = c(result_H$result_CIRE[,1],result_H$result_CIRE[,2]),
                              H = rep(H,2),
                              Type = factor(rep(c("FDR","AP"),each=length(H)),levels=c("FDR","AP")),
                              method = rep("SDS",2*length(H)),
                              transform = rep("CIRE",2*length(H)))
dat_result_poly <- data.frame(power = c(result_H$result_poly[,1],result_H$result_poly[,2]),
                              H = rep(H,2),
                              Type = factor(rep(c("FDR","AP"),each=length(H)),levels=c("FDR","AP")),
                              method = rep("SDS",2*length(H)),
                              transform = rep("poly",2*length(H)))

dat_result_H <- rbind(dat_result_indicator,dat_result_CIRE,dat_result_poly)
baseline <- data.frame(Type = "FDR",value = alpha)


library(ggplot2)
library(grid)
library(latex2exp)
require("RColorBrewer")
wd <- 1 
fig_H <- ggplot(dat_result_H, aes(x=H, y=power, color = transform,linetype= transform, shape = transform)) +
  geom_point(aes(shape = transform),size = 2*wd)+
  geom_line(aes(linetype = transform),size = wd)+
  theme(legend.position = "top")+
  xlab(TeX('$H$'))+
  ylab("power")+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(data=baseline,aes(yintercept=value),colour='darkgrey',linetype =1,size=1)

fig_H + facet_grid(Type~., scale="free_y")


####################signal amplitude mu######################################################################
#result: different methods vs mu
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
cl=makeCluster(25)
registerDoParallel(cl)

Xd <- "norm"
Sc <- 3; d <-3
m <- 10
p <- 20
n <- 1000
H <- 4          
rho <- 0.7
struc <- "AR1"
eps <- "norm"
debiased <- "OLS"

mu <- seq(0.5,1.5,0.25)

result_temp <-  matrix(0,length(mu),7)
colnames(result_temp) <- c("FDR","AP","NFD","ND","NN","Pa","time")
rownames(result_temp) <- mu
result_mu <- list(result_indicator = result_temp,result_CIRE = result_temp,result_poly = result_temp,
                  result_MCH_SIR = result_temp,
                  result_debiasedlasso = result_temp, 
                  result_MX_knockoff = result_temp)
result_mu_sd <- result_mu

for(k in 1:length(mu)){
  
  tic <- proc.time()
  temp_indicator <- Process_SDS(n,p,d,m,mu[k],struc,rho,eps,Xd,Sc,H,trans="indicator",sp,debiased,threshold = "-",alpha,BB)
  result_mu$result_indicator[k,] <- temp_indicator$result_mean
  result_mu_sd$result_indicator[k,] <- temp_indicator$result_sd
  
  temp_CIRE <- Process_SDS(n,p,d,m,mu[k],struc,rho,eps,Xd,Sc,H,trans="CIRE",sp,debiased,threshold = "-",alpha,BB)
  result_mu$result_CIRE[k,] <- temp_CIRE$result_mean
  result_mu_sd$result_CIRE[k,] <- temp_CIRE$result_sd
  
  temp_poly <- Process_SDS(n,p,d,m,mu[k],struc,rho,eps,Xd,Sc,H,trans="poly",sp,debiased,threshold = "-",alpha,BB)
  result_mu$result_poly[k,] <- temp_poly$result_mean
  result_mu_sd$result_poly[k,] <- temp_poly$result_sd
  
  temp_MCH_SIR <- Process_MCH_SIR(n,p,d,m,mu[k],struc,rho,eps,Xd,Sc,H,alpha,BB)
  result_mu$result_MCH_SIR[k,] <- temp_MCH_SIR$result_mean
  result_mu_sd$result_MCH_SIR[k,] <- temp_MCH_SIR$result_sd
  
  temp_MX_knockoff <- Process_MX_knockoff(n,p,d,m,mu[k],struc,rho,eps,Xd,Sc,H,threshold = "-",alpha,BB)
  result_mu$result_MX_knockoff[k,] <- temp_MX_knockoff$result_mean
  result_mu_sd$result_MX_knockoff[k,] <- temp_MX_knockoff$result_sd
  
  tic2 <- proc.time()-tic
  
  print(paste("mu",mu[k],"-p",p,"-time:",tic2[3],sep=""))
  
}
result_mu

stopCluster(cl)

dat_result_indicator <- data.frame(power = c(result_mu$result_indicator[,1],result_mu$result_indicator[,2]),
                                   mu = rep(mu,2),
                                   Type = factor(rep(c("FDR","AP"),each=length(mu)),levels=c("FDR","AP")),
                                   method = rep("SDS_inducator",2*length(mu)),
                                   transform = rep("indicator",2*length(mu)))

dat_result_CIRE <- data.frame(power = c(result_mu$result_CIRE[,1],result_mu$result_CIRE[,2]),
                              mu = rep(mu,2),
                              Type = factor(rep(c("FDR","AP"),each=length(mu)),levels=c("FDR","AP")),
                              method = rep("SDS_CIRE",2*length(mu)),
                              transform = rep("CIRE",2*length(mu)))

dat_result_poly <- data.frame(power = c(result_mu$result_poly[,1],result_mu$result_poly[,2]),
                              mu = rep(mu,2),
                              Type = factor(rep(c("FDR","AP"),each=length(mu)),levels=c("FDR","AP")),
                              method = rep("SDS_poly",2*length(mu)),
                              transform = rep("poly",2*length(mu)))
dat_result_MCH_SIR <- data.frame(power = c(result_mu$result_MCH_SIR[,1],result_mu$result_MCH_SIR[,2]),
                                 mu = rep(mu,2),
                                 Type = factor(rep(c("FDR","AP"),each=length(mu)),levels=c("FDR","AP")),
                                 method = rep("MCH_SIR",2*length(mu)),
                                 transform = rep("none",2*length(mu)))
dat_result_MX_knockoff <- data.frame(power = c(result_mu$result_MX_knockoff[,1],result_mu$result_MX_knockoff[,2]),
                                     mu = rep(mu,2),
                                     Type = factor(rep(c("FDR","AP"),each=length(mu)),levels=c("FDR","AP")),
                                     method = rep("MX_knockoff",2*length(mu)),
                                     transform = rep("none",2*length(mu)))
dat_result_mu <- rbind(dat_result_indicator,dat_result_CIRE,dat_result_poly,
                       dat_result_MCH_SIR,dat_result_MX_knockoff)
baseline <- data.frame(Type="FDR",value=alpha)


library(ggplot2)
library(grid)
library(latex2exp)
require("RColorBrewer")
wd <- 1 
fig_mu <- ggplot(dat_result_mu, aes(x=mu, y=power, color = method,linetype= method, shape = method)) +
  geom_point(aes(shape = method),size = 2*wd)+
  geom_line(aes(linetype = method),size = wd)+
  theme(legend.position = "top")+
  xlab(TeX('$\\mu$'))+
  ylab("power")+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(data=baseline,aes(yintercept=value),colour='darkgrey',linetype =3,size=1)

fig_mu + facet_grid(Type~., scale="free_y")



########result:for different sample size n#######
#p>n: marginalSIR+BH,MX-knockoff,SDS-Debiased.
#when p is large, use threshold+

library(parallel)
library(foreach)
library(iterators)
library(doParallel)
cl <- makeCluster(39)
registerDoParallel(cl)

source("SDR_selection_function.R",local=TRUE)
source("lasso_inference.R",local=TRUE)


BB <- 100
Xd <- "norm"  
#Sc <- 11;d <- 1 
#Sc <- 22;d <- 2 
Sc <- 33;d <- 3 
H <- 4         #4,5,6,8,10
rho <- 0.5
struc <- "AR1"
eps <- "norm"
mu <- 1
p <- 1000
m <- 10
n <- c(400,500,800)


result_temp <-  matrix(0,length(n),7)
colnames(result_temp) <- c("FDR","AP","NFD","ND","NN","Pa","time")
rownames(result_temp) <- n
result_n <- list(result_SDS_HD_indicator = result_temp,
                 result_SDS_HD_indicator_plus = result_temp,
                 result_debiasedlasso = result_temp,result_debiasedlasso_plus = result_temp,
                 result_IM1_SIR = result_temp,result_IM2_SIR = result_temp,
                 result_MCH_SIR_HD = result_temp,
                 result_MX_knockoff = result_temp,result_MX_knockoff_plus = result_temp)
result_n_sd <- result_n

for(k in 1:length(n)){
  
  tic <- proc.time()
  
  #SDS-high dimensional, plus
  temp_SDS_HD_indicator_plus <- Process_SDS_HD(n[k],p,d,m,mu,struc,rho,eps,Xd,Sc,H,trans="indicator",sp,alpha,threshold="+",BB)
  result_n$result_SDS_HD_indicator_plus[k,] <-  temp_SDS_HD_indicator_plus$result_mean
  result_n_sd$result_SDS_HD_indicator_plus[k,] <- temp_SDS_HD_indicator_plus$result_sd
  
  #SDS-debiased lasso
  temp_debiasedlasso <- Process_SDS(n[k],p,d,m,mu,struc,rho,eps,Xd,Sc,H,trans="indicator",sp,debiased="debiased_lasso",alpha,threshold="-",BB)
  result_n$result_debiasedlasso[k,] <- temp_debiasedlasso$result_mean
  result_n_sd$result_debiasedlasso[k,] <- temp_debiasedlasso$result_sd
  
  temp_debiasedlasso_plus <- Process_SDS(n[k],p,d,m,mu,struc,rho,eps,Xd,Sc,H,trans="indicator",sp,debiased="debiased_lasso",alpha,threshold="+",BB)
  result_n$result_debiasedlasso_plus[k,] <- temp_debiasedlasso_plus$result_mean
  result_n_sd$result_debiasedlasso_plus[k,] <- temp_debiasedlasso_plus$result_sd
  #
  #independent marginal SIR, 2016
  temp_IM1_SIR <- Process_IM_SIR(n[k],p,d,m,mu,struc,rho,eps,Xd,Sc,H,c=0.05,alpha,BB)
  result_n$result_IM1_SIR[k,] <- temp_IM1_SIR$result_mean
  result_n_sd$result_IM1_SIR[k,] <- temp_IM1_SIR$result_sd
  # 
  temp_IM2_SIR <- Process_IM_SIR(n[k],p,d,m,mu,struc,rho,eps,Xd,Sc,H,c=0.5,alpha,BB)
  result_n$result_IM2_SIR[k,] <- temp_IM2_SIR$result_mean
  result_n_sd$result_IM2_SIR[k,] <- temp_IM2_SIR$result_sd
  
  
  #marginal SIR in high dimension: screening by lasso and cook 2004,BH
  temp_MCH_SIR_HD <- Process_MCH_SIR_HD(n[k],p,d,m,mu,struc,rho,eps,Xd,Sc,H,alpha,BB)
  result_n$result_MCH_SIR_HD[k,] <- temp_MCH_SIR_HD$result_mean
  result_n_sd$result_MCH_SIR_HD[k,] <- temp_MCH_SIR_HD$result_sd

  #Model-X knockoff
  temp_MX_knockoff <- Process_MX_knockoff(n[k],p,d,m,mu,struc,rho,eps,Xd,Sc,H,alpha,threshold="-",BB)
  result_n$result_MX_knockoff[k,] <- temp_MX_knockoff$result_mean
  result_n_sd$result_MX_knockoff[k,] <- temp_MX_knockoff$result_sd
  
  temp_MX_knockoff_plus <- Process_MX_knockoff(n[k],p,d,m,mu,struc,rho,eps,Xd,Sc,H,alpha,threshold="+",BB)
  result_n$result_MX_knockoff_plus[k,] <- temp_MX_knockoff_plus$result_mean
  result_n_sd$result_MX_knockoff_plus[k,] <- temp_MX_knockoff_plus$result_sd
  
  
  tic2 <- proc.time()-tic
  
  print(paste("d",d,"-n",n[k],"-time:",tic2[3],sep=""))
  
}
stopCluster(cl)


result_n
save(list=c("result_n"),file='result_n.Rdata')

