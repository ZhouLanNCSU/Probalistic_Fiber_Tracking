################################################
##########Simulations Thesis Use################
################################################



library("foreach")
library("doParallel")
registerDoParallel(cores=1)

Fiber_Tracking_05<-foreach(it=c(1:200),.errorhandling=c('pass') )%dopar%{
  setwd("DATA_SCRIPTS/")
  library("bayess")
  library(rlang)
  library('igraph')
  library("bayess")
  library('rlist')
  library('MASS')
  library('matrixcalc')
  library('MCMCpack')
  library('inline')
  library('RcppArmadillo')
  library('Rcpp')
  library("rstan")
  library("Boom")
  #library("profvis")
  library("mclust")
  library("plot3D")
  library("geoR")
  library("mvtnorm")
#setwd("/Volumes/Transcend/DiST")
source("dwi_fit.R")
source("dwi_basic.R")

#######################
#### Simulate data ####
#######################
source("sim-curve-new.R")
load("pre_estimation.Rdata") 

Locations_temp=unique(pre$loc)
X_index=unique(pre$loc[,1])
Y_index=unique(pre$loc[,2])[1:7]
Z_index=unique(pre$loc[,3])[c(2,3)]


index=which(Locations_temp[,1]%in%X_index &
              Locations_temp[,2]%in%Y_index & 
              Locations_temp[,3]%in%Z_index )
Locations=Locations_temp[index,]


m=10
sigma=0.1
n=16*7*2
b=1
M=33
g=t(grad.mat)
adj=get.adjacency(make_lattice(dimvector=c(16,7,2)))
adj[lower.tri(adj, diag = FALSE)]=0

log_SM=log(dwi.obs[,index])+matrix(rnorm(n*M,0,0.5),33,224)
log_S0=rep(log(S0const),length(index))




f_opt<-function(A){
  
  A_mat=diag(3)
  A_mat[1,1]=A[1]
  A_mat[2,2]=A[2]
  A_mat[3,3]=A[3]
  A_mat[2,1]=A_mat[1,2]=A[4]
  A_mat[3,1]=A_mat[1,3]=A[5]
  A_mat[3,2]=A_mat[2,3]=A[6]
  
  
  return(
    sum((sapply(1:M, function(mm) log_SM[mm,v]-log_S0[v]+g[,mm]%*%A_mat%*% g[,mm]))^2)
  )
}



A_list=list()

for(v in 1:n){
  
  A=optim(c(1,1,1,0,0,0),f_opt)$par
  A_mat=diag(3)
  A_mat[1,1]=A[1]
  A_mat[2,2]=A[2]
  A_mat[3,3]=A[3]
  A_mat[2,1]=A_mat[1,2]=A[4]
  A_mat[3,1]=A_mat[1,3]=A[5]
  A_mat[3,2]=A_mat[2,3]=A[6]
  if(is.positive.definite(A_mat)){
    A_list[[v]]=A_mat
  }else{
    A_list[[v]]=rWishart(3,diag(3)/3)
  }
}

obs=0
FAmap<-function(lam){
  lb<-mean(lam)
  fa<-sqrt(3*((lam[1]-lb)**2+(lam[2]-lb)**2+(lam[3]-lb)**2)/sum(2*lam**2))
}

AC_list=A_list



J=1
Zj=sapply(1:M,function(j)  1 )
load("True_Tensors.Rdata")
A_list=lapply(1:J, function(j) A_list)
mean_A_list=lapply(1:J,function(j) lapply(1:n, function(v)  Reduce("+",A_list[[j]][which(adj[v,]==1)])/sum(adj[v,]==1)))

AC=0
OBS=0



m=3

iters=5000

Bookeeping=list()
for (it in 1:iters){
  
  #######################################
  ##### Update DT A_v ###################
  #######################################
  for(j in 1:J){
    
    for (v in 1:n){
      OBS=1+OBS
      
      q=100
      
      A_can=rWishart(q,A_list[[j]][[v]]/q)
      A_old=A_list[[j]][[v]]
      
      can=0
      old=0
      if(sum(Zj==j)!=0){
        ####Component 1
        mean_can=log_S0[v]-b*sapply(which(Zj==j),function(mm)  g[,mm]%*%A_can%*% g[,mm] )
        mean_old=log_S0[v]-b*sapply(which(Zj==j),function(mm)  g[,mm]%*%A_old%*% g[,mm] )
        
        can=can+dmvnorm(log_SM[Zj==j,v], mean = mean_can, sigma = diag(sum(Zj==j))*sigma, log = TRUE)
        old=old+dmvnorm(log_SM[Zj==j,v], mean = mean_old, sigma = diag(sum(Zj==j))*sigma, log = TRUE)
        
      }
      
      ####Component 2
      try(can<-can+dWishart(A_can,mean_A_list[[j]][[v]]/m,m,TRUE),TRUE)
      try(old<-old+dWishart(A_old,mean_A_list[[j]][[v]]/m,m,TRUE),TRUE)
      
      ####Component 3
      u_index=which(adj[,v]==1)
      
      mean_A_list_can=lapply(u_index, function(u)  (ifelse(is.null(Reduce("+", A_list[[j]][setdiff(which(adj[u,]==1),v)])),0,   Reduce("+", A_list[[j]][setdiff(which(adj[u,]==1),v)]))
                                                    +  A_can)/sum(adj[u,]==1))
      mean_A_list_old=mean_A_list[[j]][u_index]
      
      try(can<-can+sum(sapply(1:length(u_index), function(uu)  dWishart(A_list[[j]][[u_index[uu]]],mean_A_list_can[[uu]]/m,m,TRUE))),TRUE)
      try(old<-old+sum(sapply(1:length(u_index), function(uu)  dWishart(A_list[[j]][[u_index[uu]]],mean_A_list_old[[uu]]/m,m,TRUE))),TRUE)
      
      ####Component 3
      can=can+dWishart(A_old,A_can/m,m,TRUE)
      old=old+dWishart(A_can,A_old/m,m,TRUE)
      
      
      
      prob=min(exp(can-old),1)
      Accept=sample(c(0,1),prob=c(1-prob,prob),size = 1)
      
      if(Accept==1){
        AC=1+AC
        #cat("Acceptance Rate",AC/OBS,"\n")
        A_list[[j]][[v]]<-A_can
        for(uu in 1:length(u_index)){
          if(length(u_index)!=0){
            mean_A_list[[j]][[u_index[[uu]]]]<-mean_A_list_can[[uu]]
          }
        }
      }
      
      
      
    }
    
    
  }
  
  
  
  #######################################
  ##### Update DOF ######################
  #######################################
  m_can=exp(rnorm(1,log(m),1))
  
  if(m_can<3){m_can=3}
  if(m_can>50){m_can=50}
  
  can=sum(sapply(1:J,function(j) sapply(1:n,function(i) ifelse(is.matrix(mean_A_list[[j]][[i]]), dWishart(A_list[[j]][[i]], mean_A_list[[j]][[i]]/m_can,m_can,TRUE),0 ))))
  old=sum(sapply(1:J,function(j) sapply(1:n,function(i) ifelse(is.matrix(mean_A_list[[j]][[i]]), dWishart(A_list[[j]][[i]], mean_A_list[[j]][[i]]/m,m_can,TRUE),0 ))))
  
  prob=min(exp(can-old),1)
  Accept=sample(c(0,1),prob=c(1-prob,prob),size = 1)
  
  if(Accept==1){
    
    m=m_can
    
  }
  
  
  ##############################################:
  #####          VARIANCE (Gibbs)        #######:
  ##############################################:
  
  R      <- as.numeric(log_SM)-log_S0%x%rep(1,M)+as.numeric(
    sapply(1:n,function(v) b*sapply(1:M,function(mm)  g[,mm]%*%A_list[[Zj[mm]]][[v]]%*% g[,mm] )))
  a      <- n*M/2+0.01
  bb      <- t(R)%*%R/2+0.01
  sigma <- 1/rgamma(1,a,bb) 
  
  ##############################################:
  #####          Probabilities (Gibbs)   #######:
  ##############################################:
  
  probability=rdirichlet(sapply(1:J, function(j) 50+sum(Zj==j))  )
  
  ##############################################:
  #####          Latent Variable ZJ   #######:
  ##############################################:
  
  # for(mm in 1:M){
  #   log_p=sapply(1:J,function(j) log(probability[j])+sum(sapply(1:n,function(v) dnorm(log_SM[mm,v],log_S0[v]-b*g[,mm]%*%A_list[[j]][[v]]%*%g[,mm],sqrt(sigma),log=TRUE) )))
  #   
  #   prob=1/(1+exp(log_p[2]-log_p[1]))
  #   Zj[mm]=sample(c(1,2),prob=c(prob,1-prob),size=1)
  #   
  # }
  # print(Zj)
  
  
  Bookeeping[[it]]=A_list[[1]]
  
}

AB_list=Bookeeping[[it-1]]

result=list(SpDiST=AB_list,LS=AC_list)

}


save.image("Fiber_Tracking_01.Rdata")