rm(list = ls())
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

#setwd('C:/Users/zl527/Dropbox/LAN_ALL/VIP/DTI-ILD')
load("A_list.Rdata")
x<-scan("directions15.txt")
k<-length(x)/3
#############################file
nx<-50
ny<-20
z<-scan("part50x20.csv",sep=",")
#####################################
y<-array(z,c(ny,16,nx))
par(mfrow=c(4,4))
for (i in 1:16){
  image(t(y[,i,]))
}
n<-length(z)/(k+1)
g<-t(matrix(x,k,3)) ####gvalue
Sdata<-matrix(0,k+1,n) #### S obs
iobs<-0
for (ix in 1:nx){
  for (iy in 1:ny){
    iobs<-iobs+1
    Sdata[,iobs]<-y[iy,,ix]
  }
}
b=1 #### bvalue


M=15
adj=get.adjacency(make_lattice(dimvector=c(20,50)))
adj[lower.tri(adj, diag = FALSE)]=0
log_SM=log(Sdata)[2:(M+1),]
log_S0=log(Sdata)[1,]
m=10
sigma=0.1
n=1000


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



J=1
Zj=sapply(1:M,function(j)  1 )
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
      
      q=1000
      
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
  
  
  
  # #######################################
  # ##### Update DOF ######################
  # #######################################
  # m_can=exp(rnorm(1,log(m),1))
  # 
  # if(m_can<3){m_can=3}
  # if(m_can>50){m_can=50}
  # 
  # can=sum(sapply(1:J,function(j) sapply(1:n,function(i) ifelse(is.matrix(mean_A_list[[j]][[i]]), dWishart(A_list[[j]][[i]], mean_A_list[[j]][[i]]/m_can,m_can,TRUE),0 ))))
  # old=sum(sapply(1:J,function(j) sapply(1:n,function(i) ifelse(is.matrix(mean_A_list[[j]][[i]]), dWishart(A_list[[j]][[i]], mean_A_list[[j]][[i]]/m,m_can,TRUE),0 ))))
  # 
  # prob=min(exp(can-old),1)
  # Accept=sample(c(0,1),prob=c(1-prob,prob),size = 1)
  # 
  # if(Accept==1){
  #   
  #   m=m_can
  #   
  # }
  
  
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

# AC_list=Bookeeping[[1]]
# obs=0
# for (ix in 1:nx){
#   for (iy in 1:ny){
#     obs=obs+1
#     #AC_list[[obs]]=aa[[obs]]
#     wire3d(ellipse3d(AC_list[[obs]]/(sum(diag(AC_list[[obs]])))*ifelse(FAmap(eigen(AC_list[[obs]])$values)==0,0.01,  FAmap(eigen(AC_list[[obs]])$values))
#                      , centre=c(ix,iy,0),level=0.2),box=FALSE,axes=FALSE,col="yellow")
#     
#   }
# }


ttt<-function(Mean){
  
  if(Mean[1]<0){
    
    Mean=Mean
  }else{Mean=-Mean}
 
  return(Mean)
}


xlim <- c(0 , 51)
ylim <- c(0, 21)
plot(0, type = "n", xlim = xlim, ylim = ylim, 
     main = "Arrows,  type = 'curved'")

for (it in 1:5){
  print(it)
  if(it %% 1==0){
    obs=0
    for (ix in 1:nx){
      for (iy in 1:ny){
        obs=obs+1
        eee=eigen(Bookeeping[[it]][[obs]][1:2,1:2])
        #if( FAmap(eee$values)>0){
          Mean=eee$vectors[,1]
        
          Arrows(ix,iy,ix+Mean[1]/2,iy+Mean[2]/2,arr.width=0.03,arr.length = 0.001) 
          Arrows(ix,iy,ix-Mean[1]/2,iy-Mean[2]/2,arr.width=0.03,arr.length = 0.001) 
        #}
        
    
      }
    }
    
    
    
  }
  
  
}



data=cbind(seq(18,28,0.01),sapply(1:1001,function(i) TT_all[[i]]/500))
data=data.frame(data)
colnames(data)<-c("Degrees","Probability")
ggplot(data, aes(Degrees, Probability)) +
  geom_point() +
  geom_smooth(se=FALSE,method = "loess")+
ylab("The Probability of Pattern B")









