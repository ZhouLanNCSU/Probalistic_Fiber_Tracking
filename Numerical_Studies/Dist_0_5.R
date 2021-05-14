
library("foreach")
library("doParallel")
registerDoParallel(cores=16)

result<-foreach(it=c(1:200),.errorhandling=c('pass') )%dopar%{
  setwd("DATA_SCRIPTS/")
# need ncpu
ncpu=1

################
#### Source ####
################
source("dwi_fit.R")
source("dwi_basic.R")

#######################
#### Simulate data ####
#######################
source("sim-curve-new.R")

####################################
#### fit voxel level estimation ####
####################################
result=list()

load("pre_estimation.Rdata") 

Locations_temp=unique(pre$loc)
X_index=unique(pre$loc[,1])
Y_index=unique(pre$loc[,2])[1:7]
Z_index=unique(pre$loc[,3])[c(2,3)]


index=which(Locations_temp[,1]%in%X_index &
              Locations_temp[,2]%in%Y_index & 
              Locations_temp[,3]%in%Z_index )



  
  
pre <- v.est(dwi.obs=dwi.obs[,index]*exp(matrix(rnorm(224*33,0,0.5),33,224)), sigma=sigma, S0=S0const, b=1, grad.mat=grad.mat,
             braingrid=braingrid[,1:16,1:7,2:3], cpus=ncpu, opt=list())
v.res <- v.smooth.bw(pre=pre, braingrid=braingrid[,1:16,1:7,2:3], len=50,
                     range1=c(0.005,0.3), range2=c(0.005,0.3), cpus=ncpu,
                     xy.truncrange=100, z.truncrange=100, thres.v0.w=0.2, K=5,
                     cv.method="mad") ##mcv
v.obj <- update.v.obj(pre=pre, res=v.res)



v.obj

}


save.image("sim0_5.Rdata")