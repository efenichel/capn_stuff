# -----------------------------------------------------------------------------
# File Name: system_fns.R
# Purpose: Groundwater Function Definitions
# Author:  Eli Fenichel Edits: Ethan Addicott
# Date:  11/06/19
# Notes: Data for capN are as follows
# [1] crop share parameters
# [2] crop amoount parameters
# [3] water withdrawal constant alpha
# [4] water withdrawal coefficient beta (on wdep_)
# [5] water withdrawal coefficients on acres planted and acres planted squared
# [6] water withdrawal coefficients on acres planted to each crop
# [7] water withdrawal coefficients on acres planted squared
# [8] crop prices
# [9] cost of each crop per acre
#------------------------------------------------------------------------------

# ####functions from Fenichel et al. 2016
#####################################################################
#Area planting decision model
cropFwater<- function(water, param){
  
  cropSharePar<-param[[1]]
  cropAmountPar<-param[[2]]
  
  contr <- exp(cropSharePar$crop.share.a+water*cropSharePar$crop.share.b)
  h <- sum(contr)
  contr <- c(contr, 1)
  
  cropout <- matrix(0,nrow=length(contr),ncol=6)
  
  for(j in 1:6){
    cropout[,j]<-contr/(1+h)*cropAmountPar[,j]
  }
  
  rowSums(t(cropout))
  
}


#####################################################################
#Derivative of area planting decision model
DcropFwater<- function(water, param){
  
  cropSharePar<-param[1][[1]]
  cropAmountPar<-param[2][[1]]
  
  contr <- exp(cropSharePar$crop.share.a+water*cropSharePar$crop.share.b)
  h <- sum(contr)
  hp <- sum(contr*cropSharePar$crop.share.b)
  
  contr <- c(contr, 1)
  bv <- c(cropSharePar$crop.share.b,0)
  
  Dcropout <- matrix(0,nrow=length(contr),ncol=6)
  
  for(j in 1:6){
    Dcropout[,j]<- cropAmountPar[,j]*contr*(bv+bv*h-hp)/((1+h)^2)
  }
  
  rowSums(t(Dcropout))
  
}

#####################################################################
#Water withdrawl function

Wwd<-function(water,param){
  
  #alpha<-as.data.frame(lapply(param[[4]],identity))
  alpha<-param[[3]]
  #beta<-as.data.frame(lapply(param[[5]],identity))
  beta<-param[[4]]
  Gamma<-param[[5]]
  
  AP1<-cropFwater(water,param) #area planted
  AP <- matrix(0,nrow=1,ncol=6) #normalized area planted
  for(j in 1:6){
    AP[j]<-AP1[j]/sum(AP1)
  }
  
  #area planted and arrea planted squared
  z<-matrix(0,nrow=1,ncol=10)
  z[1]=AP1[1]
  z[2]=(AP1[1])^2
  z[3]=AP1[2]
  z[4]=(AP1[2])^2
  z[5]=AP1[3]
  z[6]=(AP1[3])^2
  z[7]=AP1[4]
  z[8]=(AP1[4])^2
  z[9]=AP1[5]
  z[10]=(AP1[5])^2
  
  wx<-exp(alpha+water*beta+sum(z*Gamma))/sum(AP1) #projects based on regession model
  min(wx,water) #imposes no shorting water
  
}
#####################################################################
#Derivative of water withdrawl function
WwdDs1<-function(water,param){
  Wwd <- Wwd(water,param)
  
  beta<-param[[4]][[1]]
  Gamma1<-param[[6]]
  Gamma2<-param[[7]]
  
  AP1<-cropFwater(water,param) #area planted
  APprime1 <- DcropFwater(water,param) #derivative of the area planted function
  
  temp<-sum(APprime1[1:5]*Gamma1)+sum((2*APprime1[1:5]*AP1[1:5])*Gamma2)
  Wwd*(beta+temp)
  
}
#####################################################################
#profit function
profit <- function(water,param){
  
  prices<-t(param[[8]])
  costCropAcreX<-param[[9]]
  
  AP1<-cropFwater(water,param) #area planted
  AP <- matrix(0,nrow=1,ncol=6) #normalized area planted
  for(j in 1:6){
    AP[j]<-AP1[j]/sum(AP1)
  }
  temp<-c(as.numeric(prices[1:5,1]-costCropAcreX),0)
  sum(AP[1:6]*(temp[1:6]))
  
  
}

#####################################################################
#profit prime s function
ProfDs1 <- function(water,param){
  
  prices<-t(param[[8]])
  costCropAcreX<-param[[9]]
  
  AP1<-cropFwater(water,param) #area planted
  APprime1 <- DcropFwater(water,param) #derivative of the area planted function
  
  
  APprime1pa <- matrix(0,nrow=1,ncol=6) #normalized area planted
  for(j in 1:6){
    APprime1pa[j]<-APprime1[j]/sum(AP1)
  }
  
  temp<-c(as.numeric(prices[1:5,1]-costCropAcreX),0)
  sum(APprime1pa[1:6]*(temp[1:6]))
  
  
}

#####################################################################
#change in groundwater with respect to time
sdot <- function(water,recharge,param){
  recharge/12 - Wwd(water,param)
}
######################################################################