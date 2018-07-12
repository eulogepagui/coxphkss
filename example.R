#####################################################
## median bias reduction for Cox regression      ####
#####################################################

source("coxphkss.R")

dataCens <- data.frame(time=c(6,7,1,2,3,4,5), 
                       status=c(1,1,1,0,1,1,0), 
                       x=c(0,2,1,1,1,0,0), 
                       sex=c(0,0,0,0,1,1,1)) 

dataNoCens <- data.frame(time=c(6,7,1,2,3,4,5), 
                         status=c(1,1,1,1,1,1,1), 
                         x=c(0,2,1,1,1,0,0), 
                         sex=c(0,0,0,0,1,1,1)) 



####################################
##  maximum likelihood method fit ##
####################################

ml<-coxph(Surv(time, status) ~ x + sex, data=dataCens ) 
summary(ml)
ml2<-coxph(Surv(time, status) ~ x + sex, data=dataNoCens ) 
summary(ml2)


####################################
##  mean bias reduced fit (Firth) ##
####################################

meanBR<-coxphf(Surv(time, status) ~ x + sex, data=dataCens ) 
summary(meanBR)
meanBR2<-coxphf(Surv(time, status) ~ x + sex, data=dataNoCens ) 
summary(meanBR2)

######################################
##  median bias reduced fit (Firth) ##
######################################

medianBR <- coxphkss(beta=c(0,0),data=dataCens,type = "AS_median")
medianBR
medianBR2 <- coxphkss(beta=c(0,0),data=dataNoCens,type = "AS_median")
medianBR2




