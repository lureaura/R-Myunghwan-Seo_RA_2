
#########################################   INPUT   #########################################
# r1   : simulation number
# sam  : sample size
# ap0(alpha) : threshold model coefficients
# gm0(gamma) : threshold parameters
# rho0 : Candidates for AR coefficient value(vector or scalar)
# tht0 : Candidates for MA coefficient value(vector or scalar)
# b    : value of block length(scalar)  
# lag : value of lag(scalar)
# mu   : level of the model (constant, should be 0 if alpha == 0 )
# mod.type   : "size" or "power" (Size : table1 / Power : table1)

#########################################   OUTPUT   ########################################
# p-value of each cases according to above parameters



###############################    PARAMETER VALUE SETTING   ################################
r1 <-50
sam <- 250
gm0 <- c(4,8)
rho0 <- c(0,-0.5,0.5,0,0)
tht0 <- c(0,0,0,-0.5,0.5)
b <- 6  #6 or 8 in the paper
lag<-3   #3 or 6 in the paper
mu <- 0 
mod.type<-"power" #or "power"

###################################      SIMULATION     #####################################
#getwd()
setwd("C:/R")
source("Functions_unit_simple.R")
ptm <- proc.time()
simulation2.b1(r1, sam, gm0, rho0, tht0, b, lag, mu, mod.type)
time<- proc.time() - ptm
for(i in 1:3){
  hour<-min<-sec<-NULL
  min<-floor(as.numeric(time[i])/60);  sec<-round(as.numeric(time[i])-min*60,5)
  hour<-floor(min/60); min<-min-hour*60
  time[i]<-paste0(hour, "h ", min, "m ", sec, "s")
}
time

