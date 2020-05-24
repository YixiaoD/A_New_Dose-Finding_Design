omia3 = function(xy)
{
  maxdose=3
  ssize=xy[1]
  scenario=xy[2]
  if (scenario ==1) {
    k0=c(0.1,0.2,0.3)
    k1=c(0.72,0.32,0.07)
    k2=c(0.09,0.24,0.07)
    k3=c(0.09,0.24,0.56) 
  }
  if (scenario ==2) {
    k0=c(0.15,0.3,0.45)
    k1=c(0.6375,0.42,0.0275)
    k2=c(0.17,0.21,0.11)
    k3=c(0.0425,0.07,0.4125) 
  }
  if (scenario ==3) {
    k0=c(0.2,0.4,0.7)
    k1=c(0.56,0.06,0.03)
    k2=c(0.16,0.18,0.03)
    k3=c(0.08,0.36,0.24) 
  }
  
  cumpr=cbind(k0,k0+k1,k0+k1+k2,1)
  count = 0
  respmat=matrix(0,maxdose,4)
  phattox=rep(0,maxdose)
  phatsucc=rep(0,maxdose)
  phatnrnt=rep(0,maxdose)
  trials=rep(0,maxdose)
  dosenum=1 #first group assigned to dose1
  maxdose1=maxdose
  mindose1=1
  
  while (count < ssize) {
    otvet=rank(c(runif(1),cumpr[dosenum,1:3]))[1]
    respmat[dosenum,otvet]=respmat[dosenum,otvet]+1
    otvet1=otvet
    otvet=rank(c(runif(1),cumpr[dosenum,1:3]))[1]
    respmat[dosenum,otvet]=respmat[dosenum,otvet]+1
    otvet2=otvet
    otvet=rank(c(runif(1),cumpr[dosenum,1:3]))[1]
    respmat[dosenum,otvet]=respmat[dosenum,otvet]+1
    otvet3=otvet
    trials[dosenum]=trials[dosenum]+3
    phattox=respmat[,1]/(trials+0.001)
    
    if (phattox[dosenum]>.5) {dosenum=max(dosenum-1,mindose1)} else {
      if ((otvet1==1 & otvet2==1) || (otvet1==1 & otvet3==1) || (otvet2==1 & otvet3==1)) {dosenum=max(dosenum-1,mindose1)}
      if ((otvet1==2 & otvet2==2) || (otvet1==2 & otvet3==2) || (otvet2==2 & otvet3==2)) {dosenum=min(dosenum+1,maxdose1)}
    }
    
    # if ((otvet1==1 & otvet2==1) || (otvet1==1 & otvet3==1) || (otvet2==1 & otvet3==1)) {dosenum=max(dosenum-1,mindose1)}
    # if ((otvet1==2 & otvet2==2) || (otvet1==2 & otvet3==2) || (otvet2==2 & otvet3==2))  {dosenum=min(dosenum+1,maxdose1)}
    # if (phattox[dosenum]>.5) {dosenum=max(dosenum-1,mindose1)}
    count=count+3
  }
  
  #Ivanova's definition of optimal dose
  N=ssize/maxdose
  ind=ifelse(trials>=N,1,0)
  toltox=ifelse(phattox<0.5,1,0)
  besttrt=order((respmat[,3]+respmat[,4])*toltox*ind/(trials+0.00001))[maxdose]

  # #O'Qquigley's definition of optimal dose
  # toltox=ifelse(phattox<0.5,1,0)
  # besttrt=order(toltox*((respmat[,2]/(trials+0.001))+2*respmat[,3]/(trials+0.001)+3*respmat[,4]/(trials+0.001)))[maxdose]

  c(trials/sum(trials),sum(trials),besttrt)
}

#simulation 10000 times
result = matrix(0,10000,5)
for (i in 1:10000) {
  result[i,]=omia3(c(30,1))
}


a=hist(result[,5])
a$counts
mean(result[,1]);mean(result[,2]);mean(result[,3])



