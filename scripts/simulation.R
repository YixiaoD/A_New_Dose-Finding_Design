omia3 = function(xy)
{
  maxdose=3
  ssize=xy[1]
  scenario=xy[2]
  if (scenario==1)	{
    prtox=c(0.1,0.2,0.3)	
    prsucc=c(0.18,0.48,0.63)	#k>=2:0.18,0.48,0.63 k=3:0.09,0.24,0.56
  }
  if (scenario==2)	{
    prtox=c(0.15,0.3,0.45)	
    prsucc=c(0.2125,0.28,0.5225)  #k>=2:0.2125,0.28,0.5225 k=3:0.0425,0.07,0.4125
  }
  if (scenario==3)	{
    prtox=c(0.2,0.4,0.7)
    prsucc=c(0.24,0.54,0.27)  #k>=:0.24,0.54,0.27 k=3:0.08,0.36,0.24
  }

  cumpr=cbind(prtox,prsucc+prtox,1)
  count = 0
  respmat=matrix(0,maxdose,3)
  phattox=rep(0,maxdose)
  phatsucc=rep(0,maxdose)
  phatnrnt=rep(0,maxdose)
  trials=rep(0,maxdose)
  dosenum=2
  stopdose=rep(0,maxdose)
  stop=0
  maxdose1=maxdose
  mindose1=1
  
  while (count < ssize) {
    otvet=rank(c(runif(1),cumpr[dosenum,1:2]))[1]
    respmat[dosenum,otvet]=respmat[dosenum,otvet]+1
    otvet1=otvet
    otvet=rank(c(runif(1),cumpr[dosenum,1:2]))[1]
    respmat[dosenum,otvet]=respmat[dosenum,otvet]+1
    otvet2=otvet
    otvet=rank(c(runif(1),cumpr[dosenum,1:2]))[1]
    respmat[dosenum,otvet]=respmat[dosenum,otvet]+1
    otvet3=otvet
    trials[dosenum]=trials[dosenum]+3
    phattox=respmat[,1]/(trials+0.001)
    phatsucc=respmat[,2]/(trials+0.001)
    phatnrnt=respmat[,3]/(trials+0.001)
    if (phattox[dosenum]>.3) {dosenum=max(dosenum-1,mindose1)} else {
      if (otvet1*otvet2==1 || otvet2*otvet3==1 || otvet1*otvet3==1) {dosenum=max(dosenum-1,mindose1)}
      if (otvet1*otvet2==9 || otvet2*otvet3==9 || otvet1*otvet3==9)  {dosenum=min(dosenum+1,maxdose1)}
    }
    count=count+3
  }
  
  N=ssize/maxdose
  ind=ifelse(trials>=N,1,0)
  besttrt=order(respmat[,2]*ind/(trials+0.00001))[maxdose]
  c(trials/sum(trials),sum(trials),sum(respmat[,2]),besttrt)
}

result = matrix(0,10000,6)
for (i in 1:10000) {
  result[i,]=omia3(c(24,3))
}

a=hist(result[,6])
a$counts
mean(result[,2])




