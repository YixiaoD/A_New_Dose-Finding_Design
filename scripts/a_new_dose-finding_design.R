
######## Splus code to run one trial from Ivanova, Biometrics 2003, 59, 1003-1009.
######## with  subjects assinged in groups of 3
######## I hope this is working... 


omia3=function(xy)
{	
  model=xy[1]
  maxdose=4
  stage=xy[2]
  e1=xy[3]
  e2=xy[4]
  stoptype=xy[5]
  alp=c(0.05, 0.10, 0.2, 0.3)
  alp1=c(.7,.5,.3,0.1)
#  ra=xi.alpha0
  if (model==1)	{prtox=c(.06,.17,.25,.3)	#pr(tox), rk
    prsucc=c(.2,.7,.6,.5)	#pr(resp & no tox), pk
  }
  if (model==2)	{prtox=c(.13,.3,.4,.5)	
    prsucc=c(.7,.5,.5,.4)
  }
  if (model==3)	{prtox=c(0,.05,.15,.3)
    prsucc=c(.1,.3,.7,.5)
  }
  if (model==4)	{prtox=c(0,0,.1,.14)	
    prsucc=c(.2,.3,.5,.7)
  }
  
  #prtox_c(.05,.1,.2,.4)
  #prsucc_1-prtox-c(.25,0.7,.35,.45)  
  
  p0=.4
  p1=.6
  cumpr=cbind(prtox,prsucc+prtox,1)
  count <- 0
  respmat=matrix(0,maxdose,3) 
  trials=rep(0,maxdose)
  #y_rep(0,maxdose)
  dosenum=2
  stopdose=rep(0,maxdose)
  stop=0
  maxdose1=maxdose
  mindose1=1
  #  random walk
  while (count < stage && stop==0) {
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
    phattox=respmat[,1]/(trials+.001)
    y=respmat[,3]
    phatnrnt=respmat[,3]/(trials+.001)
    if (otvet1*otvet2==1 || otvet2*otvet3==1 || otvet1*otvet3==1) {dosenum=max(dosenum-1,mindose1)}
    if (otvet1*otvet2==9 || otvet2*otvet3==9 || otvet1*otvet3==9)  {dosenum=min(dosenum+1,maxdose1)}
    if (phattox[dosenum]>.3) {dosenum=max(dosenum-1,mindose1)}
    if (stoptype==1) {#  SPRT
      wd=respmat[dosenum,2]*log(p1/p0)+(trials[dosenum]-respmat[dosenum,2])*log((1-p1)/(1-p0))
      if (wd<log(e2/(1-e1)))  {stopdose[dosenum]=sum(trials)
        if (dosenum==1) {mindose1=2}
        if (dosenum==4) {maxdose1=3}
      }
      if (wd>log((1-e2)/e1)) {stop=1; besttrt=dosenum}
    }
    if (stoptype>5) {if (max(trials)>=stoptype) {	stop=1
      ind7=ifelse(trials>=7,1,0)
      besttrt=order(respmat[,2]*ind7/(trials+0.00001))[maxdose]
    }
    }
    count=count+3
    #print(round(c(count,otvet1,otvet2,dosenum),2))
  }
  if (stop==0) {	ind7=ifelse(trials>=min(7,max(trials)-1),1,0)
    besttrt=order(respmat[,2]*ind7/(trials+0.00001))[maxdose]
  }							   
  c(trials/sum(trials),sum(trials),sum(respmat[,2]),besttrt)
}

