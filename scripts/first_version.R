#  Compute the isotonic regression of numeric vector 'x', with
#  weights 'wt', with respect to simple order.  The pool-adjacent-
#  violators algorithm is used.  Returns a vector of the same length
#  as 'x' containing the regression.
#  02 Sep 1994 / R.F. Raubertas
pava <- function (x, wt=rep(1,length(x)))
{
  n <- length(x)
  if (n <= 1) return (list(estim=x,levelsets = 1))
  if (any(is.na(x)) || any(is.na(wt))) {
    stop ("Missing values in 'x' or 'wt' not allowed")
  }
  lvlsets <- (1:n)
  repeat {
    viol <- (as.vector(diff(x)) < 0)  # Find adjacent violators
    if (!(any(viol))) break
    
    i <- min( (1:(n-1))[viol])        # Pool first pair of violators
    lvl1 <- lvlsets[i]
    lvl2 <- lvlsets[i+1]
    ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
    x[ilvl] <- sum(x[ilvl]*wt[ilvl]) / sum(wt[ilvl])
    lvlsets[ilvl] <- lvl1
  }
  list( estim = x, levelsets = lvlsets)
}
# example 
# pava(c(22.5,23.33,20.833,24.25),wt=c(3,3,3,2))


#  One Simulation
one_sml <- function(scenario,ssize,Tau,m,delta)
{
  maxdose <- nrow(scenario)
  p0 <- scenario[,1]; p1 <- scenario[,2]; p2 <- scenario[,3]; p3 <- scenario[,4]
  cumpr=cbind(p0,p0+p1,p0+p1+p2,1)
  respmat <- matrix(0,maxdose,4)
  trials <- rep(0,maxdose)
  dosenum <- 1
  count <- 0
  phat0 <- rep(0,maxdose)
  phat123 <- matrix(0,maxdose,3)
  
  while (count < ssize) {
    
    # simulate the result for one patient
    otvet <- rank(c(runif(1),cumpr[dosenum,1:3]))[1]
    respmat[dosenum,otvet] <- respmat[dosenum,otvet]+1  
    trials[dosenum] <-trials[dosenum]+1
    
    # isotonic estimation of toxicity and response probability
    dwp <- which(trials!=0)     # doses have at least one patient
    phat0[dwp] <- pava((respmat/trials)[dwp,1])$estim
    for (j in dwp) {
      if (respmat[j,2]>=respmat[j,4]) {phat123[j,] <- pava((respmat/trials)[j,4:2])$estim[3:1]}
      else {phat123[j,] <- pava((respmat/trials)[j,2:4])$estim}
    }
    
    # decide next dosenum
    if (trials[dosenum]<m) {
      dosenum <- dosenum
    } else if (length(which((1-phat0-phat123[,1]-phat123[,2])>delta))>0) {
      dosenum <- min(which((1-phat0-phat123[,1]-phat123[,2])>delta))
    } else if (phat123[dosenum,1]+phat123[dosenum,2]>phat0[dosenum]) {
      dosenum <- min(dosenum+1,maxdose)
    } else if (phat123[dosenum,1]+phat123[dosenum,2]<=phat0[dosenum]) {
      diff=phat0-phat123[,1]-phat123[,2]
      if (length(which(diff==0))>0) {dosenum <- min(which(diff==0))}
      else {
        for (k in 1:max((dosenum-1),1)) {
          if (diff[k]<0 & diff[k+1]>0) {
            if (abs(diff[k])>diff[k+1]) {dosenum <- min(k+1,maxdose); break}
            else {dosenum <- k; break}
          }
        }
      }
    }
    
    if (length(phat0)>maxdose) {
      while (phat0[dosenum]>Tau & dosenum>1) {dosenum <- dosenum-1} 
    }
    
    count=count+1
  }
  
  # find the optimal dose
  toltox=ifelse(phat0<0.5,1,0)
  score <- toltox*(phat123[,1]+2*phat123[,2]+3*phat123[,3])
  bestdose <- order(score)[maxdose]
  (bestdose)
}


#  Simulation Results
scenario31 = matrix(c(0.1,0.72,0.09,0.09,0.2,0.32,0.24,0.24,0.3,0.07,0.07,0.56),byrow=TRUE,ncol=4)
scenario32 = matrix(c(0.15,0.6375,0.17,0.0425,0.3,0.42,0.21,0.07,0.45,0.0275,0.11,0.4125),byrow=TRUE,ncol=4)
scenario33 = matrix(c(0.2,0.56,0.16,0.08,0.4,0.06,0.18,0.36,0.7,0.03,0.03,0.24),byrow=TRUE,ncol=4)

scenario41 = matrix(c(0.1,0.72,0.09,0.09,0.2,0.32,0.24,0.24,0.3,0.07,0.07,0.56,0.4,0.06,0.06,0.48),byrow=TRUE,ncol=4)
scenario42 = matrix(c(0.15,0.6375,0.17,0.0425,0.3,0.42,0.21,0.07,0.45,0.0275,0.11,0.4125,0.6,0.02,0.06,0.32),byrow=TRUE,ncol=4)
scenario43 = matrix(c(0.2,0.56,0.16,0.08,0.4,0.06,0.18,0.36,0.7,0.03,0.03,0.24,0.8,0.02,0.02,0.16),byrow=TRUE,ncol=4)

scenariob1 = matrix(c(0.1,0.00,0.45,0.45,0.2,0.32,0.24,0.24,0.3,0.07,0.28,0.35,0.4,0.4,0.1,0.1),byrow=TRUE,ncol=4)
scenariob2 = matrix(c(0.15,0.55,0.2,0.1,0.3,0.00,0.35,0.35,0.45,0.2,0.2,0.15,0.6,0.4,0.00,0.00),byrow=TRUE,ncol=4)
scenariob3 = matrix(c(0.0,0.6,0.4,0.0,0.1,0.7,0.2,0.0,0.7,0.1,0.1,0.1,0.8,0.2,0.0,0.0),byrow=TRUE,ncol=4)

scenario61 = matrix(c(0.1,0.5,0.4,0.0,0.2,0.3,0.3,0.2,0.3,0.1,0.25,0.35,0.4,0,0.1,0.5,0.5,0,0,0.5,0.6,0,0,0.4),byrow=TRUE,ncol=4)
scenario62 = matrix(c(0,1,0,0,0.1,0.9,0,0,0.2,0.7,0.1,0,0.3,0.4,0.2,0.1,0.4,0.1,0.2,0.3,0.5,0,0,0.5),byrow=TRUE,ncol=4)
scenario63 = matrix(c(0.15,0.4,0.35,0.1,0.3,0.3,0.25,0.15,0.45,0.1,0.15,0.3,0.5,0,0,0.5,0.75,0,0,0.25,0.9,0,0,0.1),byrow=TRUE,ncol=4)
scenario64 = matrix(c(0.3,0.3,0.2,0.2,0.35,0.3,0.2,0.15,0.4,0.15,0.15,0.3,0.45,0.05,0.1,0.4,0.5,0,0,0.5,0.55,0,0,0.45),byrow=TRUE,ncol=4)


result = matrix(0,10000)
for (i in 1:10000) {
  result[i]=one_sml(scenario64,ssize=30,Tau=0.5,m=3,delta=0.7)
}
a=hist(result)
a$counts

# 还需要将DOSE level数目增加
