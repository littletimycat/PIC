######### Simulation Study ###########################################################

##### Notation ######################
## n: sample size
## T: failure time
## S(t): time-dependent covariate
## Z: p-dim time-independent covariate vector
## delta: indicator of observing exact T
## L: left endpoint of censoring interval for T if delta = 0
## R: right endpoint of censoring interval for T if delta = 0
## beta: regression coefficient for S(t) in Cox model
## gamma: p-dim regression coefficient vector for Z in Cox model
## lambda: baseline hazard function in Cox model
## tau: length of study
## pt: proportion of exact obs among failures
## nsim: number of simulations


n =200   # sample size
tau = 0.8                                   # length of study
pt = 0.2                                           # proportion of exact obs among failures

beta = 0.5
gamma = -0.5
lambda = function(t) {return(1/(2+t))} # why it is 1/(2+t)
p = length(gamma)



max.iter = 5000                                   # max number of iterations in EM
max.tol = 1e-4                                    # max tolerance of parameter differences in EM
threshold = 0.01                                  # threshold for rescaling absolute differences in EM

beta.initial = 0                                  # initial values of parameters in EM
gamma.initial =  rep(0,p)

hn = 5/sqrt(n)                                    # tuning parameter for variance estimation

q = rep(NA, n)

data = generate.data(n, beta, gamma, lambda, tau, pt)
data$Rs = data$L*(data$R==Inf)*(1-data$delta1 - data$delta2)+data$T*data$delta1
data$Rs[(data$R<Inf)&(data$delta1==0)] = data$R[(data$R<Inf)&(data$delta2==1)]


tk = sort(unique(c(data$L*(1-data$delta1),data$R*(data$delta2),data$T*data$delta1)))
tk = tk[-c(1,length(tk))]
nt = length(tk)

lambda.initial = rep(1/nt,nt)
S.all = matrix(data$S, nrow=n, ncol=nt)

#using logistic regression to estimate pi
ind = which(data$delta1==data$delta2)

# we need to do modification about the 

fit = glm(data$eta[ind] ~ data$Z[ind], data = data, family = binomial())
zz = data$Z
q[ind] =  (exp(fit$coefficients[1]+fit$coefficients[2]*zz[ind]))/(1+exp(fit$coefficients[1]+fit$coefficients[2]*zz[ind]))
q[-ind] = 1


#estimation part
res = estbg(beta.initial,gamma.initial,lambda.initial,S.all,tk,data,max.tol,max.iter,threshold)
ress = eestbg(beta.initial,gamma.initial,lambda.initial,S.all,tk,data,max.tol,max.iter,threshold,q)

nq= res$qf 
nqq = ress$qf 

bg.est = c(res$betaf,res$gammaf)
bgg.est= c(ress$betaf,ress$gammaf)
#
etaa = data$eta
# new log-likelihood function is 
wei_logistic = function(alpha,u){
  l_star = 0
  for(k in ind){
    l_star = l_star + u[k]*(etaa[k]*(alpha[1] + alpha[2]*zz[k])-log(1+exp(alpha[1]+alpha[2]*zz[k])))
  }
  l_star = -l_star
}

#the thing is that every time, the q1
# we need to generate 200 times bootstrap i.e. 200 times new exp(1)
nboot=200
q1 = rep(NA, n)
bggb_est = matrix(0, ncol = 2, nrow = nboot)
see.est = rep(NA,2)
for(j in 1:nboot){
  u = rexp(n,1)
  est_res = nlm(wei_logistic, c(-1,0.5), u)
  wls = est_res$estimate#weighted_logistic_est
  q1[ind] =  c(exp(wls[1]+wls[2]*zz[ind])/(1+exp(wls[1]+wls[2]*zz[ind])))
  q1[-ind] = 1
  ressb = eestbgu(ress$betaf,ress$gammaf,ress$lambdaf,S.all,tk,data,max.tol,max.iter,threshold,q1, u)
  bggb_est[j,1] = ressb$betaf
  bggb_est[j,2] = ressb$gammaf
  #conv.flagb[j] = bool((nqb[j]>0)&(nqb[j]<max.iter));
  #cout<< j <<endl;
  #print(j)
}

#get standard deviation
see.est[1]= sd(bggb_est[,1])
see.est[2] = sd(bggb_est[,2])

bias_ipw = bgg.est - c(beta,gamma)
bias_full = bg.est - c(beta,gamma)

indd = rep(NA,2)
indd[1] = ifelse(bgg.est[1]-1.96*see.est[1]<=beta & bgg.est[1]+1.96*see.est[1]>=beta,1,0)
indd[2] = ifelse(bgg.est[2]-1.96*see.est[2]<=gamma & bgg.est[2]+1.96*see.est[2]>=gamma,1,0)

output = t(rbind(bias_ipw, bias_full,see.est,indd))
rownames(output) = c("beta","gamma")
colnames(output) = c("IPW","Full sample","ese","cp_ind")
write(output, file = filename, sep=",", ncol=length(output),append = length(output))



#conv.flag = as.numeric((nq>0)&(nq<max.iter)&(nqq>0)&(nqq<max.iter))
#n.conv = sum(conv.flag)  # number of converged replicates

#index = (conv.flag==1)&(apply(see.est,1,sum)!=0)
#bias = apply(bgg.est[index,],2,mean)-c(beta,gamma)
#sdd = apply(bg.est[index,],2,sd) # sd for the full sample
#ssd = apply(bgg.est[index,],2,sd)
#ese = apply(see.est[index,],2,mean)
#ess = apply(se.est[index,],2,mean)
#upper = t(bgg.est[index,]+1.96*see.est[index,])
#lower = t(bgg.est[index,]-1.96*see.est[index,])
#cp = apply(as.matrix((c(beta,gamma)<=upper)&(c(beta,gamma)>=lower)),1,mean)

#est.results = as.data.frame(cbind(bias,ssd,sdd,ese,ess,cp))
#row.names(est.results) = c("beta",rep("gamma",p))
#colnames(est.results) = c("bias", "ssd", "fullsample_ssd","ese","fullsample_ese", "cp")

#print("=====================bootstrap samples:===================")
#nboot 
#print("====================sample size is:===================")
#n 
#print("===================right censoring rate is:=======================")
#right.censor.avg
#print("=====================relaive efficiency is about:=====================")
#(sdd/ssd)^2
#print("=======================Estimating results:===================")
#round(est.results,3)  # simulation results
