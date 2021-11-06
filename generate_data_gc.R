

generate.data = function(n, beta, gamma, lambda, tau, pt)# why we need lambda
{
  U1 = runif(n, min=0, max=1)
  #U2 = runif(n, min=0, max=1)
  U3 = runif(n, min=0, max=1)
  q = rep(NA, n)
  
  
  S = U1+2*U3
  Z = matrix(0, nrow=n, ncol=1)
  Z[,1] = rbinom(n,1,0.5)
  # Z[,2] = rnorm(n, mean=0, sd=1)
  # need eta to denote selected obs in RC data 
  T = L = R = delta1 = delta2 = rep(0,n)   #T is the event time
  eta = rep(0, n)
  
  for(i in 1:n)
  {
    temp1 = runif(1, min=0, max=1)
    T[i] = (exp(-log(1-temp1)*exp(-beta*S[i]-sum(gamma*Z[i,])))-1)*2 # how?
    
    temp2 = 0
    while(temp2[length(temp2)]<tau) {# whats the meaning in while() ? ans: ensures the last element less than censoring
      temp2 = c(temp2, temp2[length(temp2)]+runif(1,min=0.1,max=tau/2)) # why use 0.1 as left limit, we can change tau/2 by tau/3
    }
    #temp3 is the censoring variable
    temp3 = c(temp2[-length(temp2)], tau) #last element of temp2 should greater than 6 and hence removed, -length(temp2), and append tau
    
    L[i] = max(temp3[temp3<T[i]])
    
    if(T[i]<=tau) {
      if(rbinom(1,size=1,prob=0.1)) {
        R[i] = Inf #add more case of right censoring
        delta1[i] = delta2[i] = 0 
      } else {
        R[i] = min(temp3[temp3>=T[i]])
        delta2[i] = 1
        delta1[i] = 0
      }
    } else {
      R[i] = Inf
      delta1[i] = delta2[i] = 0 
    }
    
    if(R[i]<Inf) {
      delta1[i] = rbinom(1, size=1, prob=pt)
      delta2[i] = 1 - delta1[i]
    }
    
    # the probability of missing is defined as
    #only defined for non-case, need revise
    q[i] = exp(0.3-2*Z[i,1])/(1+exp(0.3-2*Z[i,1]))
    #print(q[i])
    
    eta[i] = delta1[i] + delta2[i] + (1-delta1[i]-delta2[i])*rbinom(1,1,q[i])
  }
  
  return(list(delta = delta1, delta1=delta1,delta2 = delta2,eta = eta,T=T,L=L,R=R,S=S,Z=Z))
}
