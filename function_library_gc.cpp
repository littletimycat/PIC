#include <RcppArmadillo.h>
#include <iostream>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat WWmean(double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data) {
  int n = Shat.n_rows, nt = Shat.n_cols;
  mat Z = data["Z"];
  colvec delta1 = data["delta1"]; // censoring indicator
  colvec delta2 = data["delta2"]; // censoring indicator
  colvec T = data["T"];// exact happening time
  colvec L = data["L"];// interval censoring left
  colvec R = data["R"];//interval censoring right
  colvec gz = Z*gamma; // ? covariate times coefficient
  colvec eta = data["eta"];
  mat res = zeros(n,nt);
  for(int i=0; i<n; i++) {
    if(delta1[i]==1) {
      for(int k=0; k<nt; k++) {
        bool index = (tk[k]==T[i])&(is_finite(Shat(i,k))); // Shat(i,k) 
        if(index) {
          res(i,k) = 1;
        }
      }
    } else if(is_finite(R[i])) {
      colvec num = lambda%trans(exp(gz[i]+beta*Shat.row(i)));
      double temp = 0;
      for(int k=0; k<nt; k++) {
        bool index = (tk[k]>L[i])&(tk[k]<=R[i])&(is_finite(Shat(i,k)));
        if(index) {
          res(i,k) = num[k];
          temp += num[k];
        }
      }
      if(temp!=0){
        res.row(i) = res.row(i)/(1-exp(-temp));
      }
    }
  }
  return res;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
// need to introduce the weight
colvec llambdak(double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data, colvec qq) {
  int n = Shat.n_rows, nt = Shat.n_cols;
  colvec delta1 = data["delta1"];
  colvec delta2 = data["delta2"];
  mat Z = data["Z"];
  colvec Rs = data["Rs"];
  colvec eta = data["eta"];
  mat W = WWmean(beta,gamma,lambda,Shat,tk,data);
  colvec gz = Z*gamma;
  colvec l = zeros(nt);
  colvec wt = delta1 + delta2 + eta%(1-delta1-delta2)/qq;
  for(int k=0; k<nt; k++) {
    colvec e = exp(gz+beta*Shat.col(k));
    double temp1 = 0, temp2 = 0;
    for(int i=0; i<n; i++) {
      bool index = (tk[k]<=Rs[i])&(is_finite(Shat(i,k)));
      if(index) {
        temp1 += wt[i]*W(i,k);
        temp2 += wt[i]*e[i];
      }
    }
    if(temp2>0) {
      l[k] = temp1/temp2;
    }
  }
  return l;
}




// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
colvec sscorebg(double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data, colvec qq) {
  int n = Shat.n_rows, nt = Shat.n_cols;
  int p = gamma.n_elem;
  colvec delta1 = data["delta1"];
  colvec delta2 = data["delta2"];
  colvec eta = data["eta"];
  mat Z = data["Z"];
  colvec Rs = data["Rs"];
  mat W = WWmean(beta,gamma,lambda,Shat,tk,data);
  colvec gz = Z*gamma;
  colvec lbg = zeros(p+1);
  colvec wt = delta1 + delta2 + eta%(1-delta1-delta2)/qq;
  for(int k=0; k<nt; k++) {
    
    colvec e = exp(gz+beta*Shat.col(k));
    colvec eb = e%Shat.col(k);// exp(S_jk)*S_jk
    mat eg = e%Z.each_col();// whats means
    double temp1 = 0, temp2 = 0;
    colvec temp3 = zeros(p);
    
    for(int i=0; i<n; i++) {
      //wt[i] = delta1[i] + delta2[i] +  eta[i]*(1-delta1[i]-delta2[i])/0.3;
      bool index = (tk[k]<=Rs[i])&(is_finite(Shat(i,k)));
      if(index) {
        temp1 += wt[i]*e[i];
        temp2 += wt[i]*eb[i];
        temp3 += wt[i]*trans(eg.row(i));
      }
    }
    double tempb = 0;
    colvec tempg = zeros(p);
    for(int i=0; i<n; i++) {
      bool index = (tk[k]<=Rs[i])&(is_finite(Shat(i,k)));
      if(index) {
        //wt[i] = delta1[i] + delta2[i] +  eta[i]*(1-delta1[i]-delta2[i])/0.3;
        tempb += wt[i]*W(i,k)*(Shat(i,k)-temp2/temp1);
        tempg += wt[i]*W(i,k)*(trans(Z.row(i))-temp3/temp1);
      }
    }
    lbg[0] += tempb;
    lbg.subvec(1,p) += tempg;
  }
  return lbg;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat hhessianbg(double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data, colvec qq) {
  int n = Shat.n_rows, nt = Shat.n_cols;
  int p = gamma.n_elem;
  mat Z = data["Z"];
  colvec Rs = data["Rs"];
  colvec delta1 = data["delta1"];
  colvec delta2 = data["delta2"];
  colvec eta = data["eta"];
  colvec wt = delta1 + delta2 + eta%(1-delta1-delta2)/qq;
  mat W = WWmean(beta,gamma,lambda,Shat,tk,data);
  colvec gz = Z*gamma;
  mat lbg = zeros(p+1,p+1);
  
  for(int k=0; k<nt; k++) {
    colvec e = exp(gz+beta*Shat.col(k));// 
    colvec eb = e%Shat.col(k);// b represents for S, time dependent covariate, take 1st order derivative
    //for beta
    mat eg = e%Z.each_col();// g represents for gamma, time independent covariate, take 1st order derivative
    //for gamma 
    colvec ebb = eb%Shat.col(k); // second derivative for beta
    mat ebg = eb%Z.each_col(); // take derivative wrt both beta and gamma
    cube egg = zeros(p,p,n); // second derivative for gamma
    
    for(int i=0; i<n; i++) {
      egg.slice(i) = e[i]*trans(Z.row(i))*Z.row(i);// hessian for Z, how?
    }
    
    double temp1 = 0, temp2 = 0, temp4 = 0;
    rowvec temp3 = zeros<rowvec>(p), temp5 = zeros<rowvec>(p);
    mat temp6 = zeros(p,p);
    
    for(int i=0; i<n; i++) {
      bool index = (tk[k]<=Rs[i])&(is_finite(Shat(i,k)));
      if(index) {
        temp1 += wt[i]*e[i];
        temp2 += wt[i]*eb[i];
        temp3 += wt[i]*eg.row(i);
        temp4 += wt[i]*ebb[i];
        temp5 += wt[i]*ebg.row(i);
        temp6 += wt[i]*egg.slice(i); // add weight in this way is correct?
      }
    }
    
    double tempbb = 0;
    rowvec tempbg = zeros<rowvec>(p);
    mat tempgg = zeros(p,p);
    
    for(int i=0; i<n; i++) {
      bool index = (tk[k]<=Rs[i])&(is_finite(Shat(i,k)));
      if(index) {
        tempbb += wt[i]*W(i,k)*(-temp4/temp1-pow(temp2/temp1,2));
        tempbg += wt[i]*W(i,k)*(-temp5/temp1-temp2*temp3/pow(temp1,2));
        tempgg += wt[i]*W(i,k)*(-temp6/temp1-trans(temp3)*temp3/pow(temp1,2));
      }
    }
    lbg(0,0) += tempbb;
    lbg.submat(0,1,0,p) += tempbg;
    lbg.submat(1,1,p,p) += tempgg;
  }
  lbg.submat(1,0,p,0) = trans(lbg.submat(0,1,0,p));
  
  return lbg;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List eestbg(double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data, double maxtol, double maxiter, double threshold, colvec qq) {
  int p = gamma.n_elem;
  int nt = lambda.n_elem;
  
  colvec beta_hat = zeros(maxiter+2);
  mat gamma_hat = zeros(p,maxiter+2);
  mat lambda_hat = zeros(nt,maxiter+2);
  
  beta_hat[0] = beta-1, beta_hat[1] = beta;
  gamma_hat.col(0) = gamma-1, gamma_hat.col(1) = gamma;
  lambda_hat.col(0) = lambda-1, lambda_hat.col(1) = lambda;
  
  int qf = 0;
  for(int q=1; q<=maxiter; q++) {
    colvec temp1 = zeros(p+nt+1);
    temp1[0] = beta_hat[q];
    temp1.subvec(1,p) = gamma_hat.col(q);
    temp1.subvec(p+1,p+nt) = cumsum(lambda_hat.col(q));
    
    colvec temp2 = zeros(p+nt+1);
    temp2[0] = beta_hat[q-1];
    temp2.subvec(1,p) = gamma_hat.col(q-1);
    temp2.subvec(p+1,p+nt) = cumsum(lambda_hat.col(q-1));
    
    colvec temp3 = abs(temp1-temp2);
    for(int i=0; i<=(p+nt); i++) {
      colvec temp4 = zeros(2);
      temp4[0] = temp1[i], temp4[1] = temp2[i];
      if(min(abs(temp4))>threshold) {
        temp3[i] = temp3[i]/min(abs(temp4));
      }
    }// why
    
    if((max(temp3.subvec(0,p))+max(temp3.subvec(p+1,p+nt)))>maxtol) {
      
      lambda_hat.col(q+1) = llambdak(beta_hat[q],gamma_hat.col(q),lambda_hat.col(q),Shat,tk,data,qq);
      
      colvec sbg = sscorebg(beta_hat[q],gamma_hat.col(q),lambda_hat.col(q+1),Shat,tk,data,qq);
      
      mat hbg = hhessianbg(beta_hat[q],gamma_hat.col(q),lambda_hat.col(q+1),Shat,tk,data,qq);
      
      colvec update = inv(hbg)*sbg;
      beta_hat[q+1] = beta_hat[q]-update[0];
      gamma_hat.col(q+1) = gamma_hat.col(q)-update.subvec(1,p);
      qf += 1;
    } else {
      break;
    }
    
  }
  List res;
  res["betaf"] = beta_hat[qf+1];
  res["gammaf"] = gamma_hat.col(qf+1);
  res["lambdaf"] = lambda_hat.col(qf+1);
  res["qf"] = qf;
  
  return res;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double lloglik(double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data, colvec qq) {
  int n = Shat.n_rows, nt = Shat.n_cols;
  mat Z = data["Z"];
  colvec delta1 = data["delta1"];
  colvec T = data["T"];
  colvec L = data["L"];
  colvec R = data["R"];
  colvec eta = data["eta"];
  colvec gz = Z*gamma;
  double logl = 0;
  for(int i=0; i<n; i++) {
    colvec le = lambda%trans(exp(gz[i]+beta*Shat.row(i)));
    if(delta1[i]==1) {
      double tempt1 = 0, tempt2 = 0, ct = 0;
      for(int k=0; k<nt; k++) {
        bool indext1 = (tk[k]<=T[i])&(is_finite(Shat(i,k)));
        bool indext2 = (tk[k]==T[i])&(is_finite(Shat(i,k)));
        if(indext1) {
          tempt1 += le[k];
          ct += 1;
        }
        if(indext2) {
          tempt2 += le[k];
        }
      }
      if(ct>0){
        logl += log(tempt2)-tempt1;//dont know why
      }
    } else {
      double templ = 0, tempr = 0, cr = 0;
      for(int k=0; k<nt; k++) {
        bool indexl = (tk[k]<=L[i])&(is_finite(Shat(i,k)));
        bool indexr = (tk[k]<=R[i])&(is_finite(Shat(i,k)));
        if(indexl) {
          templ += le[k];
        }
        if(indexr) {
          tempr += le[k];
          cr += 1;
        }
      }
      double ll = exp(-templ);
      double lr = exp(-tempr);
      if(!is_finite(R[i])){// for right censored obs, R_i is infinity
        logl += eta[i]*log(ll)/qq[i];// add weight here, since the weight for case is 1, so omit before and 
        // only add here
      } else if(cr>0){
        logl += log(ll-lr);
      }
    }
  }
  return logl;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
//check the code/formula first
//try to use bpptstrap if it under estimate the true variability
// check the condition of pl
double pplfun(colvec para, colvec lambda, mat Shat, colvec tk, List data, double maxtol, double maxiter, double threshold, colvec qq) {
  int p = para.n_elem-1;
  double beta = para[0];
  colvec gamma = para.subvec(1,p);
  int nt = lambda.n_elem;
  mat plambda_hat = zeros(nt,maxiter+2);
  plambda_hat.col(0) = lambda-1, plambda_hat.col(1) = lambda;
  int qlf = 0;
  for(int q=1; q<=maxiter; q++) {
    colvec temp1 = cumsum(plambda_hat.col(q));
    colvec temp2 = cumsum(plambda_hat.col(q-1)); //compute capital lambda
    colvec temp3 = abs(temp1-temp2);
    for(int i=0; i<nt; i++) {
      colvec temp4 = zeros(2);
      temp4[0] = temp1[i], temp4[1] = temp2[i];
      if(min(abs(temp4))>threshold) {
        temp3[i] = temp3[i]/min(abs(temp4));
      }
    }
    if(max(abs(temp3))>maxtol) {
      plambda_hat.col(q+1) = llambdak(beta,gamma,plambda_hat.col(q),Shat,tk,data,qq);
      qlf += 1;
    } else {
      break;
    }
  }
  double pl = lloglik(beta,gamma,plambda_hat.col(qlf+1),Shat,tk,data,qq);
  
  return pl;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
colvec ssebgfun(colvec para, colvec lambda, mat Shat, colvec tk, List data, double maxtol, double maxiter, double threshold, double h, colvec qq) {
  int p = para.n_elem;
  mat covmat = zeros(p,p);
  double temp1 = pplfun(para,lambda,Shat,tk,data,maxtol,maxiter,threshold,qq);
  for(int i=0; i<p; i++) {
    colvec ei = zeros(p);
    ei[i] = 1;
    double temp2 = pplfun(para+h*ei,lambda,Shat,tk,data,maxtol,maxiter,threshold,qq);
    for(int j=0; j<p; j++) {
      colvec ej = zeros(p);
      ej[j] = 1;
      double temp3 = pplfun(para+h*ej,lambda,Shat,tk,data,maxtol,maxiter,threshold,qq);
      double temp4 = pplfun(para+h*ei+h*ej,lambda,Shat,tk,data,maxtol,maxiter,threshold,qq);
      covmat(i,j) = (temp1-temp2-temp3+temp4)/pow(h,2);
    }
  }
  mat covinv = -inv(covmat);
  colvec sebg = sqrt(covinv.diag());
  
  return sebg;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
colvec llambdaku(double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data, colvec qq, colvec u) {
  int n = Shat.n_rows, nt = Shat.n_cols;
  colvec delta1 = data["delta1"];
  colvec delta2 = data["delta2"];
  mat Z = data["Z"];
  colvec Rs = data["Rs"];
  colvec eta = data["eta"];
  mat W = WWmean(beta,gamma,lambda,Shat,tk,data);
  colvec gz = Z*gamma;
  colvec l = zeros(nt);
  colvec wt = (delta1 + delta2 + eta%(1-delta1-delta2)/qq)%u;
  for(int k=0; k<nt; k++) {
    colvec e = exp(gz+beta*Shat.col(k));
    double temp1 = 0, temp2 = 0;
    for(int i=0; i<n; i++) {
      bool index = (tk[k]<=Rs[i])&(is_finite(Shat(i,k)));
      if(index) {
        temp1 += wt[i]*W(i,k);
        temp2 += wt[i]*e[i];
      }
    }
    if(temp2>0) {
      l[k] = temp1/temp2;
    }
  }
  return l;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
colvec sscorebgu(double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data, colvec qq, colvec u) {
  int n = Shat.n_rows, nt = Shat.n_cols;
  int p = gamma.n_elem;
  colvec delta1 = data["delta1"];
  colvec delta2 = data["delta2"];
  colvec eta = data["eta"];
  mat Z = data["Z"];
  colvec Rs = data["Rs"];
  mat W = WWmean(beta,gamma,lambda,Shat,tk,data);
  colvec gz = Z*gamma;
  colvec lbg = zeros(p+1);
  colvec  wt = (delta1 + delta2 + eta%(1-delta1-delta2)/qq)%u;
  for(int k=0; k<nt; k++) {
    
    colvec e = exp(gz+beta*Shat.col(k));
    colvec eb = e%Shat.col(k);// exp(S_jk)*S_jk
    mat eg = e%Z.each_col();// whats means
    double temp1 = 0, temp2 = 0;
    colvec temp3 = zeros(p);
    
    for(int i=0; i<n; i++) {
      //wt[i] = delta1[i] + delta2[i] +  eta[i]*(1-delta1[i]-delta2[i])/0.3;
      bool index = (tk[k]<=Rs[i])&(is_finite(Shat(i,k)));
      if(index) {
        temp1 += wt[i]*e[i];
        temp2 += wt[i]*eb[i];
        temp3 += wt[i]*trans(eg.row(i));
      }
    }
    double tempb = 0;
    colvec tempg = zeros(p);
    for(int i=0; i<n; i++) {
      bool index = (tk[k]<=Rs[i])&(is_finite(Shat(i,k)));
      if(index) {
        //wt[i] = delta1[i] + delta2[i] +  eta[i]*(1-delta1[i]-delta2[i])/0.3;
        tempb += wt[i]*W(i,k)*(Shat(i,k)-temp2/temp1);
        tempg += wt[i]*W(i,k)*(trans(Z.row(i))-temp3/temp1);
      }
    }
    lbg[0] += tempb;
    lbg.subvec(1,p) += tempg;
  }
  return lbg;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat hhessianbgu(double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data, colvec qq, colvec u) {
  int n = Shat.n_rows, nt = Shat.n_cols;
  int p = gamma.n_elem;
  mat Z = data["Z"];
  colvec Rs = data["Rs"];
  colvec delta1 = data["delta1"];
  colvec delta2 = data["delta2"];
  colvec eta = data["eta"];
  colvec wt = (delta1 + delta2 + eta%(1-delta1-delta2)/qq)%u;
  mat W = WWmean(beta,gamma,lambda,Shat,tk,data);
  colvec gz = Z*gamma;
  mat lbg = zeros(p+1,p+1);
  
  for(int k=0; k<nt; k++) {
    colvec e = exp(gz+beta*Shat.col(k));// 
    colvec eb = e%Shat.col(k);// b represents for S, time dependent covariate, take 1st order derivative
    //for beta
    mat eg = e%Z.each_col();// g represents for gamma, time independent covariate, take 1st order derivative
    //for gamma 
    colvec ebb = eb%Shat.col(k); // second derivative for beta
    mat ebg = eb%Z.each_col(); // take derivative wrt both beta and gamma
    cube egg = zeros(p,p,n); // second derivative for gamma
    
    for(int i=0; i<n; i++) {
      egg.slice(i) = e[i]*trans(Z.row(i))*Z.row(i);// hessian for Z, how?
    }
    
    double temp1 = 0, temp2 = 0, temp4 = 0;
    rowvec temp3 = zeros<rowvec>(p), temp5 = zeros<rowvec>(p);
    mat temp6 = zeros(p,p);
    
    for(int i=0; i<n; i++) {
      bool index = (tk[k]<=Rs[i])&(is_finite(Shat(i,k)));
      if(index) {
        temp1 += wt[i]*e[i];
        temp2 += wt[i]*eb[i];
        temp3 += wt[i]*eg.row(i);
        temp4 += wt[i]*ebb[i];
        temp5 += wt[i]*ebg.row(i);
        temp6 += wt[i]*egg.slice(i); // add weight in this way is correct?
      }
    }
    
    double tempbb = 0;
    rowvec tempbg = zeros<rowvec>(p);
    mat tempgg = zeros(p,p);
    
    for(int i=0; i<n; i++) {
      bool index = (tk[k]<=Rs[i])&(is_finite(Shat(i,k)));
      if(index) {
        tempbb += wt[i]*W(i,k)*(-temp4/temp1-pow(temp2/temp1,2));
        tempbg += wt[i]*W(i,k)*(-temp5/temp1-temp2*temp3/pow(temp1,2));
        tempgg += wt[i]*W(i,k)*(-temp6/temp1-trans(temp3)*temp3/pow(temp1,2));
      }
    }
    lbg(0,0) += tempbb;
    lbg.submat(0,1,0,p) += tempbg;
    lbg.submat(1,1,p,p) += tempgg;
  }
  lbg.submat(1,0,p,0) = trans(lbg.submat(0,1,0,p));
  
  return lbg;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List eestbgu(double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data, double maxtol, double maxiter, double threshold, colvec qq, colvec u) {
  int p = gamma.n_elem;
  int nt = lambda.n_elem;
  
  colvec beta_hat = zeros(maxiter+2);
  mat gamma_hat = zeros(p,maxiter+2);
  mat lambda_hat = zeros(nt,maxiter+2);
  
  beta_hat[0] = beta-1, beta_hat[1] = beta;
  gamma_hat.col(0) = gamma-1, gamma_hat.col(1) = gamma;
  lambda_hat.col(0) = lambda-1, lambda_hat.col(1) = lambda;
  
  int qf = 0;
  for(int q=1; q<=maxiter; q++) {
    colvec temp1 = zeros(p+nt+1);
    temp1[0] = beta_hat[q];
    temp1.subvec(1,p) = gamma_hat.col(q);
    temp1.subvec(p+1,p+nt) = cumsum(lambda_hat.col(q));
    
    colvec temp2 = zeros(p+nt+1);
    temp2[0] = beta_hat[q-1];
    temp2.subvec(1,p) = gamma_hat.col(q-1);
    temp2.subvec(p+1,p+nt) = cumsum(lambda_hat.col(q-1));
    
    colvec temp3 = abs(temp1-temp2);
    for(int i=0; i<=(p+nt); i++) {
      colvec temp4 = zeros(2);
      temp4[0] = temp1[i], temp4[1] = temp2[i];
      if(min(abs(temp4))>threshold) {
        temp3[i] = temp3[i]/min(abs(temp4));
      }
    }// why
    
    if((max(temp3.subvec(0,p))+max(temp3.subvec(p+1,p+nt)))>maxtol) {
      
      lambda_hat.col(q+1) = llambdaku(beta_hat[q],gamma_hat.col(q),lambda_hat.col(q),Shat,tk,data,qq,u);
      
      colvec sbg = sscorebgu(beta_hat[q],gamma_hat.col(q),lambda_hat.col(q+1),Shat,tk,data,qq,u);
      
      mat hbg = hhessianbgu(beta_hat[q],gamma_hat.col(q),lambda_hat.col(q+1),Shat,tk,data,qq,u);
       
      colvec update = inv(hbg)*sbg;
      beta_hat[q+1] = beta_hat[q]-update[0];
      gamma_hat.col(q+1) = gamma_hat.col(q)-update.subvec(1,p);
      qf += 1;
    } else {
      break;
    }
    
  }
  List res;
  res["betaf"] = beta_hat[qf+1];
  res["gammaf"] = gamma_hat.col(qf+1);
  res["lambdaf"] = lambda_hat.col(qf+1);
  res["qf"] = qf;
  
  return res;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
colvec Bootstrap(int m,double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data, double maxtol, double maxiter,double threshold, colvec qq){
  int n = Shat.n_rows;
  int p = 1;
  //colvec nqb = zeros(100); does not consider unconverged case
  mat bggb_est = zeros(m, p+1);
  colvec see_est = zeros(p+1);
  for(int j=0; j<m; j++){
    //mat bggb_est = zeros(100,p);
    colvec u = Rcpp::rexp(n, 1.0);
    List ressb = eestbgu(beta,gamma,lambda,Shat,tk,data, maxtol, maxiter, threshold,qq, u);
    bggb_est(j,0) = ressb["betaf"];
    bggb_est(j,p) = ressb["gammaf"];
    //nqb[j] = ressb["qf""];
    //conv.flagb[j] = bool((nqb[j]>0)&(nqb[j]<max.iter));
   //cout<< j <<endl;
  }
  
  see_est[0] = arma::stddev(bggb_est.col(0));
  see_est[p] = arma::stddev(bggb_est.col(1));
  return see_est;
}


