#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat Wmean(double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data) {
  int n = Shat.n_rows, nt = Shat.n_cols;
  mat Z = data["Z"];
  colvec delta = data["delta"]; // censoring indicator
  colvec T = data["T"];// exact happening time
  colvec L = data["L"];// interval censoring left
  colvec R = data["R"];//interval censoring right
  colvec gz = Z*gamma; // ? covariate times coefficient
  mat res = zeros(n,nt);
  for(int i=0; i<n; i++) {
    if(delta[i]==1) {
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

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
colvec lambdak(double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data) {
  int n = Shat.n_rows, nt = Shat.n_cols;
  mat Z = data["Z"];
  colvec Rs = data["Rs"];
  mat W = Wmean(beta,gamma,lambda,Shat,tk,data);
  colvec gz = Z*gamma;
  colvec l = zeros(nt);
  for(int k=0; k<nt; k++) {
    colvec e = exp(gz+beta*Shat.col(k));
    double temp1 = 0, temp2 = 0;
    for(int i=0; i<n; i++) {
      bool index = (tk[k]<=Rs[i])&(is_finite(Shat(i,k)));
      if(index) {
        temp1 += W(i,k);
        temp2 += e[i];
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
colvec scorebg(double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data) {
  int n = Shat.n_rows, nt = Shat.n_cols;
  int p = gamma.n_elem;
  mat Z = data["Z"];
  colvec Rs = data["Rs"];
  mat W = Wmean(beta,gamma,lambda,Shat,tk,data);
  colvec gz = Z*gamma;
  colvec lbg = zeros(p+1);
  for(int k=0; k<nt; k++) {
    colvec e = exp(gz+beta*Shat.col(k));
    colvec eb = e%Shat.col(k);
    mat eg = e%Z.each_col();
    double temp1 = 0, temp2 = 0;
    colvec temp3 = zeros(p);
    for(int i=0; i<n; i++) {
      bool index = (tk[k]<=Rs[i])&(is_finite(Shat(i,k)));
      if(index) {
        temp1 += e[i];
        temp2 += eb[i];
        temp3 += trans(eg.row(i));
      }
    }
    double tempb = 0;
    colvec tempg = zeros(p);
    for(int i=0; i<n; i++) {
      bool index = (tk[k]<=Rs[i])&(is_finite(Shat(i,k)));
      if(index) {
        tempb += W(i,k)*(Shat(i,k)-temp2/temp1);
        tempg += W(i,k)*(trans(Z.row(i))-temp3/temp1);
      }
    }
    lbg[0] += tempb;
    lbg.subvec(1,p) += tempg;
  }
  return lbg;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat hessianbg(double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data) {
  int n = Shat.n_rows, nt = Shat.n_cols;
  int p = gamma.n_elem;
  mat Z = data["Z"];
  colvec Rs = data["Rs"];
  mat W = Wmean(beta,gamma,lambda,Shat,tk,data);
  colvec gz = Z*gamma;
  mat lbg = zeros(p+1,p+1);
  for(int k=0; k<nt; k++) {
    colvec e = exp(gz+beta*Shat.col(k));
    colvec eb = e%Shat.col(k);
    mat eg = e%Z.each_col();
    colvec ebb = eb%Shat.col(k);
    mat ebg = eb%Z.each_col();
    cube egg = zeros(p,p,n);
    for(int i=0; i<n; i++) {
      egg.slice(i) = e[i]*trans(Z.row(i))*Z.row(i);
    }
    double temp1 = 0, temp2 = 0, temp4 = 0;
    rowvec temp3 = zeros<rowvec>(p), temp5 = zeros<rowvec>(p);
    mat temp6 = zeros(p,p);
    for(int i=0; i<n; i++) {
      bool index = (tk[k]<=Rs[i])&(is_finite(Shat(i,k)));
      if(index) {
        temp1 += e[i];
        temp2 += eb[i];
        temp3 += eg.row(i);
        temp4 += ebb[i];
        temp5 += ebg.row(i);
        temp6 += egg.slice(i);
      }
    }
    double tempbb = 0;
    rowvec tempbg = zeros<rowvec>(p);
    mat tempgg = zeros(p,p);
    for(int i=0; i<n; i++) {
      bool index = (tk[k]<=Rs[i])&(is_finite(Shat(i,k)));
      if(index) {
        tempbb += W(i,k)*(-temp4/temp1-pow(temp2/temp1,2));
        tempbg += W(i,k)*(-temp5/temp1-temp2*temp3/pow(temp1,2));
        tempgg += W(i,k)*(-temp6/temp1-trans(temp3)*temp3/pow(temp1,2));
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
List estbg(double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data, double maxtol, double maxiter, double threshold) {
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
    }
    if((max(temp3.subvec(0,p))+max(temp3.subvec(p+1,p+nt)))>maxtol) {
      lambda_hat.col(q+1) = lambdak(beta_hat[q],gamma_hat.col(q),lambda_hat.col(q),Shat,tk,data);
      colvec sbg = scorebg(beta_hat[q],gamma_hat.col(q),lambda_hat.col(q+1),Shat,tk,data);
      mat hbg = hessianbg(beta_hat[q],gamma_hat.col(q),lambda_hat.col(q+1),Shat,tk,data);
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
double loglik(double beta, colvec gamma, colvec lambda, mat Shat, colvec tk, List data) {
  int n = Shat.n_rows, nt = Shat.n_cols;
  mat Z = data["Z"];
  colvec delta = data["delta"];
  colvec T = data["T"];
  colvec L = data["L"];
  colvec R = data["R"];
  colvec gz = Z*gamma;
  double logl = 0;
  for(int i=0; i<n; i++) {
    colvec le = lambda%trans(exp(gz[i]+beta*Shat.row(i)));
    if(delta[i]==1) {
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
        logl += log(tempt2)-tempt1;
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
      if(!is_finite(R[i])){
        logl += log(ll);
      } else if(cr>0){
        logl += log(ll-lr);
      }
    }
  }
  return logl;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double plfun(colvec para, colvec lambda, mat Shat, colvec tk, List data, double maxtol, double maxiter, double threshold) {
  int p = para.n_elem-1;
  double beta = para[0];
  colvec gamma = para.subvec(1,p);
  int nt = lambda.n_elem;
  mat plambda_hat = zeros(nt,maxiter+2);
  plambda_hat.col(0) = lambda-1, plambda_hat.col(1) = lambda;
  int qlf = 0;
  for(int q=1; q<=maxiter; q++) {
    colvec temp1 = cumsum(plambda_hat.col(q));
    colvec temp2 = cumsum(plambda_hat.col(q-1));
    colvec temp3 = abs(temp1-temp2);
    for(int i=0; i<nt; i++) {
      colvec temp4 = zeros(2);
      temp4[0] = temp1[i], temp4[1] = temp2[i];
      if(min(abs(temp4))>threshold) {
        temp3[i] = temp3[i]/min(abs(temp4));
      }
    }
    if(max(abs(temp3))>maxtol) {
      plambda_hat.col(q+1) = lambdak(beta,gamma,plambda_hat.col(q),Shat,tk,data);
      qlf += 1;
    } else {
      break;
    }
  }
  double pl = loglik(beta,gamma,plambda_hat.col(qlf+1),Shat,tk,data);
  
  return pl;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
colvec sebgfun(colvec para, colvec lambda, mat Shat, colvec tk, List data, double maxtol, double maxiter, double threshold, double h) {
  int p = para.n_elem;
  mat covmat = zeros(p,p);
  double temp1 = plfun(para,lambda,Shat,tk,data,maxtol,maxiter,threshold);
  for(int i=0; i<p; i++) {
    colvec ei = zeros(p);
    ei[i] = 1;
    double temp2 = plfun(para+h*ei,lambda,Shat,tk,data,maxtol,maxiter,threshold);
    for(int j=0; j<p; j++) {
      colvec ej = zeros(p);
      ej[j] = 1;
      double temp3 = plfun(para+h*ej,lambda,Shat,tk,data,maxtol,maxiter,threshold);
      double temp4 = plfun(para+h*ei+h*ej,lambda,Shat,tk,data,maxtol,maxiter,threshold);
      covmat(i,j) = (temp1-temp2-temp3+temp4)/pow(h,2);
    }
  }
  mat covinv = -inv(covmat);
  colvec sebg = sqrt(covinv.diag());
  
  return sebg;
}

