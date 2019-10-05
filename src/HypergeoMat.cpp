// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include "Dico.h"
using namespace std;
using namespace Rcpp;


Dico DictParts(int m, int n){
  unordered_map<int,int> D;
  arma::Mat<int> Last = {{0,m,m}};
  int fin = 0;
  for(int i = 0; i < n; i++){
    arma::Mat<int> NewLast = arma::Mat<int>(0,3);
    for(int j = 0; j < Last.n_rows; j++){
      int manque = Last(j,1);
      int l = min(manque, Last(j,2));
      if(l > 0){
        D[Last(j,0)] = fin + 1;
        for(int r = 0; r < l; r++){
          arma::Row<int> newrow = {fin+r+1, manque-r-1, r+1};
          NewLast = arma::join_cols(NewLast, newrow);
        }
        fin += l;
      }
    }
    Last = NewLast;
  }
  Dico out;
  out.dict = D; out.last = fin;
  return out;
}


IntegerVector cleanPart(IntegerVector kappa){
  return kappa[kappa>0];
}


NumericVector dualPartition(IntegerVector kappa, int to = -1){
  kappa = cleanPart(kappa);
  int l = kappa.size();
  if(l == 0){
    return NumericVector(0);
  }
  int l0 = to == -1 ? kappa(0) : to;
  NumericVector out = NumericVector(l0);
  out(0) = (double)l;
  for(int i = 1; i < l0 ; i++){
    int s = 0.0;
    for(int j = 0; j < l; j++){
      if(kappa(j) > i){
        s += 1.0;
      }
    }
    out(i) = s;
  }
  return out;
}


NumericVector sequence(int start, int end){
  int lout = end-start+1;
  NumericVector out = NumericVector(lout);
  for(int i = 0; i < lout; i++){
    out(i) = (double)(i+start);
  }
  return out;
}


double product(NumericVector v){
  double out = 1;
  for(int i = 0; i < v.size(); i++){
    out *= v(i);
  }
  return out;
}


double betaratio(IntegerVector kappa, IntegerVector mu, int k, double alpha){
  double t = (double)k - alpha * (double)mu(k-1);
  NumericVector v;
  if(k > 1){
    NumericVector mu_dbl = as<NumericVector>(mu);
    NumericVector ss = sequence(1, k-1);
    v = alpha * mu_dbl[Range(0,k-2)] - ss + t;
  }else{
    v = NumericVector(0);
  }
  NumericVector u;
  if(k > 0){
    NumericVector kappa_dbl = as<NumericVector>(kappa);
    NumericVector sss = sequence(1, k);
    u = alpha * kappa_dbl[Range(0,k-1)] - sss + t + 1.0;
  }else{
    u = NumericVector(0);
  }
  int l = mu(k-1) - 1;
  NumericVector w;
  if(l > 0){
    NumericVector muPrime = dualPartition(mu, l);
//    NumericVector muPrime = as<NumericVector>(muDual);
    NumericVector lrange = sequence(1,l);
    w = muPrime - alpha*lrange - t;
  }else{
    w = NumericVector(0);
  }
  double prod1 = product(u / (u + alpha - 1.0));
  double prod2 = product((v + alpha) / v);
  double prod3 = product((w + alpha) / w);
  return alpha * prod1 * prod2 * prod3;
}


double T_(double alpha, NumericVector a, NumericVector b, IntegerVector kappa){
  int lkappa = kappa.size();
  if(lkappa == 0 || kappa(0) == 0){
    return 1.0;
  }
  lkappa -= 1;
  int kappai = kappa(lkappa);
  double kappai_dbl = (double)kappai;
  double i = (double)lkappa;
  double c = kappai_dbl - 1.0 - i/alpha;
  double prod1_den = product(b + c);
  if(prod1_den == 0.0){
    return 0.0;
  }
  double d = kappai_dbl * alpha - i - 1.0;
  NumericVector e;
  if(kappai > 1){
    NumericVector s = sequence(1, kappai-1);
    NumericVector kappaPrime = dualPartition(kappa, kappai-1);
    e = kappaPrime - alpha*s + d;
  }else{
    e = NumericVector(0);
  }
  NumericVector g = e + 1.0;
  NumericVector f;
  if(lkappa > 0){
    NumericVector kappa_dbl = as<NumericVector>(kappa)[Range(0,lkappa-1)];
    NumericVector ss = sequence(1, lkappa);
    f = alpha*kappa_dbl - ss - d;
  }else{
    f = NumericVector(0);
  }
  NumericVector h = f + alpha;
  NumericVector l = h * f;
  double prod1_num = product(a + c);
  double prod2 = product((g - alpha) * e / g / (e + alpha));
  double prod3 = product((l-f) / (l+h));
  return prod1_num/prod1_den * prod2 * prod3;
}


void jack(double alpha, NumericVector x, unordered_map<int,int> dico, int k,
          double beta, int c, int t, IntegerVector mu, NumericMatrix& jarray,
          IntegerVector kappa, int nkappa){
  int i0 = k > 1 ? k : 1;
  int i1 = mu.size();
  for(int i = i0; i <= i1; i++){
    int u = mu(i-1);
    if(mu.size() == i || u > mu(i)){
      double gamma = beta * betaratio(kappa, mu, i, alpha);
      IntegerVector muP = clone(mu);
      muP(i-1) = u-1;
      muP = cleanPart(muP);
      int nmuP = 0;
      for(int j = 0; j < muP.size(); j++){
        nmuP = dico.at(nmuP) + muP(j) - 1;
      }
      if(muP.size() >= i && u > 1){
        jack(alpha, x, dico, i, gamma, c + 1, t, muP, jarray, kappa, nkappa);
      }else{
        if(nkappa > 1){
          if(muP.size() > 0){ // && muP(0)>0){
            jarray(nkappa-1, t-1) += gamma * jarray(nmuP-1, t-2) *
              pow(x(t-1),c+1);
          }else{
            jarray(nkappa-1, t-1) += gamma * pow(x(t-1),c+1);
          }
        }
      }
    }
  }
  if(k == 0){
    if(nkappa > 1){
      jarray(nkappa-1, t-1) += jarray(nkappa-1, t-2);
    }
  }else{
    int nmu = 0;
    for(int i = 0; i < mu.size(); i++){
      nmu = dico.at(nmu) + mu(i) - 1;
    }
    jarray(nkappa-1, t-1) += beta * pow(x(t-1),c) * jarray(nmu - 1, t-2);
  }
}


double summation(NumericVector a, NumericVector b, NumericVector x,
                 unordered_map<int,int> dico, int n, double alpha, int i,
                 double z, int j, IntegerVector kappa, NumericMatrix &jarray){
  if(i == n){
    return 0.0;
  }
  int lkappa = kappa.size();
  int lkappaP = lkappa + 1;
  int kappai = 1;
  double s = 0.0;
  while((i>0 || kappai<=j) &&
        (i==0 || ((lkappa==0 || kappai <= kappa(lkappa-1)) && kappai <= j))){
    IntegerVector kappaP(lkappa+1);
    for(int k=0; k < lkappa; k++){
      kappaP(k) = kappa(k);
    }
    kappaP(lkappa) = kappai;
    int nkappaP = 0;
    for(int k = 0; k < lkappaP; k++){
      nkappaP = dico.at(nkappaP) + kappaP(k) - 1;
    }
    z = z * T_(alpha, a, b, kappaP);
    if(nkappaP > 1 && (lkappaP == 1 || kappaP(1) == 0)){
      jarray(nkappaP-1, 0) = x(0) * (1.0 + alpha * (double)(kappaP(0)-1)) *
        jarray(nkappaP-2, 0);
    }
    for(int t = 2; t <= n; t++){
      jack(alpha, x, dico, 0, 1.0, 0, t, kappaP, jarray, kappaP, nkappaP);
    }
    s += z * jarray(nkappaP-1, n-1);
    if(j > kappai && i <= n){
      s += summation(a, b, x, dico, n, alpha, i+1, z, j-kappai, kappaP, jarray);
    }
    kappai += 1;
  }
  return s;
}


double summationI(NumericVector a, NumericVector b, double x, int n,
                  double alpha, int i, double z, int j, IntegerVector kappa){
  int lkappa = kappa.size();
  int kappai = 1;
  double s = 0.0;
  while((i>0 || kappai<=j) && (i==0 || kappai <= kappa(i-1) && kappai <= j)){
    IntegerVector kappaP(lkappa+1);
    for(int k=0; k < lkappa; k++){
      kappaP(k) = kappa(k);
    }
    kappaP(lkappa) = kappai;
    double t = T_(alpha, a, b, kappaP);
    z = z * x * t * ((double)(n-i) + alpha * (double)(kappai-1));
    if(j > kappai && i <= n){
      s += summationI(a, b, x, n, alpha, (i+1), z, (j-kappai), kappaP);
    }
    s += z;
    kappai += 1;
  }
  return s;
}


double hypergeoI(int m, double alpha, NumericVector a, NumericVector b,
                 int n, double x){
  return 1.0 + summationI(a, b, x, n, alpha, 0, 1.0, m, IntegerVector(0));
}


// [[Rcpp::export]]
double Rcpp_hypergeomPFQ(int m, NumericVector a, NumericVector b,
                         NumericVector x, double alpha){
  if(is_true(all(x == x(0)))){
    return hypergeoI(m, alpha, a, b, x.size(), x(0));
  }
  int n = x.size();
  Dico dict = DictParts(m, n);
  NumericMatrix jarray(dict.last, n);
  NumericVector xx = cumsum(x);
  for(int j = 0; j < n; j++){
    jarray(0,j) = xx(j);
  }
  IntegerVector emptyPart = IntegerVector(0);
  double s =
    summation(a, b, x, dict.dict, n, alpha, 0, 1.0, m, emptyPart, jarray);
  return 1.0 + s;
}

