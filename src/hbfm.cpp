#include <Rcpp.h>
using namespace Rcpp;

// C++ code for speeding up calculations in hbfm
//

// [[Rcpp::export]]
NumericMatrix lPOWa(NumericMatrix a, NumericMatrix tl, NumericVector p) {
  int ng = a.nrow();  //number of genes
  int nc = tl.ncol(); //number of cells
  int nf = a.ncol();  //number of factors
  
  NumericMatrix xnew(nc, ng);
  
  for (int i = 0; i < ng; i++){
    for (int j = 0; j < nc; j++){
      xnew(j,i) = 1;
      for(int k = 0; k < nf; k++){
        xnew(j,i) *= exp(-1*p[k]/2*abs(a(i,k))) * pow(tl(k,j),a(i,k));
      }
    }
  }
  return xnew;
}



// [[Rcpp::export]]
NumericVector lVECa(NumericMatrix a, NumericVector tl, NumericVector p) {
  int ng = a.ncol();  //number of genes
  int nf = a.nrow();  //number of factors
  
  NumericVector xnew(ng);
  
  for (int i = 0; i < ng; i++){
    xnew(i) = 1;
    for(int k = 0; k < nf; k++){
      xnew(i) *= exp(-1*p[k]/2*abs(a(k,i))) * pow(tl(k),a(k,i));
    }
  }
  return xnew;
}



// [[Rcpp::export]]
NumericMatrix margLik(int nh, NumericMatrix a, NumericVector b, NumericVector p, NumericMatrix Y) {
  int ng = Y.nrow();  //number of genes
  int nc = Y.ncol(); //number of cells
  int nf = a.ncol(); //number of factors
  
  NumericMatrix xlam(nf, nc);
  NumericMatrix xmean(ng, nc);
  double xmu;
  //double xml[ng][nc];
  
  for (int h = 0; h < nh; h++){
    //sample lambdas
    for (int i = 0; i < nc; i++){
      for (int k = 0; k < nf; k++){
        xlam(k,i) = R::rlnorm(0.0, sqrt(p[k]));
      }
    }
    
    //calculate log marginal likelihood
    NumericMatrix xnew = lPOWa(a, xlam, p); //nc by ng matrix
    for (int g = 0; g < ng; g++){
      for (int i = 0; i < nc; i++){
        xmu = b[g]*xnew(i,g);
        xmean(g,i) += R::dpois(Y(g,i), xmu, true)/nh;
      }//end i loop
    }//end g loop
  }//end h loop
  
  return xmean;
}


