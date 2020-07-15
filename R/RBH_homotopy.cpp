#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int Backtrack(IntegerVector RCV,
              int end) {
  // Compute the position of the last non-negative entries in the rejection-corrected vector up to location "end". -1 outputed if all entries are negative. 
  // Inputs: 
  //   RCV: rejection-corrected vector. k-th element #{i: p_i <= alpha a_k / n A} - a_k
  //   end: the last position to look at. It is not necessarily the length of RCV
  // 
  int ind;

  for (int i = 1; i <= end; ++i) {
    ind = end - i;
    if (RCV[ind] >= 0){
      return ind;
    }
  }
  return -1;
}

// [[Rcpp::export]]
IntegerVector RejsBH(IntegerVector posit,
                     IntegerVector sgn, 
                     IntegerVector RCV,
                     IntegerVector avals) {
  int m = posit.size();
  int n = RCV.size();
  int nrejs;
  int ind;
  int ind0;
  IntegerVector RBH(m + 1);

  ind0 = Backtrack(RCV, n);
  if (ind0 < 0) {
    nrejs = 0;
  } else {
    nrejs = avals[ind0] + RCV[ind0];
  }
  RBH[0] = nrejs;
  for (int i = 0; i < m; ++i) {
    ind = posit[i];
    RCV[ind] += sgn[i];
    if (ind > ind0) {
      if (RCV[ind] == 0){
        nrejs = avals[ind];
        ind0 = ind;
      }
    } else if (ind == ind0) {
      if (RCV[ind] < 0) {
        ind0 = Backtrack(RCV, ind);
        nrejs = avals[ind0] + RCV[ind0];
      } else {
        nrejs += sgn[i];
      }
    }
    RBH[i + 1] = nrejs;
  }
  return RBH;
}
