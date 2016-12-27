#include <Rcpp.h>
using namespace Rcpp;


//' Sample 1 observation from cell probabilities that are columns of a matrix
//'
//' Takes a matrix in which rows sum to one. For each row, performs a
//' single multinomial draw from amongst the columns, weighted by their values in that column
//'
//' @param M a matrix whose rows are reals, each summing to one
//'
//' @return a vector length = \code{nrow(M)} of indices, with each element being
//' the column that was chosen in that row's sampling
// [[Rcpp::export]]
IntegerVector samp_from_mat(NumericMatrix M) {
  int C = M.ncol();
  int R = M.nrow();
  int r,c, res;
  double cumul, rando;
  IntegerVector out(R);
  NumericVector rando_vec(R);

  rando_vec = runif(R);

  for(r = 0; r < R; r++) {
    cumul = 0.0;
    rando = rando_vec[r];
    for(c = 0; c < C; c++) {
      cumul += M(r, c);
      res = c + 1L;
      if(cumul >= rando) {
        break;
      }
    }
    out[r] = res;
  }
  return(out);
}
