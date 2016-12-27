#include <Rcpp.h>
using namespace Rcpp;

//' Compute pairwise relationship measures between all individuals in source and one individual in target
//'
//' More explanation later.
//'
//' @param S "source", a matrix whose rows are integers, with NumInd-source rows and NumLoci columns, with each
//' entry being a a base-0 representation of the genotype of the c-th locus at the r-th individual.
//' These are the individuals you can think of as parents if there is directionality to the
//' comparisons.
//' @param T "target",  a matrix whose rows are integers, with NumInd-target rows and NumLoci columns, with each
//' entry being a a base-0 representation of the genotype of the c-th locus at the r-th individual.
//' These are the individuals you can think of as offspring if there is directionality to the
//' comparisons.
//' @param t the index (base-1) of the individual in T that you want to compare against
//' everyone on S.
//' @param values the vector of genotype specific values.  See the probs field of \code{\link{flatten_ckmr}}.
//' @param nGenos a vector of the number of genotypes at each locus
//' @param base0_locus_starts the base0 indexes of the starting positions of each locus in probs.
//'
//' @return a data frame with columns "ind" (the base-1 index of the individual in S),
//' "value" (the value extracted, typically a log likelihood ratio), and "num_loc" (the
//' number of non-missing loci in the comparison.)
//' @export
// [[Rcpp::export]]
DataFrame comp_ind_pairwise(IntegerMatrix S, IntegerMatrix T, int t, NumericVector values, IntegerVector nGenos, IntegerVector Starts) {
  int L = S.ncol();
  int nS = S.nrow();
  int i,j,n;
  int tG, sG;
  double sum;
  IntegerVector popi(nS);
  NumericVector val(nS);
  IntegerVector nonmiss(nS);  // for the number of non_missing loci


  for(i=0;i<nS;i++) {
    sum = 0.0;
    n = 0;
    for(j=0;j<L;j++) {
      tG = T(t - 1, j);
      sG = S(i, j);
      if(tG >=0 && sG >= 0) {  // if both individuals have non-missing data
        n++;
        sum += values[Starts[j] + nGenos[j] * sG + tG];
      }
    }
    popi[i] = i + 1;
    val[i] = sum;
    nonmiss[i] = n;
  }

  return(DataFrame::create(Named("ind") = popi,
                           Named("value") = val,
                           Named("num_loc") = nonmiss));
}
