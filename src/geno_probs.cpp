#include <Rcpp.h>
using namespace Rcpp;


//' compute the matrix X_l given allele frequencies and kappa
//'
//' Ha! I did this initially using R straight up and even tried to make
//' it reasonably vectorized, but the damn thing took forever and was
//' incredible space-inefficient. Ergo, we are going after it using
//' Rcpp.  The genotypes are ordered in rows and columns in the way
//' as described in the paper.
//' @param p vector of allele frequencies
//' @param kappa a 3-vector of the Cotterman coefficients
// [[Rcpp::export]]
NumericMatrix make_matrix_X_l(NumericVector p, NumericVector kappa) {
  int nA = p.size(); // number of alleles
  int a,b,c,d;  // allele indexes of the first (a) and second (b) alleles at individual 1, and (c) and (d) are similar for individual 2
  int i, j;  // indexes of the first and second individual's genotype
  int nG = nA * (nA + 1) / 2.0;  // number of genotypes
  double k0 = kappa[0];
  double k1 = kappa[1];
  double k2 = kappa[2];
  double gp1, gp2;  // for marginal genotype probs of individual 1 and 2
  double tp;  // for the PO transmission prob
  NumericMatrix ret(nG, nG); // allocate memory to be returned

  for(a = 0, i=0; a < nA; a++) {
    for(b = a; b < nA; b++, i++) {
      for(c = 0, j = 0; c < nA; c++) {
        for(d = c; d < nA; d++, j++) {
          // compute the marginal genotype probs
          gp1 = p[a] * p[b] * (2 - (a == b));
          gp2 = p[c] * p[d] * (2 - (c == d));

          // I could write this more quickly with arithmetic like (a==b) * ((c==d) * ...)
          // but that would be harder to read and understand when I came back to it.  I
          // don't think there will be any speed diff.  Note that for this calculation I do
          // not assume that a<=b and c<=d.
          if(a==b)  { // for all cases in which indiv 1 is homozygous
            if(c==d) {
              if(c==a) {
                tp = p[a];
              }
              else {
                tp = 0.0;
              }
            }
            else {  // for the case when indiv 2 is heterozygous
              if(c!=a && d!=a) {  // indiv 2 has no alleles in common with indiv 1
                tp = 0.0;
              }
              else if (c==a) {  // by this point, indiv2 must have exactly one shared allele with indiv 1, so it is one or the other here
                tp = p[d];
              }
              else {
                tp = p[c];
              }
            }
          }
          else {   // all cases in which indiv 1 is heterozygous
            if(d==c) {  // if indiv 2 is homozygous
              if(a!=c && b!=c) {  // they share no alleles
                tp = 0.0;
              }
              else if(a==c) {  // then they must share exactly 1 allele
                tp = p[a] / 2.0;
              }
              else {
                tp = p[b] / 2.0;
              }
            }
            else {  // indiv 2 is heterozygous
              if( (a==c && b==d) || (a==d && b==c)) {  // indiv 1 and 2 have the same heterozygous genotype
                tp = (p[a] + p[b]) / 2.0;
              }
              else if(c==a || c==b) {  // allele c in indiv 2 is the one that is shared
                tp = p[d] / 2.0;
              }
              else if(d==a || d==b) { // allele d in indiv 2 is the one that is shared
                tp = p[c] / 2.0;
              }
              else {  // in this case, no alleles are shared
                tp = 0.0;
              }
            }
          }
          ret(i,j) = k0 * gp1 * gp2  +  k1 * gp1 * tp   +   k2 * gp1 * (i==j);
        }
      }
    }
  }

  return(ret);
}
