#include <Rcpp.h>
#include <algorithm>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

namespace {

int geno_index(int a, int b, int n_alleles) {
  if(a > b) {
    std::swap(a, b);
  }
  return 2 + (a - 1) * (n_alleles + 2) - (a * (a + 1) / 2) + (b - a);
}

int sample_allele(NumericVector probs) {
  double u = R::runif(0.0, 1.0);
  double cumul = 0.0;
  int n = probs.size();

  for(int i = 0; i < n; i++) {
    cumul += probs[i];
    if(cumul >= u) {
      return i + 1;
    }
  }

  return n;
}

int offset(int indiv, int homolog, int locus, int n_indiv, int n_loci) {
  return locus + n_loci * (homolog + 2 * indiv);
}

class GeneDropper {
public:
  GeneDropper(
    IntegerVector kid_,
    IntegerVector pa_,
    IntegerVector ma_,
    IntegerVector n_alleles_,
    NumericVector map_pos_,
    IntegerVector chrom_id_,
    IntegerVector chrom_start_,
    IntegerVector chrom_end_,
    List freqs_,
    int min_crossovers_
  ) : kid(kid_),
      pa(pa_),
      ma(ma_),
      n_alleles(n_alleles_),
      map_pos(map_pos_),
      chrom_id(chrom_id_),
      chrom_start(chrom_start_),
      chrom_end(chrom_end_),
      freqs(freqs_),
      min_crossovers(min_crossovers_),
      n_indiv(kid_.size()),
      n_loci(n_alleles_.size()),
      hap(n_indiv * 2 * n_loci),
      done(n_indiv) {

    for(int i = 0; i < n_indiv; i++) {
      row_by_id[kid[i]] = i;
    }
  }

  void reset() {
    std::fill(hap.begin(), hap.end(), 0);
    std::fill(done.begin(), done.end(), false);
  }

  void ensure(int indiv) {
    if(done[indiv]) {
      return;
    }

    int father = pa[indiv];
    int mother = ma[indiv];

    if(father == 0 && mother == 0) {
      sample_founder(indiv);
      done[indiv] = true;
      return;
    }

    if(father == 0 || mother == 0) {
      stop("Every non-founder in the pedigree must have both parents specified.");
    }

    int father_row = parent_row(father);
    int mother_row = parent_row(mother);

    ensure(father_row);
    ensure(mother_row);

    transmit(father_row, indiv, 0);
    transmit(mother_row, indiv, 1);

    done[indiv] = true;
  }

  int genotype(int indiv, int locus) const {
    int a = hap[offset(indiv, 0, locus, n_indiv, n_loci)];
    int b = hap[offset(indiv, 1, locus, n_indiv, n_loci)];
    return geno_index(a, b, n_alleles[locus]);
  }

private:
  IntegerVector kid;
  IntegerVector pa;
  IntegerVector ma;
  IntegerVector n_alleles;
  NumericVector map_pos;
  IntegerVector chrom_id;
  IntegerVector chrom_start;
  IntegerVector chrom_end;
  List freqs;
  int min_crossovers;
  int n_indiv;
  int n_loci;
  std::vector<int> hap;
  std::vector<bool> done;
  std::unordered_map<int, int> row_by_id;

  int parent_row(int parent_id) const {
    auto it = row_by_id.find(parent_id);
    if(it == row_by_id.end()) {
      stop("Parent ID not found in pedigree.");
    }
    return it->second;
  }

  void sample_founder(int indiv) {
    for(int locus = 0; locus < n_loci; locus++) {
      NumericVector p = freqs[locus];
      hap[offset(indiv, 0, locus, n_indiv, n_loci)] = sample_allele(p);
      hap[offset(indiv, 1, locus, n_indiv, n_loci)] = sample_allele(p);
    }
  }

  void transmit(int parent, int child, int child_homolog) {
    int n_chrom = chrom_start.size();

    for(int chrom = 0; chrom < n_chrom; chrom++) {
      int start = chrom_start[chrom] - 1;
      int end = chrom_end[chrom] - 1;
      double left = map_pos[start];
      double right = map_pos[end];
      double span = right - left;
      int current_homolog = R::runif(0.0, 1.0) < 0.5 ? 0 : 1;
      std::vector<double> crossovers;

      if(span > 0.0) {
        int n_cross = R::rpois(span / 100.0);
        if(n_cross < min_crossovers) {
          n_cross = min_crossovers;
        }
        crossovers.reserve(n_cross);
        for(int i = 0; i < n_cross; i++) {
          crossovers.push_back(R::runif(left, right));
        }
        std::sort(crossovers.begin(), crossovers.end());
      }

      int next_crossover = 0;
      for(int locus = start; locus <= end; locus++) {
        while(next_crossover < static_cast<int>(crossovers.size()) &&
              crossovers[next_crossover] < map_pos[locus]) {
          current_homolog = 1 - current_homolog;
          next_crossover++;
        }

        hap[offset(child, child_homolog, locus, n_indiv, n_loci)] =
          hap[offset(parent, current_homolog, locus, n_indiv, n_loci)];
      }
    }
  }
};

} // namespace

//' Gene-drop true linked genotypes for the two focal individuals in a pedigree
//'
//' @param kid integer vector of pedigree IDs.
//' @param pa integer vector of paternal IDs, with 0 for founders.
//' @param ma integer vector of maternal IDs, with 0 for founders.
//' @param n_alleles integer vector with the number of alleles at each locus.
//' @param map_pos numeric vector of map positions in centiMorgans.
//' @param chrom_id integer chromosome ID per locus.
//' @param chrom_start one-based start locus index for each chromosome.
//' @param chrom_end one-based end locus index for each chromosome.
//' @param freqs list of allele-frequency vectors, one per locus.
//' @param reps number of replicate gene drops.
//' @param min_crossovers minimum number of crossovers in the typed marker span
//' of each chromosome for each meiosis.
//' @return A list with integer matrices \code{indiv1} and \code{indiv2}.
// [[Rcpp::export]]
List linked_gene_drop_true_genotypes(
    IntegerVector kid,
    IntegerVector pa,
    IntegerVector ma,
    IntegerVector n_alleles,
    NumericVector map_pos,
    IntegerVector chrom_id,
    IntegerVector chrom_start,
    IntegerVector chrom_end,
    List freqs,
    int reps,
    int min_crossovers = 0
) {
  if(kid.size() != pa.size() || kid.size() != ma.size()) {
    stop("kid, pa, and ma must have the same length.");
  }
  if(n_alleles.size() != map_pos.size() || n_alleles.size() != chrom_id.size()) {
    stop("n_alleles, map_pos, and chrom_id must have the same length.");
  }
  if(freqs.size() != n_alleles.size()) {
    stop("freqs must have one component per locus.");
  }
  if(reps < 1) {
    stop("reps must be positive.");
  }
  if(min_crossovers < 0) {
    stop("min_crossovers must be non-negative.");
  }

  int n_loci = n_alleles.size();
  IntegerMatrix indiv1(reps, n_loci);
  IntegerMatrix indiv2(reps, n_loci);

  std::unordered_map<int, int> row_by_id;
  for(int i = 0; i < kid.size(); i++) {
    row_by_id[kid[i]] = i;
  }
  if(row_by_id.find(1) == row_by_id.end() || row_by_id.find(2) == row_by_id.end()) {
    stop("Pedigree must include focal individuals with Kid IDs 1 and 2.");
  }
  int indiv1_row = row_by_id[1];
  int indiv2_row = row_by_id[2];

  GeneDropper dropper(
    kid, pa, ma, n_alleles, map_pos, chrom_id, chrom_start, chrom_end,
    freqs, min_crossovers
  );

  for(int rep = 0; rep < reps; rep++) {
    dropper.reset();
    dropper.ensure(indiv1_row);
    dropper.ensure(indiv2_row);

    for(int locus = 0; locus < n_loci; locus++) {
      indiv1(rep, locus) = dropper.genotype(indiv1_row, locus);
      indiv2(rep, locus) = dropper.genotype(indiv2_row, locus);
    }
  }

  return List::create(
    _["indiv1"] = indiv1,
    _["indiv2"] = indiv2
  );
}
