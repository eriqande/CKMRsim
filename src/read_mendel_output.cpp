#include <Rcpp.h>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace Rcpp;


//' pick the genotypes out of the Mendel output pedigree file to use to compute Q
//'
//' Not sure exactly how I am going to do this, as we need to include genotyping error on there as well.
//' Crap.  But at least I have a start on how to parse this monstrous file. This assumes that the focal pair
//' of individuals are labeled 1 and 2.
//' @param Input the path to the Mendel output file to read in.
//' @param NumA the number of alleles at each locus
//' @examples
//' read_mendel_outped("/Users/eriq/Desktop/mendel-example-Ped.out")
// [[Rcpp::export]]
List read_mendel_outped(CharacterVector Input, IntegerVector NumA) {
  int i,j, slash;
  int aa,bb, a, b, A;
  std::string tempstr;
  std::string line;
  std::string word;
  std::string fname = Rcpp::as<std::string>(Input);
  int NumLoc = NumA.size();
  int geno;
  int Num1 = 0, Num2 = 0;
  int indiv;


  std::ifstream infile (fname);

  /* first we scan through the file and count lines from individuals "1" and "2" */
  while (std::getline(infile, line)) {
    std::stringstream ss(line);

    for(j=0;j<3;j++) ss >> word;  // get to the indiv identifier on the line.  This is screwed up because at reps < 1000 there is a space before the comma...
    if(word == "1") Num1++;
    if(word == "2") Num2++;
  }
  if(Num1 != Num2) stop("Different numbers of 1's and 2's in mendel output.");

  /* then allocate memory for return genotype values and reset the file */
  Rcpp::Rcout << "Have scanned " << Num1 << " pairs in Mendel output file." << std::endl;
  IntegerMatrix geno1(Num1, NumLoc);  // Return as matrix with Num1 rows and NumLoc columns.
  IntegerMatrix geno2(Num2, NumLoc);
  Num1 = -1;  // reset these to capture data into the matrices
  Num2 = -1;

  infile.clear();
  infile.seekg(0, infile.beg);

  while (std::getline(infile, line)) {
    std::stringstream ss(line);

    for(j=0;j<3;j++) ss >> word;  // get to the indiv identifier on the line
    if(word == "1" || word == "2") {
      if(word == "1") {
        Num1++;
        indiv = 1;
      }
      if(word == "2") {
        Num2++;
        indiv = 2;
      }

      //Rcpp::Rcout << word << "  ";
      for(i=0;i<7;i++) ss >> word;  // skip over all the other columns to get to the genotypes
      for(i=0;i<NumLoc;i++) {
        ss >> word;  // read the genotype
        slash = word.find("/");
        aa = std::stoi(word.substr(1, slash - 1));
        bb = std::stoi(word.substr(slash + 1));
        if(aa <= bb) {
          a = aa;
          b = bb;
        } else {
          a = bb;
          b = aa;
        }
        A = NumA[i];
        geno = 2 + (a - 1) * (A + 2) - ( a * (a + 1) / 2) + (b - a);
        if(indiv == 1) geno1(Num1, i) = geno;
        if(indiv == 2) geno2(Num2, i)  = geno;

        //Rcpp::Rcout << word << " " << a << "-" << b << "  " << geno << "      ";
      }
      //Rcpp::Rcout << "\n";
    }
  }
  List ret;
  ret["indiv1"] = geno1;
  ret["indiv2"] = geno2;
  return(ret);
}
