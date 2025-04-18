% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_close_matching_genotypes.R
\name{find_close_matching_genotypes}
\alias{find_close_matching_genotypes}
\title{Return every pair of individuals that mismatch at no more than max_miss loci}
\usage{
find_close_matching_genotypes(LG, CK, max_mismatch)
}
\arguments{
\item{LG}{a long genotypes data frame.}

\item{CK}{a ckmr object created from the allele frequencies computed from LG.}

\item{max_mismatch}{maximum allowable number of mismatching genotypes betwen the pairs.}
}
\value{
a data frame with columns:
\describe{
\item{indiv_1}{the id (from the rownames in S) of the firt member of the pair}
\item{indiv_2}{the id (from the rownames in S) of the second individual of the pair}
\item{num_mismatch}{the number of loci at which the pair have mismatching genotypes}
\item{num_loc}{the total number of loci that are non-missing in both individuals of the pair.}
}
}
\description{
This is used for identifying duplicate individuals/genotypes in large
data sets. I've specified this in terms of the max number of missing loci because
I think everyone should already have tossed out individuals with a lot of
missing data, and then it makes it easy to toss out pairs without even
looking at all the loci, so it is faster for all the comparisons.
}
