% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ge_model_TGIE.R
\name{ge_model_TGIE}
\alias{ge_model_TGIE}
\title{implements a simple true-genotype-independent genotyping error model}
\usage{
ge_model_TGIE(L, epsilon = 0.01)
}
\arguments{
\item{L}{required locus specific information}

\item{epsilon}{the rate at which genotypes are incorrectly observed.}
}
\description{
implements a simple true-genotype-independent genotyping error model
}
\examples{
L <- list(
 freqs = c(A = 0.6, B = 0.4),
 geno_freqs = c(`A-A` = 0.36, `A-B` = 0.48, `B-B` = 0.16)
)
ge_model_TGIE(L, epsilon = 0.02)
}
