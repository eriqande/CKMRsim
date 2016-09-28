



#' A list of pedigrees specifying different pairwise relationships for MORGAN's genedrop program.
#'
#' A list of pedigrees, each specified as a data frame with five columns.
#'
#'
#' @format A named list of pedigrees, each one specified in a data frame with five columns in the following order:
#' \describe{
#'   \item{Kid}{alphanumeric character specifier of the kid's name.  Must be <15 characters.}
#'   \item{Pa}{alphanumeric character specifier of the dad's name.  Must be <15 characters.}
#'   \item{Ma}{alphanumeric character specifier of the mom's name.  Must be <15 characters.}
#'   \item{Sex}{Integer specifier of sex of the individual: 0 - unknown, 1 - male, 2 - female.}
#'   \item{Observed}{Integer specifying if the individual is observed.  Typically there will be two observed
#'                   individuals.  Those should probably be 1 and 2. }
#' }
"pedigrees"



#' A matrix of Cotterman coefficients for a variety of pairwise relationships
#'
#' A matrix with named rows (relationships) and named columns (the kappas).
#' Eventually each row should go into a list along with the pedigree information
#' as well.  But for now I wanted to get it down somewhere.
#'
#'
#' @format A matrix.  The rows are relationships (see the CKMRsim paper for abbreviations).
#' The columns are the kappas.
"kappas"



#' An example of how to specify markers mapped onto chromosomes
#'
#' This is a named list.  Each component is named with the number of the chromosome.
#' In this case, those numbers are: "4", "9", "14", "19", and "21". Each  component
#' is itself a named list that holds the position in centiMorgans and the allele frequencies
#' of the alleles at the marker.
#'
#'
#' @format A named list of chromosomes upon which mapped markers are available.  The component
#' corresponding to a chromome is itself a list which may or may not be named.  If it is named
#' then its names are the names of the loci. The loci should be in sorted order along the chromosome.
#' Each locus within a chromosome is represented as a list with the following components.
#' \describe{
#'   \item{pos}{Position in centiMorgans along the chromosome. This is a little nonstandard, but you
#'   can think of the first locus as having a certain rate of recombination with the very beginning of
#'   the chromosome, and then pos is the cumulative sum of the inter-locus recombination fractions. }
#'   \item{freqs}{The allele frequencies at the locus}
#' }
"markers_on_map"



#' A long format way of specifying markers.
#'
#' This is a tbl_df data frame (a la dplyr) that provides an example of how to
#' specify markers.  This is particularly useful for pumping stuff into Mendel.
#' And it is way easier than having everything in a big named list.
#'
#'
#' @format A data frame with necessary columns:
#' \describe{
#'   \item{Chrom}{The name of the chromosome.  This should always be an integer (or maybe "X" or "Y")}
#'   \item{Locus}{The name of the locus.  This can be any string.}
#'   \item{Pos}{The position of the locus in \emph{bases} along the chromosome.}
#'   \item{Allele}{The name of allele.  Any string is OK}
#'   \item{Freq}{The frequency of the allele in the population}
#'   \item{AlleIdx}{The index of the allele.  These have to go from 1 up to the number of
#'   alleles at the locus.  Order is immaterial. }
#'   \item{LocIdx}{The index of the locus within the chromosome.  This goes from 1 up to the
#'   number of loci on the chromosome.  These should be ordered according to Pos}
#'
#' }
"long_markers"


#' 40 dinucleotide markers for a quick test data set
#'
#' This was made by running the example in \code{\link{reindix_markers}}.
#' See \code{\link{long_markers}} for more info.
"markers40"


#' 100 microhaplotype markers for a quick test data set
#'
#' This was made by running some code in scratch/create-data.R.
#' See \code{\link{long_markers}} for more info about the columns.
#' Note that Chrom is a character here (since they are unmapped) and
#' are all the same.  Pos is actually ordered the way I want them to be.
#' That turns out to be necessary for proper sorting with \code{\link{reindex_markers}}.
"microhaps"


#' Microhaplotype data from kelp rockfish \emph{Sebastes atrovirens}
#'
#' These are data from 165 amplicons sequenced in two different runs of an Illumina
#' MiSeq machine at the SWFSC lab in Santa Cruz.  Genotypes of roughly 150 individuals
#' were called an resolved into haplotypes using the fact that each short read comes
#' from a single chromosome.
#' @format  This tbl_df-ed data frame includes 825 alleles/haplotypes from 165 genomic regions and
#' is in the format of \code{\link{long_markers}}. However since true positions
#' of these markers in the genome are not known, it is instructive to see from these data
#' how to insert them into the format of \code{\link{long_markers}}.  Simply the column
#' \code{Chrom} simply has "GTseq" in it, denoting that these are markers obtained from
#' a GTseq procedure.  \code{Locus} has the name of each locus in it.  The column
#' \code{Pos} has a simple index for each locus.  They don't denote genomic positions, but they
#' are useful for sorting things if you run the data through \code{\link{reindex_markers}}.  The
#' \code{Allele} column gives the sequences of each allele/haplotype observed in the population.
"sebastes"


#' Microhaplotype data from kelp rockfish \emph{Sebastes atrovirens} with pretend linkage positions
#'
#' These are data from 165 amplicons sequenced in two different runs of an Illumina
#' MiSeq machine at the SWFSC lab in Santa Cruz.  Genotypes of roughly 150 individuals
#' were called an resolved into haplotypes using the fact that each short read comes
#' from a single chromosome. These were randomly assigned to some positions on some random
#' chromosomes for testing linkage stuff.
#' @format  This tbl_df-ed data frame includes 825 alleles/haplotypes from 165 genomic regions and
#' is in the format of \code{\link{long_markers}}. However since true positions
#' of these markers in the genome are not known I just pretended that we knew some (randomly generated)
#' locations.   The column
#' \code{Chrom} Is randomly given  \code{Locus} has the name of each locus in it.  The column
#' \code{Pos} has a random position (in base pairs) for each locus.  The
#' \code{Allele} column gives the sequences of each allele/haplotype observed in the population.
"linked_mhaps"
