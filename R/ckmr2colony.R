



#' make a rudimentary colony input file
#'
#' This takes a long format data frame of CKMR frequencies and a corresponding
#' one of genotypes of individuals (gotta have LocIdx and alle1 and alle2 columns,
#' which give the integer representation of the alleles).  Then it spits out
#' a colony file with allele freqs either specified (by default) or not.
#'
#' The idea is that you use this to follow up on clusters of pairs that look like
#' they might be sibs of some sort.  This spits out a Colony2.dat file that can
#' be run with Colony.  To change different settings and options for the run, you
#' just hand edit that output file.
#'
#' @param CKF a data frame like long_markers that has the allele freqs and indexes, etc.
#' @param CKG a data frame of genotype data.  It needs to have these columns:
#' \describe{
#' \item{id}{the identifier of the individual}
#' \item{Locus}{the names of the loci}
#' \item{LocIdx}{the index of the locus (should correspond exactly to CKF)}
#' \item{alle1}{The integer index of the allelic type of the first gene copy in the indiv.}
#' \item{alle2}{The integer index of the allelic type of the second gene copy in the indiv.}
#' }
#' @param freqs_known set to FALSE if you want Colony to estimate the allele freqs.  Otherwise
#' it will just use the freqs from CKF (This makes a lot of sense if you have 10,000 individuals
#' but you don't want to run all of them in Colony)
#' @param outfile  path to the output file to send this to
#' @param err_rate the per-locus rate of mis-genotyped loci assumed for the analysis.  Recycles to the
#' proper length.
#' @param drop_rate the per-locus rate of assumed for the analysis. Recycles to the proper length.
#' @export
ckmr2colony <- function(CKF, CKG, freqs_known = TRUE, outfile = "Colony2.dat", err_rate = c(0.01), drop_rate = c(0.01)) {
  NumLoci <- max(CKF$LocIdx)
  cd <- colony_template()
  cd[3, 1] <- length(unique(CKG$id))  # number of individuals
  cd[4, 1] <- NumLoci

  if(freqs_known == FALSE) {
    cd[14, 1] <- 0
  }

  freqs <- afreq_str(CKF)

  # now make a nice wide frame of the genotypes, denote NAs with O's if they are present
  # and if they are absent as well
  wide_genos <- CKG %>%
    dplyr::mutate(a1 = ifelse(is.na(alle1), 0, alle1),
           a2 = ifelse(is.na(alle2), 0, alle2),
           genostr = paste(a1, a2, sep = " ")) %>%
    dplyr::select(id, LocIdx, genostr) %>%
    tidyr::spread(data = ., key = LocIdx, value = genostr, fill = "0 0")

  # now start writing stuff out to the file
  write.table(cd[1:14,], row.names = F, col.names = F, quote = F, sep = "\t", file = outfile)

  # if allele frequencies are known, spew those out
  if(freqs_known == TRUE) {
    cat(freqs, file = outfile, append = TRUE, sep = "\n")
  }

  # then spit the rest of the template out
  write.table(cd[15:21,], row.names = F, col.names = F, quote = F, sep = "\t", file = outfile, append = TRUE)

  # now spew the locus names out
  tmp <- CKG %>%
    dplyr:: group_by(Locus, LocIdx) %>%
    dplyr::tally() %>%
    dplyr::arrange(LocIdx)

  cat(tmp$Locus, file = outfile, append = TRUE, sep = " ", eol = "\n")

  # and spew out the codominant 0's
  cat(rep(0, length.out = NumLoci), sep = " ", eol = "\n", file = outfile, append = TRUE)

  # and then the genotyping error rates
  cat(rep(err_rate, length.out = NumLoci), sep = " ", eol = "\n", file = outfile, append = TRUE)
  cat(rep(drop_rate, length.out = NumLoci), sep = " ", eol = "\n", file = outfile, append = TRUE)

  # now spew out the genotypes
  write.table(wide_genos, file = outfile, append = TRUE, row.names = F, col.names = F, quote = F, sep = "\t")

  # then, finally spew the postamble out
  cat("\n\n0.0  0.0\n0  0\n\n0\n\n0\n\n0\n\n0\n\n0\n\n0\n\n0\n\n0", file = outfile, append = TRUE)
}





#### These are some non-exported helper functions ####


colony_template <- function() {
  comments <- c(
    "!   data set name",
    "!   outfile prefix",
    "!   Number of offspring in the sample",
    "!   Number of loci",
    "!   Seed for random number generator",
    "!   0/1=Not updating/updating allele frequency",
    "!   2/1=Dioecious/Monoecious species",
    "!   0/1=No inbreeding/inbreeding",
    "!   0/1=Diploid species/HaploDiploid species",
    "!   0/1=Polygamy/Monogamy for males & females",
    "!   0/1=Clone inference =No/Yes",
    "!   0/1=Scale full sibship=No/Yes",
    "!   0/1/2/3=No/Weak/Medium/Strong sibship prior; 4=optimal sibship prior",
    "!   0/1=Unknown/Known population allele frequency",
    "!   Number of runs",
    "!   1/2/3/4=short/medium/long/very long run",
    "!   0/1=Monitor method by Iterate#/Time in second",
    "!   Monitor interval in Iterate# / in seconds",
    "!   Windows version",
    "!   0/1/2=PairLikelihood score/Fulllikelihood/FPLS",
    "!   0/1/2/3=Low/Medium/High/Very high precision with Fulllikelihood")

  defaults <- c(
    "ckmr",
    "'output'",
    "NA",
    "NA",
    "4999",
    "1",
    "2",
    "0",
    "0",
    "0 0",
    "0",
    "1",
    "0",
    "1",
    "1",
    "1",
    "1",
    "10000",
    "0",
    "1",
    "1")

  data.frame(values = defaults, comments = comments, stringsAsFactors = FALSE)
}


#' prints out the integer allele names and their freqs into alternating
#' strings in a vector, so all you have to do is cat it with sep = "backslash n" to
#' get it to print to the colony file.
afreq_str <- function(CKF) {
  cl <- split(CKF, CKF$LocIdx)

  nums <- unname(unlist(lapply(cl, nrow)))
  numstr <- paste(nums, collapse = " ")

  freqs <- unlist(lapply(cl, function(x) {
    idxstr <- paste(sprintf("%d", x$AlleIdx), collapse = " ")
    freqstr <- paste(sprintf("%.8f", x$Freq), collapse = " ")
    c(idxstr, freqstr)
  }))

  c(numstr, freqs)
}
