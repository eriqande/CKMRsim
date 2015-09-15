

# functions to interface with the GENEDROP program in the package
# MORGAN, to simulate linked marker data.

#' convert a pedigree into control files used by morgan to drop genes
#'
#' @param ped A data frame with five columns: Kid, Pa, Ma, Sex, Observed.
#' The entries in the first three columns are strings or numbers, and in the
#' remaining two are integers: Sex: 1 - male, 2 - female, 0 - unknown.
#' Observed: 1 - observed, 2 - not observed. Note the founders in the pedigree
#' have to appear as kids in a record that has 0 and 0 for the Pa and Ma.
#' @param markers A list of chromosomes with markers upon them in the format of
#' \code{\link{markers_on_map}}
#' @param Dir  the directory where the output files should be placed.
#' @param seed a vector of 2 integers, each of them positive but less than 5e08  The
#' first one has to be odd.
#' @examples
#'    data(pedigrees)
#'    ped_df2morgan_par_and_ped(pedigrees$FS)
ped2morgan_input <- function(ped, markers, seed, Dir = ".") {
  parf <- file.path(Dir, "morgan_parfile")
  pedf <- file.path(Dir, "morgan_pedfile")
  markerseedf <- file.path(Dir, "marker.seed")

  # check the seeds:
  if(any(seeds < 0)) stop("Both seeds have to be greater than 0")
  if(any(seeds > 5e08)) stop("Both seeds have to be less than or equal to 5e08")
  if(seeds[1] %% 2 != 1) stop("The first seed has to be odd")

  ### fill the seedfile
  seedy <- paste("0x", as.hexmode(seeds), sep = "", collapse = " ")
  cat(paste("set marker seeds", seedy), file = markerseedf, sep = "\n")


  ### make the pedfile
  cat("************\n", file = pedf)
  write.table(ped, row.names = FALSE, col.names = FALSE, quote = FALSE, file = pedf, append = TRUE)

  ### make the parfile
  # first the preambly-stuff
  cat("input pedigree record names 3 integers 2\nset printlevel 5\noutput pedigree chronological", file = parf, sep = "\n")
  cat(paste("input seed file \"", markerseedf, "\"", sep = ""), file = parf, append = TRUE, sep = "\n")
  cat(paste("output overwrite seed file \"", markerseedf, "\"", sep = ""), file = parf, append = TRUE, sep = "\n")
  cat("output marker seeds only", file = parf, append = TRUE, sep = "\n")
  cat(paste("input pedigree file \"", pedf, "\"", sep = ""), file = parf, append = TRUE, sep = "\n")
  cat(paste("input pedigree size", nrow(ped)), file = parf, append = TRUE, sep = "\n")
  cat("input pedigree record observed present", file = parf, append = TRUE, sep = "\n")

  # next the information about the markers to be simulated
  # simulation commands for each chromosome
  cat(paste("simulate chrom", names(markers), "markers"), file = parf, append = TRUE, sep = "\n")

  # map positions for each chromosome
  dump <- lapply(names(markers), function(x) {
    posstring <- paste(diff(sapply(markers[[x]], function(y) y$pos)), collapse = " ")
    cat(paste("map chrom", x, "marker dist", posstring), file = parf, append = TRUE, sep = "\n")
  })

  # allele freqs for each marker on each chromosome
  dump <- lapply(names(markers), function(x) {
    xxx <<- 0
    lapply(markers[[x]], function(y) {
      xxx <<- xxx + 1
      freqstring <- paste(y$freqs, collapse = " ")
      cat(paste("set chrom", x, "markers", xxx, "allele freqs", freqstring), file = parf, append = TRUE, sep = "\n")
    })
  })


}
