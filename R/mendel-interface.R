

#' convert a long-format data frame of loci into a vector of lines for a Mendel def file
#'
#' This makes a vector of text that has carriage returns and things in it as necessary.
#' So, to include it in a def file, if the return from this function is called
#' D you can just do cat(D, sep = "", file = myDef.txt), for example.
#'
#' Because it appears that Mendel requires unique locus names we will name them
#' Chrom-underscore-Loc.
#'
#' @param df  A data frame in the format of \code{\link{long_markers}}.
#' @examples
#' data(long_markers)
#' D <- markers2mendel_def_lines(long_markers)
#' cat(D[1:100], sep = "")
markers2mendel_def_lines <- function(df) {
  tmp <- df %>%
    ungroup() %>%
    arrange(Chrom, LocIdx, AlleIdx) %>%
    group_by(Chrom, LocIdx) %>%
    mutate(LocLine = paste(Chrom, "_", LocIdx, ",", Chrom, ",", n_distinct(AlleIdx), sep = ""),
           LocScrub = ifelse(AlleIdx == 1, paste(LocLine,"\n"), ""),
           AlleLine = paste(AlleIdx, ",", Freq, "\n", sep = ""),
           TextVec = paste(LocScrub, AlleLine, sep = "   ")
    )

  tmp$TextVec
}


#' create a vector of Mendel formatted map file strings
#'
#' This makes a vector of text that has carriage returns and things in it as necessary.
#' So, to include it in a map file, if the return from this function is called
#' D you can just do cat(D, sep = "", file = myMap.txt), for example.
#'
#' @param df  A data frame in the format of \code{\link{long_markers}}.
#' @examples
#' data(long_markers)
#' D <- markers2mendel_map_lines(long_markers)
#' cat(D[1:100], sep = "")
markers2mendel_map_lines <- function(df) {
  tmp <- df %>%
    arrange(Chrom, LocIdx, AlleIdx) %>%
    group_by(Chrom) %>%
    mutate(ender = ifelse(LocIdx == max(LocIdx), "\n\n", "\n")) %>%
    ungroup %>%
    group_by(Chrom, LocIdx) %>%
    summarize(map_lines = paste(Chrom, "_", LocIdx, ", ", Pos, ender , sep = "")[1])

  tmp$map_lines
}


#' for making a pedigree file with no observed genotypes on it
pedigree2mendel_ped_file <- function(df, name, filename = NA) {


  df2 <- df %>%
    tbl_df %>%
    mutate(pedname = name,
           MZtwin_stat = "",
           Panew = ifelse(Pa == 0, "", Pa),
           Manew = ifelse(Ma == 0, "", Ma),
           Sexnew = ifelse(Sex == 0, "", c("M", "F")[Sex]),
           TwinStatus = ""
           ) %>%
    select(pedname, Kid, Panew, Manew, Sexnew, TwinStatus)

    if(is.na(filename)) {
      return(df2)
    }

    write.table(df2, file = filename, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")

    invisible(NULL)
}


#' prepend an ID to a mendel input or output file
#'
#' This is just a little helper fucntion to ensure consistency in file naming.
#' @param ID The ID to be prepended.
#' @param suff The final part of the file name (like Def.in, etc.)
id_prepend <- function(ID, suff) {
  paste(ID, "-", suff, sep = "")
}


#' return a list of values for a Mendel definitions file for gene dropping
#'
#' @param ID the name/prefix to be given to all the files involved in this run.
mendel_control_list <- function(ID, Reps, Seed) {
  ret <- list(
    # input files
    DEFINITION_FILE = id_prepend(ID, "Def.in"),
    MAP_FILE = id_prepend(ID, "Map.in"),
    PEDIGREE_FILE = id_prepend(ID, "Ped.in"),

    # output files
    NEW_PEDIGREE_FILE = id_prepend(ID, "Ped.out"),
    OUTPUT_FILE = id_prepend(ID, "Mendel.out"),
    SUMMARY_FILE = id_prepend(ID, "Summary.out"),

    # analysis options
    ANALYSIS_OPTION = "Gene_dropping",
    REPETITIONS = Reps,
    SEED = Seed,
    MODEL = 2,
    KEEP_FOUNDER_GENOTYPES = "False",
    MISSING_DATA_PATTERN = "None",
    GENE_DROP_OUTPUT = "Unordered",
    MAP_DISTANCE_UNITS = "bp"
  )

  ret
}




#' Write all necessary files for Mendel to do a gene-dropping run
#'
#' Writes all files in the directory Dir but does not launch mendel.
#' @param ID  The ID or identifier of the run.  This string will be prepended to all input and output files.
#' @param Reps The number of replicate data sets to simulate.
#' @param Seed The random seed to use.  This must be a single integer between 1 and 30,000.
#' @param Markers A data frame in the format of \code{\link{long_markers}} that gives information
#' about marker positions and allele frequencies.
#' @param Pedigree A data frame formatted like any component of \code{\link{pedigrees}} which holds the
#' pedigree to drop genes upon.
#' @param Dir The directory in which to write the files.  It must be created already.  Defaults to the
#' current working directory.
#' @examples
#' # make a temp directory to place things into:
#' tmpDir <- tempdir()
#' write_all_mendel_files("mendel-example", 10, 1234, long_markers, pedigrees$HS, Dir = tmpDir)
write_all_mendel_files <- function(ID, Reps, Seed, Markers, Pedigree, Dir = ".") {
  # get the lines
  DL <- markers2mendel_def_lines(Markers)
  ML <- markers2mendel_map_lines(Markers)
  MCL <- mendel_control_list(ID, Reps, Seed)

  # write the pedigree file
  pedigree2mendel_ped_file(df = Pedigree, name = "xyz", filename = file.path(Dir, id_prepend(ID, "Ped.in")))

  # write the control file
  cat(paste(names(MCL), MCL, sep = " = "), sep = "\n", file = file.path(Dir, id_prepend(ID, "Control.in")))

  # write the definitions file
  cat(DL, sep = "", file = file.path(Dir, id_prepend(ID, "Def.in")))

  # write the map file
  cat(ML, sep = "", file = file.path(Dir, id_prepend(ID, "Map.in")))
}


#' slurp up the pedigree output file from a Mendel genedrop and get the needed genotypes
#'
#' This assumes that the target individuals are numbered 1 and 2, and it pitches the data
#' for everyone else, then it converts the genotypes into a big long format data frame.
#' @param Pedfile  The name of the pedigree output file
mendel_outped2genos <- function(Pedfile) {
  x <- readLines(Pedfile)

  x2 <- stringr::str_replace(x, "^.*xyz *, *", "") # make another vector with the crap pulled off the front

  x3 <- x[stringr::str_detect(x2, "^[12] ")] # keep only the lines referring to the target individuals

  # I stopped here because readLines is so slow!
}

#' read genotypes of target individuals out of Mendel output pedigree file.
#'
#' this is a version that uses awk and sed that I am sure will be much faster than
#' the standard R readLines approach.  But only works if those tools are installed
#' on the system.
#' @param OutPed  path to the mendel output pedigree file
#' @param Def path to the mendel defs file.  This is used to get the number of alleles at
#' each locus, so, the model loci in the mendel run must be equal to the loci listed in Def.
#' @return Returns a tbl_df with four columns: 1) Rep, a base-1 index of which simulated replicate
#' it is; 2) Indiv, the index of the target individual (typically 1 or 2); 3) Locus, the index
#' of the locus; 4) Geno, the genotype of the individual as a single number.
#' @examples
#' # first prepare files for mendel simulation and run it in tmpDir
#' # with naming prefix "mendel-example"
#' example(run_mendel)
#'
#' # then grab the results
#' results <- fast_mendel_outped2genos(file.path(tmpDir, "mendel-example-Ped.out"), file.path(tmpDir, "mendel-example-Def.in"))
#'
#' # finally show the results:
#' results
fast_mendel_outped2genos <- function(OutPed, Def) {

  # I couldn't get pipe() to work, so I will just write to a tempfile and read back in
  tf <- tempfile()
  scrpath <- system.file("script", "slurp_outped.sh", package = "CKMRsim")
  command <- paste(scrpath, OutPed, Def, ">", tf)
  system(command)
  ret <- read.table(tf, stringsAsFactors = FALSE, colClasses = "integer") %>%
    tbl_df %>%
    setNames(c("Rep", "Indiv", "Locus", "Geno")) %>%
    mutate(Rep = Rep + 1L)
  ret
}



#' Given the operating system, check for mendel program in the typical location
#'
mendelBin <- function() {
  if(Sys.info()['sysname'] == "Darwin") {  # we are on a Mac
    path <- "/Applications/Mendel/mendel"
  } else {
    stop("Sorry, we don't know where to expect the mendel binary on an operating system of type ", Sys.info()['sysname'])
  }
  if(!file.exists(path)) stop("Didn't find mendel binary where we expected it at ", path)
  path
}

#' Run the mendel binary in directory Dir using control file Control
#'
#' @param Dir the directory in which to run mendel
#' @param Control the "control file" with which to run mendel.  If not an absolute path
#' it should be given relative to the Dir directory.
#' @examples
#' # first prepare the input files and define tmpDir
#' example(write_all_mendel_files)
#'
#' # then run mendel
#' run_mendel(tmpDir, "mendel-example-Control.in")
run_mendel <- function(Dir, Control) {
  COMM <- paste("cd",
                Dir,
                ";",
                mendelBin(),
                "-c",
                Control
                )
  system(COMM)
}
