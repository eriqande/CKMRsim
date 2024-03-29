

# The mendel package does not seem to be reliably distributed any longer.
# The web page is down, and has been for years.   Here just set up a
# system to download the binaries directly into the R library.


#' Download the mendel binary and install it where CKMRpop expects it
#'
#' This checks the operating system and installs the correct version
#' of mendel for you operating system (Mac, Linux, or Windows).  To install
#' the binary it downloads it from Eric's mbin repo on GitHub.  It also
#' installs in the same directory an intel libiomp DLL on Windows.
#' @param Dir the directory to install mendel into.  Because of restrictions
#' on functions writing to the user's home filespace, this is set, by default,
#' to a temporary directory.  But to really use this function to install mendel,
#' this parameter must be set to `system.file(package = "CKMRsim")`.
#' @export
#' @return No return value.  Called for side effect of installing the 'mendel' binary.
#' @examples
#' \dontrun{
#' install_mendel(Dir = system.file(package = "CKMRsim"))
#' }
install_mendel <- function(
    Dir = tempfile()
) {

  if(Dir != system.file(package = "CKMRsim")) {
    message("\n*** Note: To properly install mendel, the function install_mendel() must be called like this: ***\n\n    install_mendel(Dir = system.file(package = \"CKMRsim\"))
\n*** The current invocation of install_mendel() will not properly install it. ***")
  }

  # first check the OS
  Sys <- Sys.info()["sysname"]

  if(!(Sys %in% c("Darwin", "Linux", "Windows"))) {
    stop(paste("mendel binary not available for operating system ", Sys, collapse = ""))
  }

  # have a variable to hold the base GitHub address:
  Git_base <- "https://github.com/eriqande/mbin/raw/main/bin/"


  Git_full <- paste(Git_base, Sys, "/mbin.zip", sep = "")

  # get the destination path and create the directory
  Dest_dir <- file.path(Dir, "bin")
  dir.create(Dest_dir, showWarnings = FALSE, recursive = TRUE)

  # record destination file name
  Dest_file <- file.path(Dest_dir, "mbin.zip")

  # now, download the file and save to the correct destination:
  utils::download.file(url = Git_full, destfile = Dest_file)

  # unzip the contents (these are in a file called mendel)
  unzip(Dest_file, exdir = Dest_dir)

  # get the path to the result
  Bin_file <- file.path(Dest_dir, "mendel")

  # all sorts of Windows stuff here:
  #  1. We need to make the Bin_file name Mendel.exe
  #  2. We put a dll in the same location and hope that works.
  if(Sys == "Windows") {
    Bin_file <- file.path(Dest_dir, "Mendel.exe")
  }

  # finally, change the file permissions to be user and group executable and writable
  # and world readable
  Sys.chmod(Bin_file, mode = "0774", use_umask = FALSE)

  # and remove the zip
  file.remove(Dest_file)

}

