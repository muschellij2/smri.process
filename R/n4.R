
#' N4 Bias Field Correction
#'
#' @param x List of processed filenames
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#'
#' @return List of filenames
#' @export
#' @importFrom plyr llply
#' @importFrom extrantsr bias_correct
n4_raw = function(
  x,
  verbose = TRUE,
  outdir = tempdir()) {

  nii_names = names(x)
  if (verbose > 0) {
    message("N4 Bias Field Correction")
  }
  #################################
  # Bias correct the images using N4
  #################################
  fnames = file.path(
    outdir,
    paste0(nii_names,
           "_N4.nii.gz"))
  names(fnames) = nii_names


  if (!all_exists(fnames)) {
    #################################
    # Swap image for remove_neck.
    # Can use swapdim in remove_neck, but
    # registration seems to work better if
    # images are in std RPI
    #################################
    n4 = mapply(function(infile, outfile) {
      bias_correct(
        file = infile,
        correction = "N4",
        outfile = outfile,
        retimg = FALSE)
    }, x, fnames, SIMPLIFY = FALSE)
  }

  n4 = llply(fnames,identity)

  return(n4)
}