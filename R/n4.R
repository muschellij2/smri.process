
#' N4 Bias Field Correction
#'
#' @param x List of processed filenames
#' @param outdir Output directory
#' @param mask Mask for Bias-field correction
#' @param verbose print diagnostic messages
#' @param suffix Name to append to the image filename
#'
#' @return List of filenames
#' @export
#' @importFrom plyr llply
#' @importFrom extrantsr bias_correct
n4_raw = function(
  x,
  verbose = TRUE,
  mask = NULL,
  outdir = tempdir(),
  suffix = "_N4") {

  nii_names = names(x)
  if (length(nii_names) != length(x)) {
    stop("x must be a named vector or named list")
  }

  if (verbose > 0) {
    message("N4 Bias Field Correction")
  }
  #################################
  # Bias correct the images using N4
  #################################
  fnames = file.path(
    outdir,
    paste0(nii_names,
           suffix,
           ".nii.gz"))
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
        verbose = verbose > 1,
        mask = mask,
        retimg = FALSE)
    }, x, fnames, SIMPLIFY = FALSE)
  }

  n4 = llply(fnames,identity)

  return(n4)
}