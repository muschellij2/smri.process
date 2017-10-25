#' Winsorize Image
#'
#' @param x List of processed filenames
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#' @param probs passed to \code{\link{robust_window}} for Winsorization
#'
#' @return List of filenames
#' @export
#' @importFrom neurobase robust_window
winsor = function(
  x,
  verbose = TRUE,
  outdir = tempdir(),
  probs = c(0, 0.999)) {

  nii_names = names(x)
  if (length(nii_names) != length(x)) {
    stop("x must be a named vector or named list")
  }

  if (verbose > 0) {
    message("Winsorizing Image")
  }
  ####################################################
  # Dropping empty dimensions
  ####################################################
  fnames = file.path(
    outdir,
    paste0(nii_names,
           "_winsor.nii.gz"))
  names(fnames) = nii_names

  if (!all_exists(fnames)) {
    x = check_nifti(x)
    x = lapply(x, robust_window, probs = probs)

    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, x, fnames)
  }
  rm_neck = lapply(
    fnames,
    identity)
  return(rm_neck)
}
