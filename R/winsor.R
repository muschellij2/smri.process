#' Winsorize Image
#'
#' @param x List of processed filenames
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#' @param non_zero Should zeroes be excluded from the calculation
#' of quantiles? Passed to \code{\link{robust_window}}
#' @param probs passed to \code{\link{robust_window}} for Winsorization
#' @param suffix Name to append to the image filename
#'
#' @return List of filenames
#' @export
#' @importFrom neurobase robust_window
winsor = function(
  x,
  verbose = TRUE,
  outdir = tempdir(),
  non_zero = TRUE,
  probs = c(0.001, 0.999),
  suffix = "_winsor") {

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
           suffix,
           ".nii.gz"))
  names(fnames) = nii_names

  if (!all_exists(fnames)) {
    x = check_nifti(x)
    x = lapply(x, robust_window, probs = probs,
               non_zero = non_zero)

    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, x, fnames)
  }
  rm_neck = lapply(
    fnames,
    identity)
  return(rm_neck)
}
