#' RPI Orientation of Images
#'
#' @param x List of images or character vector.  Must be named with
#' the imaging modalities.
#' @param gold_standard Gold Standard image/filename, if applicable
#' @param verbose print diagnostic messages
#'
#' @return List of reoriented output filenames
#' @export
#' @importFrom fslr rpi_orient_file is_rpi enforce_form
#' @importFrom neurobase checkimg
rpi_raw = function(
  x,
  gold_standard = NULL,
  verbose = TRUE) {

  nii_names = names(x)
  nii_names = toupper(nii_names)
  if (length(nii_names) != length(x)) {
    stop("x must be a named vector or named list")
  }
  # if (!all(c("T1", "FLAIR") %in% nii_names)) {
  #   stop("T1 and FLAIR must be in the images vector")
  # }

  x = checkimg(x)
  x = lapply(x, identity)
  names(x) = nii_names

  if (!is.null(gold_standard)) {
    gold_standard = checkimg(gold_standard)
  }

  x$GOLD_STANDARD = gold_standard

  # REMOVE NULL
  nulls = sapply(x, is.null)
  x = x[!nulls]

  nii_names = names(x)

  if (verbose) {
    message("Image Orientation to RPI")
  }
  x = lapply(
    x,
    function(fname) {
      fname = fslr::enforce_form(file = fname)
      if (!fslr::is_rpi(fname)) {
        res = rpi_orient_file(fname, verbose = verbose > 1)
      } else {
        res = list(img = fname)
      }
      return(res)
    })
  x = lapply(x, function(r) {
    r$img
  })
  return(x)
}
