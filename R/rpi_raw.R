#' @title RPI Orientation of Images
#'
#' @param x List of images or character vector.  Must be named with
#' the imaging modalities. T1 and FLAIR must be included
#' @param gold_standard Gold Standard image/filenameLesion
#' @param verbose print diagnostic messages
#'
#' @return List of reoriented output filenames
#' @export
#' @importFrom fslr rpi_orient_file
#' @importFrom neurobase checkimg
rpi_raw = function(
  x,
  gold_standard = NULL,
  gs_space = NULL,
  verbose = TRUE) {

  nii_names = names(x)
  nii_names = toupper(nii_names)
  if (length(nii_names) != length(x)) {
    stop("x must be a named vector or named list")
  }
  if (!all(c("T1", "FLAIR") %in% nii_names)) {
    stop("T1 and FLAIR must be in the images vector")
  }

  x = checkimg(x)
  x = lapply(x, identity)
  names(x) = nii_names

  if (!is.null(gold_standard)) {
    if (is.null(gs_space)) {
      stop("gs_space must be specified if gold_standard specified")
    }
    gs_space = match.arg(gs_space, choices = names(x))
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
    rpi_orient_file,
    verbose = verbose)
  return(x)
}
