
#' Reduce Image, Dropping Image Dimensions
#'
#' @param x List of processed filenames
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#'
#' @return List of filenames
#' @export
#' @importFrom extrantsr oMask
reduce_img = function(
  x,
  verbose = TRUE,
  outdir = tempdir()) {

  nii_names = names(x)
  if (length(nii_names) != length(x)) {
    stop("x must be a named vector or named list")
  }

  if (verbose > 0) {
    message("Dropping Empty Dimensions")
  }
  ####################################################
  # Dropping empty dimensions
  ####################################################
  fnames = file.path(
    outdir,
    paste0(nii_names,
           "_reduced.nii.gz"))
  names(fnames) = nii_names

  if (!all_exists(fnames)) {
    rm_neck = llply(
      x,
      function(nn){
        dd = mask_reduce(nn)
        return(dd$outimg)
      }, .progress = ifelse(verbose, "text", "none"))

    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, rm_neck, fnames)
  }
  rm_neck = lapply(
    fnames,
    identity)
  return(rm_neck)
}

#' @rdname reduce_img
#' @export
#' @importFrom extrantsr oMask
#' @importFrom neurobase dropEmptyImageDimensions
mask_reduce = function(x) {
  # cleanup 0 because more liberal
  mask = oMask(x, cleanup = 0)
  dd = dropEmptyImageDimensions(
    mask,
    keep_ind = TRUE,
    other.imgs = x)
  dd$outimg = dd$other.imgs
  dd$other.imgs = NULL

  return(dd)
}