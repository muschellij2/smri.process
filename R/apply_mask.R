#' Apply Brain mask
#'
#' @param x List of images or character vector.  Must be named with
#' the imaging modalities.
#' @param mask Brain Mask image/filename
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#' @param suffix Name to append to the image filename
#'
#' @return List of masked brains
#' @export
#' @importFrom neurobase mask_img check_nifti
apply_mask = function(
  x,
  mask,
  outdir = tempdir(),
  verbose = TRUE,
  suffix = "_brain") {

  nii_names = names(x)
  nii_names = toupper(nii_names)
  if (length(nii_names) != length(x)) {
    stop("x must be a named vector or named list")
  }

  mask = check_nifti(mask)


  if (verbose > 0) {
    message("Applying Mask")
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
    rm_neck = llply(
      x,
      mask_img,
      mask = mask,
      .progress = ifelse(verbose, "text", "none"))

    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, rm_neck, fnames)
  }
  rm_neck = lapply(
    fnames,
    identity)
  return(rm_neck)
}
