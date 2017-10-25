#' Remove Neck (Doubly)
#'
#' @param x List of processed filenames
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#' @param suffix Name to append to the image filename
#'
#' @return List of filenames
#' @export
#' @importFrom extrantsr double_remove_neck
#' @importFrom fslr mni_fname
#' @importFrom plyr llply
noneck = function(
  x,
  verbose = TRUE,
  outdir = tempdir(),
  suffix = "_noneck") {

  nii_names = names(x)
  if (length(nii_names) != length(x)) {
    stop("x must be a named vector or named list")
  }

  if (verbose > 0) {
    message("Removing necks")
  }
  ####################################################
  # Making fname stubs
  ####################################################
  fnames = file.path(
    outdir,
    paste0(nii_names,
           suffix,
           ".nii.gz"))
  names(fnames) = nii_names

  if (!all_exists(fnames)) {
    mni.template.file = fslr::mni_fname(
      mm = "1",
      brain = TRUE,
      mask = FALSE)
    mni.template.mask = fslr::mni_fname(
      mm = "1",
      brain = TRUE,
      mask = TRUE)

    noneck = llply(
      x,
      function(r) {
        r = check_nifti(r) # new version will make sure no extensions
        double_remove_neck(
          r,
          template.file = mni.template.file,
          template.mask = mni.template.mask,
          verbose = verbose > 1
        )
      },
      .progress = ifelse(verbose, "text", "none"))
    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, noneck, fnames)
  }
  noneck = lapply(
    fnames,
    identity)

  return(noneck)
}