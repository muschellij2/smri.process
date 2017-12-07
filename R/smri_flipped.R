#' SMRI Flipped Image (x-y flipped)
#'
#' Flips an image and registers the image to the original.
#' Assumes RL orientation
#' @param x List of images
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#' @param mask mask to normalize over
#' @param suffix what the append to end of the output filename
#' @param interpolator interpolation done in \code{\link{registration}}
#' @param typeofTransform type of transformed used, passed to
#' \code{\link{registration}}
#' @return List of filenames
#' @export
#'
#' @importFrom neurobase flip_img
smri_flipped = function(
  x,
  mask = NULL,
  outdir = tempdir(),
  verbose = TRUE,
  suffix = "_flip",
  typeofTransform = "SyN",
  interpolator = "Linear"
){


  nii_names = names(x)
  if (length(nii_names) != length(x)) {
    stop("x must be a named vector or named list")
  }

  if (verbose > 0) {
    message("Flipping and Registering Image")
  }
  #################################
  # Flip and Register
  #################################
  fnames = file.path(
    outdir,
    paste0(nii_names,
           suffix,
           ".nii.gz"))
  names(fnames) = nii_names

  if (!all_exists(fnames)) {

    x = check_nifti(x)
    if (!is.null(mask)) {
      mask = check_nifti(mask)
    }
    # data is already RPI
    flip_x = lapply(x, function(img) {
      min_img = min(img)
      ##########################
      # needs to be positive!!!
      ##########################
      if (min_img < 0) {
        img = img - min_img + 1
      }
      if (!is.null(mask)) {
        img = mask_img(img, mask = mask)
      }

      fimg = flip_img(
        img = img,
        x = TRUE,
        y = FALSE,
        z = FALSE)
      reg = registration(
        filename = fimg,
        template.file = img,
        typeofTransform = typeofTransform,
        interpolator = interpolator,
        verbose = verbose > 1)
      out = reg$outfile
      if (min_img < 0) {
        out = out + min_img - 1
      }
      if (!is.null(mask)) {
        out = mask_img(out, mask = mask)
      }
      return(out)
    })

    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, flip_x, fnames)
  }


  flip_x = lapply(
    fnames,
    identity)

  return(flip_x)
}

