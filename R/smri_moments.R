#' Structural MRI Moment images
#'
#' @param radius array of values for neighborhood radius (in voxels),
#' passed to \code{\link{neighborhood}}
#' @param moments Moments to create, mean, standard deviation,
#' skew, kurtosis, gradient, z-score (mean/sd)
#' @param x List of images
#' @param suffix Name to append to the image filename.  The moments will
#' be appended to this suffix.  Should likely end with an underscore (_) or
#' other separator
#' @param outdir output directory
#' @param verbose print diagnostic messages
#' @param mask brain mask for moment
#' @param ... Additional arguments to \code{\link{create_moment}}
#'
#' @importFrom extrantsr create_moment
#'
#' @return A list of filenames, separated by imaging modality, containing
#' filenames for all the moments specified
#' @export
smri_moments = function(
  x,
  mask = NULL,
  radius = rep(1, 3),
  suffix = "_",
  outdir = tempdir(),
  verbose = TRUE,
  moments = c("mn", "sd", "skew",
            "kurt", "grad", "z"),
  ...
) {


  nii_names = names(x)
  if (length(nii_names) != length(x)) {
    stop("x must be a named vector or named list")
  }



  if (verbose > 0) {
    message("Creating Moment Images")
  }


  all_fnames = lapply(nii_names, function(img_name) {
    x = file.path(
      outdir,
      paste0(img_name,
             suffix,
             moments,
             ".nii.gz"))
    names(x) = moments
    x
  })
  names(all_fnames) = nii_names

  vec_fnames = unlist(all_fnames)

  if (!all_exists(vec_fnames)) {

    for (iimg in nii_names) {

      fnames = all_fnames[[iimg]]

      if (!all_exists(fnames)) {
        if (verbose > 0) {
          message(paste0("Moments for ", iimg))
        }
        img = x[[iimg]]
        mom_imgs = create_moment(
          img,
          radius = radius,
          mask = mask,
          retimg = TRUE,
          verbose = verbose > 1
          )
        mom_imgs = mom_imgs[moments]
        fnames = fnames[moments]
        if (verbose > 0) {
          message("Writing images")
        }
        mapply(function(img, fname){
          writenii(img, filename = fname)
        }, mom_imgs, fnames)
      }
    }
  }

  return(all_fnames)

}