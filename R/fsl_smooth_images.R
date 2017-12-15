#' Smoothing sMRI data with FSL
#'
#' @param x List of images
#' @param sigma vector of sigmas for smoothers in FSL
#' @param suffix Name to append to the image filename.  The sigma will
#' be appended to this suffix.
#' @param outdir output directory
#' @param verbose print diagnostic messages
#' @param mask brain mask for smoother
#' @param ... Additional arguments to \code{\link{fsl_smooth}}
#'
#' @return List of the images, brain mask, suffix, and output directory
#' @export
#'
#' @importFrom fslr fsl_smooth
fsl_smooth_images = function(
  x,
  sigma = c(3, 5, 10, 20),
  mask = NULL,
  outdir = tempdir(),
  suffix = "_smooth",
  verbose = TRUE,
  ...
) {

  nii_names = names(x)
  nii_names = toupper(nii_names)
  if (length(nii_names) != length(x)) {
    stop("x must be a named vector or named list")
  }
  x = checkimg(x)
  x = lapply(x, identity)
  names(x) = nii_names


  nii_names = nii_names[ !(nii_names %in% "GOLD_STANDARD") ]
  x = x[nii_names]


  all_fnames = lapply(sigma, function(isigma) {
    x = file.path(
      outdir,
      paste0(nii_names,
             suffix,
             isigma,
             ".nii.gz"))
    names(x) = nii_names
    x
  })
  names(all_fnames) = sigma

  vec_fnames = unlist(all_fnames)

  if (!all_exists(vec_fnames)) {
    if (!is.null(mask)) {
      mask = checkimg(mask)
    }

    if (length(sigma) > 0) {

      for (isigma in sigma) {
        if (verbose) {
          message(paste0("Smoothing Images: Sigma = ", isigma))
        }

        ####################################################
        # Making fname stubs
        ####################################################
        fnames = all_fnames[[as.character(isigma)]]
        fnames = unlist(fnames)
        if (!all_exists(fnames)) {

          if (!is.null(mask)) {
            # New 2017May03
            smoothed_mask = fslr::fsl_smooth(
              file = mask,
              sigma = isigma,
              smooth_mask = FALSE,
              verbose = verbose > 1,
              ...)
          } else {
            smoothed_mask = NULL
          }

          ## smooth the images using fslsmooth from the fslr package
          smooth = mapply(function(infile, outfile) {
            fslr::fsl_smooth(
              file = infile,
              sigma = isigma,
              mask = mask,
              retimg = FALSE,
              smooth_mask = TRUE,
              smoothed_mask = smoothed_mask,
              outfile = outfile,
              verbose = verbose > 1,
              ...)
          }, x, fnames, SIMPLIFY = FALSE)
        }
      }
    }
  }


  return(all_fnames)

}