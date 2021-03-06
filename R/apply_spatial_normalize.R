#' Apply Spatial Normalization
#'
#' @param x List of filenames
#' @param template Resampling to \code{1x1x1} (template = none) or
#' Template to register to - either MNI or Eve template
#' @param template_fname Filename of Template
#' @param fwdtransforms Transforms to template if not resampling
#' @param suffix Name to append to the image filename
#' @param interpolator interpolation passed to \code{\link{ants_apply_transforms}}
#' @param outdir Output directory
#' @param verbose Print diagnostic messages
#' @param copy_origin Copy image origin from \code{target},
#' using \code{\link{antsCopyOrigin}}
#' @param origin the origin to set for the image if
#' \code{copy_origin = TRUE}
#' @return List of filenames of normalized data
#' @export
apply_spatial_normalize = function(
  x,
  template = c("none", "Eve", "MNI"),
  template_fname,
  fwdtransforms = NULL,
  suffix = "",
  interpolator = "lanczosWindowedSinc",
  outdir = tempdir(),
  verbose = TRUE,
  copy_origin = TRUE,
  origin = NULL
) {


  template = match.arg(template)
  native = template == "none"

  if (native) {
    app = "_resampled"
  } else {
    app = paste0("_regto", template)
  }


  suffix = paste0(suffix, app)

  nii_names = names(x)
  fnames = file.path(
    outdir,
    paste0(nii_names,
           suffix,
           ".nii.gz"))
  names(fnames) = nii_names

  if (!all_exists(fnames)) {
    if (native) {
      if (verbose > 0) {
        message("Resampling Image")
      }
      # resample to 1x1x1 for
      resampled = lapply(
        x,
        resample_to_target,
        target = template_fname,
        verbose = verbose > 1,
        interpolator = interpolator,
        copy_origin = copy_origin)

      mapply(function(img, fname){
        writenii(img, filename = fname)
      }, resampled, fnames)

    } else {
      if (verbose > 0) {
        message(paste0(
          "Apply Registration to ", template, " Template")
        )
      }
      set_the_origin = function(x, origin) {
        x = check_ants(x)
        t2 = ANTsRCore::antsImageClone(x)
        antsSetOrigin(t2, as.numeric(origin))
        return(t2)
      }
      if (copy_origin) {
        if (is.null(origin)) {
          stop(paste0(
            "Trying to apply spatial normalize to template",
            " and copying ",
            "the origin, but no origin given!"
          ))
        }
      }
      resampled = lapply(
        x,
        function(r) {
          if (copy_origin) {
            r = set_the_origin(r, origin)
          }
          ants_apply_transforms(
            moving = r,
            fixed = template_fname,
            transformlist = fwdtransforms,
            interpolator = interpolator)
        })
      mapply(function(img, fname){
        writenii(img, filename = fname)
      }, resampled, fnames)

    }
  }
  resampled = lapply(fnames, identity)

  return(resampled)
}