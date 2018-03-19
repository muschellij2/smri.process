#' SMRI Z-score compared to a Template
#'
#' @param x List of images
#' @param mask mask to place all non-mask values to zero
#' @param template Template of the space of the mean/SD images,
#' only currently accepts Eve.
#' @param mean_imgs A vector/list of images for measures of
#' centrality - usually the mean.
#' See \code{\link{download_template_img_data}}
#' @param sd_imgs A vector/list of images for measures of
#' spread - usually the standard deviation.
#' See \code{\link{download_template_img_data}}#'
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#' @param remask should \code{\link{mask_img}}
#' be applied to the output images
#' after normalization
#' @param interpolator Interpolation done, passed to
#' \code{\link{registration}}
#' @param suffix what the append to end of the output filename
#'
#' @return List of images
#' @export
#'
#' @importFrom neurobase zscore_img quantile_img
#' @importFrom WhiteStripe whitestripe_norm whitestripe whitestripe_hybrid
template_z_score = function(
  x,
  mask = NULL,
  template = "Eve",
  mean_imgs = NULL,
  sd_imgs = NULL,
  outdir = tempdir(),
  verbose = TRUE,
  remask = TRUE,
  interpolator = "lanczosWindowedSinc",
  suffix = "_ztemp") {

  template_fname = switch(
    template,
    Eve = EveTemplate::getEvePath(what = "Brain")
    # ,
    # MNI = MNITemplate::getMNIPath(what = "Brain", res = "1mm")
  )

  nii_names = names(x)
  if (length(nii_names) != length(x)) {
    stop("x must be a named vector or named list")
  }
  if (length(mean_imgs) != length(x)) {
    stop("Mean Images need to be specified and same length as x!")
  }
  if (length(sd_imgs) != length(x)) {
    stop("SD Images need to be specified and same length as x!")
  }

  t1 = x[["T1"]]
  if (is.null(t1)) {
    stop("T1 is not in x!")
  }

  if (verbose > 0) {
    message("Z-score Template Normalization")
  }
  #################################
  # Bias correct the images using N4
  #################################
  fnames = file.path(
    outdir,
    paste0(nii_names,
           suffix,
           ".nii.gz"))
  names(fnames) = nii_names

  if (verbose > 0) {
    message(paste0(
      "Registering T1 to ", template, " Template")
    )
  }

  if (!all_exists(fnames)) {

    t1_reg = registration(
      filename = t1,
      template.file = template_fname,
      skull_strip = FALSE,
      correct = FALSE,
      verbose = verbose > 1,
      typeofTransform = "SyN",
      remove.warp = FALSE,
      interpolator = interpolator)
    keeper = function(x) {
      x[ grep("Generic|Warp", x)]
    }
    t1_reg$fwdtransforms = keeper(t1_reg$fwdtransforms)
    t1_reg$invtransforms = keeper(t1_reg$invtransforms)

    mean_imgs = check_nifti(mean_imgs)
    sd_imgs = check_nifti(sd_imgs)

    norm = mapply(function(fname, mean.img, sd.img) {
      img = ants_apply_transforms(
        fixed = template_fname,
        moving = fname,
        transformlist = t1_reg$fwdtransforms,
        interpolator = interpolator)
      z.img = (img - mean.img) / sd.img
      z.img[ !is.finite(z.img)] = NA
      res_z = ants_apply_transforms(
        fixed = t1, # reversed!
        moving = z.img,
        transformlist = t1_reg$invtransforms,
        interpolator = interpolator)
      z.img
    }, x, mean_imgs, sd_imgs)

    if (remask) {
      if (!is.null(mask)) {
        mask = check_nifti(mask)
        norm = lapply(norm, mask_img, mask = mask)
      }
    }

    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, norm, fnames)
  }

  norm = lapply(
    fnames,
    identity)

  return(norm)
}