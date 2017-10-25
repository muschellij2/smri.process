
#' Normalize using Min/Max
#' @param img Image to normalize
#' @param mask Mask same dimensions as image
#' @param remask Should the image be remasked after normalizing?
#'
#' @return Object of same class as img
#' @export
#' @importFrom neurobase finite_img mask_vals
minmax = function(img, mask = NULL, remask = TRUE) {
  img = check_nifti(img, allow.array = TRUE)

  dimg = dim(img)
  if (is.null(mask)) {
    mask = array(1, dim = dimg)
  }
  vals = mask_vals(img, mask = mask)
  r = range(vals)
  img = (img - r[1]) / (r[2] - r[1])
  img = finite_img(img)
  if (remask) {
    img = mask_img(img, mask)
  }
  return(img)
}

#' SMRI Normalize
#'
#' @param x List of images
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#' @param mask mask to normalize over
#' @param ... Additional arguments to MALF
#'
#' @return List of images
#' @export
#'
#' @importFrom neurobase zscore_img quantile_img
smri_normalize = function(
  x,
  mask = NULL,
  normalization = c("z", "trimmed_z", "quantile", "minmax"),
  outdir = tempdir(),
  verbose = TRUE,
  remask = TRUE,
  ...
) {

  args = list(...)

  normalization = match.arg(normalization)
  if (normalization %in% c("z", "trimmed_z")) {
    if (normalization == "z") {
      args$trim = 0
    }
    args$mask = mask
    args$centrality = "trimmed_mean"
    args$variability = "trimmed_sd"
    norm = lapply(x, function(r) {
      args$img = r
      do.call("zscore_img", args = args)
    })
  }

  if (normalization %in% "quantile") {
    args$mask = mask
    norm = lapply(x, function(r) {
      args$img = r
      do.call("quantile_img", args = args)
    })
  }

  if (normalization %in% c("minmax")) {
    args$mask = mask
    norm = lapply(x, function(r) {
      args$img = r
      do.call("minmax", args = args)
    })
  }

  if (remask) {
    if (!is.null(mask)) {
      norm = lapply(norm, mask_img, mask = mask)
    }
  }


  return(norm)
}
