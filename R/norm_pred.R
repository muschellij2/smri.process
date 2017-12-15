#' SMRI Intensity Normalize and Create Predictors
#'
#' @param normalized The \code{normalized} element from
#' \code{\link{seg_normalize}} output
#' @param outdir Output directory
#'
#' @param normalization Normalization Method
#' @param remask should mask_img be applied to the output images
#' after normalization
#' @param radius radius of neighborhood.  See \code{\link{smri_corr_img}}
#' and \code{\link{smri_moments}}
#' @param sigma vector of sigmas for smoothers in FSL see
#' \code{\link{fsl_smooth_images}}
#' @param suffix what the append to end of the output filename
#' @param ... Additional arguments to \code{normaliztaion} function
#' @param verbose print diagnostic messages
#'
#' @return List of images
#' @export
#'
#' @importFrom neurobase zscore_img quantile_img
#' @importFrom WhiteStripe whitestripe_norm whitestripe whitestripe_hybrid
norm_predictors = function(
  normalized,
  normalization = c("trimmed_z", "z", "quantile", "minmax",
                    "whitestripe_T1", "whitestripe_T2",
                    "whitestripe_hybrid"),
  remask = TRUE,
  suffix = normalized$suffix,
  outdir = normalized$outdir,
  verbose = TRUE,
  radius = rep(1, 3),
  sigma = c(3, 5, 10, 20),
  ...
) {

  normalization = match.arg(normalization)
  n_suffix = gsub("_", "", normalization)
  xsuffix = suffix
  suffix = paste0(suffix, "_", n_suffix)

  norm = intensity_normalize(
    x = normalized$images,
    mask = normalized$brain_mask,
    normalization = normalization,
    outdir = normalized$outdir,
    remask = remask,
    suffix = suffix,
    verbose = verbose,
    ...
  )

  mom_norm = smri_moments(
    x = norm,
    radius = ,
    mask = normalized$brain_mask,
    outdir = normalized$outdir,
    suffix = paste0(suffix, "_"),
    verbose = verbose
  )

  mom_tissue = smri_moments(
    x = normalized$tissue[ c("CSF", "GM", "WM")],
    radius = radius,
    mask = normalized$brain_mask,
    outdir = normalized$outdir,
    suffix = "_",
    verbose = verbose
  )

  corr = smri_corr_img(
    x = norm,
    radius = radius,
    mask = normalized$brain_mask,
    outdir = normalized$outdir,
    suffix = paste0(suffix, "_corr"),
    verbose = verbose
  )

  fsl_smoothed = fsl_smooth_images(
    x = norm,
    sigma = sigma,
    suffix = paste0(suffix, "_fslsmooth"),
    outdir = normalized$outdir,
    mask = normalized$brain_mask,
    verbose = verbose)


  suffix = paste0(xsuffix, "_quantile")

  if (normalization != "quantile") {
    qnorm = intensity_normalize(
      x = normalized$images,
      mask = normalized$brain_mask,
      normalization = "quantile",
      outdir = normalized$outdir,
      suffix = suffix,
      verbose = verbose
    )
    names(qnorm) = paste0(names(qnorm), "_quantile")
  } else {
    qnorm = NULL
  }

  names(fsl_smoothed) = paste0("sigma_", names(fsl_smoothed))

  L = list(
    normalized = normalized,
    intensity_normalized = norm,
    normalized_smoothed = fsl_smoothed,
    normalized_moments = mom_norm,
    tissue_prob_moments = mom_tissue,
    correlation_images = corr)
  L$quantile_images = qnorm

  return(L)

}

#' Gather all predictors from norm_predictors function
#'
#' @param x Output from \code{\link{norm_predictors}}
#'
#' @return A vector of filenames.  Note hard tissue segmentations and
#' sub-structural segmentations are removed
#' @export
gather_predictors = function(x) {
  # these aren't predictors for now
  x$normalized$tissue$STRUCTURES = NULL
  x$normalized$tissue$TISSUES = NULL
  x$normalized$fast = lapply(
    x$normalized$fast,
    function(x) {
      x[!grepl("seg|mixeltype", x)]
    })

  fnames =
    c(x$intensity_normalized,
      x$normalized$fast,
      x$normalized$tissue,
      x$normalized_smoothed,
      x$normalized_moments,
      x$tissue_prob_moments,
      x$correlation_images,
      x$quantile_images)
  fnames = unlist(fnames)
  return(fnames)
}