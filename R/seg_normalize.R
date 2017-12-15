#' Spatially Normalize SMRI data and Tissue Segmentation
#'
#' @param prenormalize list from \code{\link{smri_prenormalize}}
#' @param template Resampling to \code{1x1x1} (template = none) or
#' Template to register to - either MNI or Eve template
#' @param verbose print diagnostic messages
#' @param typeofTransform Transformation to use for registration,
#' passed to \code{\link{registration}}
#' @param interpolator Interpolation done, passed to
#' \code{\link{registration}} for continuous data
#' @param dis_interpolator Interpolation done, passed to
#' \code{\link{registration}} for discrete data
#'
#' @return List of images, suffix, brain mask, gold standard,
#' and registration if applicable
#' @export
seg_normalize = function(
  prenormalize,
  template = c("none", "Eve", "MNI"),
  verbose = TRUE,
  typeofTransform = "SyN",
  interpolator = "lanczosWindowedSinc",
  dis_interpolator = "genericLabel",
  tissue_suffix = ""
) {

  tissue = t1_segment(
    t1 = prenormalize$images$T1,
    outdir = prenormalize$outdir,
    num_templates = prenormalize$num_templates)

  fast_res = multi_fast(
    x = prenormalize$images,
    bias_correct = FALSE,
    outdir = prenormalize$outdir,
    verbose = verbose)

  resampled = spatial_normalize(
    prenormalize,
    verbose = verbose,
    template = template,
    typeofTransform = typeofTransform,
    interpolator = interpolator
  )

  resampled_hard = apply_spatial_normalize(
    x = tissue[c("TISSUES", "STRUCTURES")],
    template = resampled$template,
    template_fname = resampled$template_fname,
    fwdtransforms = resampled$fwdtransforms,
    suffix = tissue_suffix,
    interpolator = dis_interpolator,
    outdir = resampled$outdir,
    verbose = verbose
  )

  resampled_probs = apply_spatial_normalize(
    x = tissue[setdiff(names(tissue), c("TISSUES", "STRUCTURES"))],
    template = resampled$template,
    template_fname = resampled$template_fname,
    fwdtransforms = resampled$fwdtransforms,
    suffix = tissue_suffix,
    interpolator = resampled$interpolator,
    outdir = resampled$outdir,
    verbose = verbose
  )

  resampled_tissue = c(
    resampled_hard,
    resampled_probs)

  resampled_fast = lapply(
    fast_res,
    function(x) {
      names(x) = nii.stub(x, bn = TRUE)
      seg = grepl("seg|mixeltype", names(x))

      dis_data = apply_spatial_normalize(
        x = x[seg],
        template = resampled$template,
        template_fname = resampled$template_fname,
        fwdtransforms = resampled$fwdtransforms,
        suffix = tissue_suffix,
        interpolator = dis_interpolator,
        outdir = resampled$outdir,
        verbose = verbose)
      dis_data = unlist(dis_data)

      con_data = apply_spatial_normalize(
        x = x[!seg],
        template = resampled$template,
        template_fname = resampled$template_fname,
        fwdtransforms = resampled$fwdtransforms,
        suffix = tissue_suffix,
        interpolator = "nearestNeighbor",
        outdir = resampled$outdir,
        verbose = verbose)
      con_data = unlist(con_data)
      x = c(dis_data, con_data)
      x = x[names(x)]
      x
    })

  resampled$tissue = resampled_tissue
  resampled$fast = resampled_fast
  prenormalize$tissue = tissue
  prenormalize$fast = fast

  L = list(
    normalized = resampled,
    native = prenormalize
  )
  return(L)
}