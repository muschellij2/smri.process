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
#' @param tissue_suffix Suffix to add onto tissue segmentation in addition
#' to the suffix for normalization
#' @param norm_outdir Output directory of spatially normalized data.  The
#' output directory of the segmentation images will still be specified
#' in \code{prenormalize$outdir}
#' @param ... arguments passed to \code{\link{t1_segment}}
#'
#' @return List of images, suffix, brain mask, gold standard,
#' and registration if applicable
#' @export
seg_normalize = function(
  prenormalize,
  template = c("none", "Eve", "MNI"),
  norm_outdir = prenormalize$outdir,
  verbose = TRUE,
  typeofTransform = "SyN",
  interpolator = "lanczosWindowedSinc",
  dis_interpolator = "genericLabel",
  tissue_suffix = "",
  ...
) {

  inverted = prenormalize$inverted
  args = list()
  if (is.null(inverted)) {
    args = list(...)
    inverted = args$inverted
    if (is.null(inverted)) {
      inverted = FALSE
    }
  }
  args$t1 = prenormalize$images$T1
  args$outdir = prenormalize$outdir
  args$num_templates = prenormalize$num_templates
  args$inverted = inverted
  args$verbose = verbose

  tissue = do.call(t1_segment, args = args)

  # tissue = t1_segment(
  #   t1 = prenormalize$images$T1,
  #   outdir = prenormalize$outdir,
  #   num_templates = prenormalize$num_templates,
  #   inverted = inverted,
  #   ...)

  fast_res = multi_fast(
    x = prenormalize$images,
    bias_correct = FALSE,
    outdir = prenormalize$outdir,
    verbose = verbose)

  if (is.null(norm_outdir)) {
    norm_outdir = prenormalize$outdir
  }
  prenormalize$outdir = norm_outdir

  resampled = spatial_normalize(
    prenormalize,
    verbose = verbose,
    template = template,
    typeofTransform = typeofTransform,
    interpolator = interpolator
  )

  if (verbose > 0) {
    msg = "Applying to Tissues/Structures"
    message(msg)
  }
  resampled_hard = apply_spatial_normalize(
    x = tissue[c("TISSUES", "STRUCTURES")],
    template = resampled$template,
    template_fname = resampled$template_fname,
    fwdtransforms = resampled$reg_to_template$fwdtransforms,
    suffix = tissue_suffix,
    interpolator = dis_interpolator,
    outdir = resampled$outdir,
    verbose = verbose
  )

  if (verbose > 0) {
    msg = "Applying to Tissue Probabilities"
    message(msg)
  }
  resampled_probs = apply_spatial_normalize(
    x = tissue[setdiff(names(tissue), c("TISSUES", "STRUCTURES"))],
    template = resampled$template,
    template_fname = resampled$template_fname,
    fwdtransforms = resampled$reg_to_template$fwdtransforms,
    suffix = tissue_suffix,
    interpolator = resampled$interpolator,
    outdir = resampled$outdir,
    verbose = verbose
  )

  resampled_tissue = c(
    resampled_hard,
    resampled_probs)

  if (verbose > 0) {
    msg = "Applying to FAST Output"
    message(msg)
  }
  resampled_fast = lapply(
    fast_res,
    function(x) {
      names(x) = nii.stub(x, bn = TRUE)
      seg = grepl("seg|mixeltype", names(x))

      dis_data = apply_spatial_normalize(
        x = x[seg],
        template = resampled$template,
        template_fname = resampled$template_fname,
        fwdtransforms = resampled$reg_to_template$fwdtransforms,
        suffix = tissue_suffix,
        interpolator = dis_interpolator,
        outdir = resampled$outdir,
        verbose = verbose)
      dis_data = unlist(dis_data)

      con_data = apply_spatial_normalize(
        x = x[!seg],
        template = resampled$template,
        template_fname = resampled$template_fname,
        fwdtransforms = resampled$reg_to_template$fwdtransforms,
        suffix = tissue_suffix,
        interpolator = resampled$interpolator,
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