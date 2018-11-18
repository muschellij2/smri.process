#' Spatially Normalize SMRI data and Tissue Segmentation
#'
#' Spatially Normalize SMRI data and Tissue Segmentation, note T1 is
#' required
#'
#' @param prenormalize list from \code{\link{smri_prenormalize}}
#' @param template Resampling to \code{1x1x1} (template = none) or
#' Template to register to - either MNI or Eve template
#' @param verbose print diagnostic messages
#' @param typeofTransform Transformation to use for registration,
#' passed to \code{\link{registration}}
#' @param segment_typeofTransform Transformation to use for
#' \code{malf} segmentation
#' @param interpolator Interpolation done, passed to
#' \code{\link{registration}} for continuous data
#' @param malf_interpolator interpolator for MALF procedure,
#' passed to \code{\link{t1_segment}}
#' @param dis_interpolator Interpolation done, passed to
#' \code{\link{registration}} for discrete data
#' @param norm_outdir Output directory of spatially normalized data.  The
#' output directory of the segmentation images will still be specified
#' in \code{prenormalize$outdir}
#' @param ... arguments passed to \code{\link{t1_segment}}
#' @param force_registration Should the MALF registrations be re-run if
#' they already exist?
#' @param copy_origin Copy image origin from image
#' being registered, using \code{\link{antsCopyOrigin}}
#'
#' @return List of images, suffix, brain mask, gold standard,
#' and registration if applicable
#' @export
# #' @param tissue_suffix Suffix to add onto tissue segmentation in addition
# #' to the suffix for normalization
seg_normalize = function(
  prenormalize,
  template = c("none", "Eve", "MNI"),
  norm_outdir = prenormalize$outdir,
  verbose = TRUE,
  typeofTransform = "SyN",
  segment_typeofTransform = typeofTransform,
  interpolator = "lanczosWindowedSinc",
  dis_interpolator = "genericLabel",
  malf_interpolator = "genericLabel",
  copy_origin = TRUE,
  force_registration = FALSE,
  ...
) {

  inverted = prenormalize$inverted
  args = list(...)
  if (is.null(inverted)) {
    inverted = args$inverted
    if (is.null(inverted)) {
      inverted = FALSE
    }
  }
  if (!"outdir" %in% names(prenormalize)) {
    stop("prenormalize$outdir must not be NULL")
  }
  if (!"num_templates" %in% names(prenormalize)) {
    stop("prenormalize$num_templates must not be NULL")
  }
  if (!"images" %in% names(prenormalize)) {
    stop("prenormalize$images must not be NULL")
  }
  imgs = prenormalize$images
  if (!"T1" %in% names(imgs)) {
    stop("prenormalize$images must have T1 element")
  }

  args$t1 = prenormalize$images$T1
  args$outdir = prenormalize$outdir
  args$num_templates = prenormalize$num_templates
  args$inverted = inverted
  args$interpolator = malf_interpolator
  args$verbose = verbose
  args$typeofTransform = segment_typeofTransform
  args$force_registration = force_registration

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
    interpolator = interpolator,
    dis_interpolator = dis_interpolator,
    copy_origin = copy_origin
  )

  if (verbose > 0) {
    msg = "Applying to Tissues/Structures"
    message(msg)
  }
  t1 = prenormalize$images$T1
  if (!is.null(t1)) {
    t1 = check_ants(t1)
    origin = antsGetOrigin(t1)
    rm(t1)
  } else {
    origin = NULL
  }
  apm = function(x, interpolator = dis_interpolator) {
    apply_spatial_normalize(
      x = x,
      template = resampled$template,
      template_fname = resampled$template_fname,
      fwdtransforms = resampled$reg_to_template$fwdtransforms,
      suffix = "", # I think should work
      interpolator = interpolator,
      outdir = resampled$outdir,
      verbose = verbose,
      copy_origin = copy_origin,
      origin = origin
    )
  }
  resampled_hard = apm(
    x = tissue[c("TISSUES", "STRUCTURES")],
    interpolator = dis_interpolator)

  if (verbose > 0) {
    msg = "Applying to Tissue Probabilities"
    message(msg)
  }
  resampled_probs =  apm(
    x = tissue[setdiff(names(tissue), c("TISSUES", "STRUCTURES"))],
    interpolator = resampled$interpolator
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

      if (any(seg)) {
        dis_data = apm(
          x = x[seg],
          interpolator = dis_interpolator
        )
        dis_data = unlist(dis_data)
      } else {
        dis_data = NULL
      }

      if (any(!seg)) {
        con_data = apm(
          x = x[!seg],
          interpolator = resampled$interpolator
        )
        con_data = unlist(con_data)
      } else {
        con_data = NULL
      }
      x = c(dis_data, con_data)
      x = x[names(x)]
      x
    })

  resampled$tissue = resampled_tissue
  resampled$fast = resampled_fast
  prenormalize$tissue = tissue
  prenormalize$fast = fast_res

  L = list(
    normalized = resampled,
    native = prenormalize
  )

  return(L)
}