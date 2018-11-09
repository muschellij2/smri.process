#' SMRI Pipeline for Preprocessing Prior to Normalization
#'
#' @param x List of images
#' @param gold_standard Gold Standard image/filename, if applicable
#' @param gs_space space the Gold Standard is located
#' @param probs passed to \code{\link{winsor}} for Winsorization
#' @param non_zero Should zeroes be excluded from the calculation
#' of quantiles? Passed to \code{\link{winsor}}.
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#' @param interpolator interpolation passed to \code{\link{reg_to_t1}}
#' @param brain_mask Brain Mask in \code{reg_space}
#' @param gs_interpolator Gold Standard interpolation
#' passed to \code{\link{reg_to_t1}}
#' @param num_templates Number of templates to use for MALF
#' @param malf_transform type of registration transformation for MALF
#' @param reg_space space to register images to
#' @param remove_negative After reorientation, N4, and registration,
#' any values < 0 are set to 0
#' @param zero_origin Should the origin be set to 0 for the
#' image to be registered to?
#' @param brain_threshold Percentage threshold for the brain MALF image
#' to be considered part of the mask
#' @param brain_malf_function Function to be passed to \code{\link{malf}},
#' and subsequently \code{\link{stat_img}} for brain extraction
#' @param cleanup Argument to \code{\link{reduce_img}} for reducing the images.
#' @param reduce If \code{TRUE}, the image size will be reduced using
#' \code{\link{reduce_img}}.
#' @param brain_extraction_method Brain extraction method, either
#' \code{\link{malf}} or \code{abp} for \code{\link{n4_skull_strip}}
#' @param ... Additional arguments to MALF
#'
#' @return List of the images, brain mask, suffix, and output directory
#' @export
#'
#' @importFrom extrantsr registration malf
#' @importFrom neurobase readnii check_nifti write_nifti
smri_prenormalize = function(
  x,
  gold_standard = NULL,
  gs_space = NULL,
  probs = c(0, 0.995),
  non_zero = TRUE,
  interpolator = "lanczosWindowedSinc",
  brain_mask = NULL,
  gs_interpolator = "genericLabel",
  num_templates = 15,
  malf_transform = "SyNAggro",
  outdir = tempdir(),
  reg_space = "T1",
  brain_extraction_method = c("abp", "malf"),
  brain_malf_function = "staple_prob",
  brain_threshold = 0.5,
  verbose = TRUE,
  remove_negative = TRUE,
  zero_origin = TRUE,
  reduce = TRUE,
  cleanup =1,
  ...
) {


  proc = bc_noneck_reduce(
    x = x,
    remove_negative = remove_negative,
    gold_standard = gold_standard,
    gs_space = gs_space,
    outdir = outdir,
    verbose = verbose,
    probs = probs,
    non_zero = non_zero,
    cleanup = cleanup,
    reduce = reduce)

  suffix = proc$suffix
  suffix = paste0(suffix, "_regto", reg_space)
  reg = reg_to_t1(
    x = proc$images,
    gs_space = gs_space,
    interpolator = interpolator,
    gs_interpolator = gs_interpolator,
    outdir = outdir,
    verbose = verbose,
    reg_space = reg_space,
    suffix = suffix,
    remove_negative = remove_negative,
    zero_origin = zero_origin)

  reg = reg$images
  rigid_registrations = reg$registrations

  gold_standard = reg$GOLD_STANDARD
  reg$GOLD_STANDARD = NULL
  gs_suffix = suffix

  inverted = NULL

  if (!is.null(brain_mask)) {
    brain_mask = check_nifti(brain_mask)
    malf_result = NULL
  } else {
    brain_extraction_method = match.arg(brain_extraction_method)
    brain_mask_file = file.path(
      outdir,
      "Brain_Mask.nii.gz")
    if (brain_extraction_method == "abp") {
      malf_result = NULL
      if (all_exists(brain_mask_file)) {
        # warning("Using brain mask file in outdir")
        brain_mask = readnii(brain_mask_file)
      } else {
        brain_mask = n4_skull_strip(
          file = reg$T1,
          template = penn115_image_fname(),
          template_mask = penn115_brain_mask_fname(),
          verbose = verbose,
          n_iter = 3)
        write_nifti(brain_mask, filename = brain_mask_file)
      }
    } else {
      brain_pct_file = file.path(
        outdir,
        "Brain_Percentages.nii.gz")
      if (all_exists(brain_mask_file, brain_pct_file)) {
        # warning("Using brain mask file in outdir")
        brain_mask = readnii(brain_mask_file)
        malf_result = readnii(brain_pct_file)
      } else {
        ind = seq(num_templates)
        templates = malf.templates::malf_images()
        images = templates$images[ind]
        masks = templates$masks[ind]

        if (verbose > 0) {
          msg = paste0(
            "Running MALF for brain mask with ", num_templates,
            " templates - this may take some time")
          message(msg)
        }
        args = list(...)
        inverted = args$inverted
        keep_regs = args$keep_regs

        malf_result = malf(
          infile = reg$T1,
          template.images = images,
          template.structs = masks,
          interpolator = "Linear",
          other_interpolator = "genericLabel",
          invert_interpolator = "genericLabel",
          typeofTransform = malf_transform,
          verbose = verbose,
          # func = "mode",
          func = brain_malf_function,
          retimg = TRUE,
          outfile = brain_pct_file,
          ...
        )
        if (!is.null(keep_regs)) {
          if (keep_regs) {
            malf_result = malf_result$outimg
          }
        }
        brain_mask = malf_result >= brain_threshold
        writenii(brain_mask, filename = brain_mask_file)
      }
    }
  }

  suffix = paste0(suffix, "_brain")
  brains = apply_mask(
    x = reg,
    mask = brain_mask,
    outdir = outdir,
    verbose = verbose,
    suffix = suffix)

  suffix = paste0(suffix, "_N4")
  n4_brains = n4_raw(
    x = brains,
    verbose = verbose,
    mask = brain_mask,
    outdir = outdir,
    suffix = suffix)

  L = list(
    images = n4_brains,
    intermediate = list(
      reduced = proc$images,
      registered = reg,
      masked = brains),
    brain_mask = brain_mask,
    suffix = suffix,
    gs_suffix = gs_suffix,
    outdir = outdir,
    rigid_registrations = rigid_registrations,
    num_templates = num_templates,
    brain_threshold = brain_threshold,
    brain_malf_function = brain_malf_function
  )
  L$brain_pct = malf_result
  L$inverted = inverted
  L$GOLD_STANDARD = gold_standard

  return(L)
}

