
#' SMRI Pipeline for Preprocessing Prior to Normalization
#'
#' @param x List of images
#' @param gold_standard Gold Standard image/filename, if applicable
#' @param gs_space space the Gold Standard is located
#' @param probs passed to \code{\link{winsor}} for Winsorization
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#' @param interpolator interpolation passed to \code{\link{reg_to_t1}}
#' @param brain_mask Brain Mask in T1 Space
#' @param gs_interpolator Gold Standard interpolation
#' passed to \code{\link{reg_to_t1}}
#' @param num_templates Number of templates to use for MALF
#' @param malf_transform type of registration transformation for MALF
#' @param ... Additional arguments to MALF
#'
#' @return List of the images, brain mask, suffix, and output directory
#' @export
#'
#' @importFrom extrantsr registration malf
#' @importFrom neurobase readnii check_nifti
smri_prenormalize = function(
  x,
  gold_standard = NULL,
  gs_space = NULL,
  probs = c(0, 0.999),
  interpolator = "lanczosWindowedSinc",
  brain_mask = NULL,
  gs_interpolator = "NearestNeighbor",
  num_templates = 15,
  malf_transform = "SyNAggro",
  outdir = tempdir(),
  verbose = TRUE,
  ...
) {


  proc = bc_noneck_reduce(
    x = x,
    gold_standard = gold_standard,
    gs_space = gs_space,
    outdir = outdir,
    verbose = verbose,
    probs = probs)

  suffix = proc$suffix
  suffix = paste0(suffix, "_regtoT1")
  reg = reg_to_t1(
    x = proc$images,
    gs_space = gs_space,
    interpolator = interpolator,
    gs_interpolator = gs_interpolator,
    outdir = outdir,
    verbose = verbose,
    suffix = suffix)

  reg = reg$images
  rigid_registrations = reg$registrations

  gold_standard = reg$GOLD_STANDARD
  reg$GOLD_STANDARD = NULL
  gs_suffix = suffix

  if (!is.null(brain_mask)) {
    brain_mask = check_nifti(brain_mask)
    malf_result = NULL
  } else {
    brain_mask_file = file.path(
      outdir,
      "Brain_Mask.nii.gz")
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
      malf_result = malf(
        infile = reg$T1,
        template.images = images,
        template.structs = masks,
        interpolator = "genericLabel",
        typeofTransform = malf_transform,
        verbose = verbose > 1,
        # func = "mode",
        func = "pct",
        retimg = TRUE,
        outfile = brain_pct_file
        #, ...
      )
      brain_mask = malf_result >= 0.5
      writenii(brain_mask, filename = brain_mask_file)
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
    rigid_registrations = rigid_registrations
  )
  L$brain_pct = malf_result
  L$GOLD_STANDARD = gold_standard

  return(L)
}

