
#' SMRI Pipeline
#'
#' @param x List of images
#' @param gold_standard Gold Standard image/filename, if applicable
#' @param gs_space space the Gold Standard is located
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
#' @return
#' @export
#'
#' @examples
#' @importFrom extrantsr registration malf
#' @importFrom neurobase readnii check_nifti
smri_pipeline = function(
  x,
  gold_standard = NULL,
  gs_space = NULL,
  outdir = tempdir(),
  verbose = TRUE,
  interpolator = "Linear",
  brain_mask = NULL,
  gs_interpolator = "NearestNeighbor",
  num_templates = 15,
  malf_transform = "SynAggro",
  ...
) {


  proc = bc_noneck_reduce(
    x = x,
    gold_standard = gold_standard,
    gs_space = gs_space,
    outdir = outdir,
    verbose = verbose)

  reg = reg_to_t1(
    x = proc$images,
    gs_space = gs_space,
    interpolator = interpolator,
    gs_interpolator = gs_interpolator,
    outdir = outdir,
    verbose = verbose)


  if (!is.null(brain_mask)) {
    brain_mask = check_nifti(brain_mask)
  } else {
    brain_mask_file = file.path(
      outdir,
      "Brain_Mask.nii.gz")
    if (file.exists(brain_mask_file)) {
      warning("Using brain mask file in outdir")
      brain_mask = readnii(brain_mask)
    } else {
      ind = seq(num_templates)
      templates = malf.templates::malf_images()
      images = templates$images[ind]
      masks = templates$masks[ind]

      malf_result = malf(
        infile = reg$T1,
        template.images = images,
        template.structs = masks,
        interpolator = "genericLabel",
        typeofTransform = malf_transform,
        verbose = verbose,
        func = "mode",
        retimg = FALSE,
        ...
      )
      brain_mask = malf_result$outimg
    }

  }


}