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
  malf_transform = "SynAggro"
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
    interpolator = "Linear",
    gs_interpolator = "NearestNeighbor",
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