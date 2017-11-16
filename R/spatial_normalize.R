#' Spatially Normalize SMRI data
#'
#' @param prenormalize list from \code{\link{smri_prenormalize}}
#' @param template Resampling to \code{1x1x1} (template = none) or
#' Template to register to - either MNI or Eve template
#' @param verbose print diagnostic messages
#' @param typeofTransform Transformation to use for registration,
#' passed to \code{\link{registration}}
#' @param interpolator Interpolation done, passed to
#' \code{\link{registration}}
#'
#' @return List of images, suffix, brain mask, gold standard,
#' and registration if applicable
#' @export
#' @importFrom EveTemplate getEvePath
#' @importFrom MNITemplate getMNIPath
#' @importFrom extrantsr registration resample_to_target resample_image
spatial_normalize = function(
  prenormalize,
  template = c("none", "Eve", "MNI"),
  verbose = TRUE,
  typeofTransform = "SyN",
  interpolator = "lanczosWindowedSinc"
) {

  template = match.arg(template)
  native = template == "none"

  if (native) {
    app = "_resampled"
  } else {
    app = paste0("_regto", template)
  }

  suffix = prenormalize$suffix
  gold_standard = prenormalize$GOLD_STANDARD
  brain_mask = prenormalize$brain_mask
  brain_pct = prenormalize$brain_pct
  outdir = prenormalize$outdir

  suffix = paste0(suffix, app)
  nii_names = names(prenormalize$images)

  fnames = file.path(
    outdir,
    paste0(nii_names,
           suffix,
           ".nii.gz"))
  names(fnames) = nii_names

  if (!is.null(brain_mask)) {
    brain_fname = file.path(
      outdir,
      paste0("Brain_Mask",
            app,
             ".nii.gz")
    )
  } else {
    brain_fname = NULL
  }

  if (!is.null(brain_pct)) {
    brain_pct_fname = file.path(
      outdir,
      paste0("Brain_Percentages",
             app,
             ".nii.gz")
    )
    brain_pct_fname = c(BRAIN_PCT = brain_pct_fname)
  } else {
    brain_pct_fname = NULL
  }

  if (!is.null(gold_standard)) {
    gs_fname = file.path(
      outdir,
      paste0("GOLD_STANDARD",
             prenormalize$gs_suffix,
             app,
             ".nii.gz")
    )
    gs_fname = c(GOLD_STANDARD = gs_fname)
  } else {
    gs_fname = NULL
  }
  all_fnames = c(fnames, BRAIN_MASK = brain_fname,
                 gs_fname, brain_pct_fname)


  # typeofTransform = "SyN"
  # interpolator = "lanczosWindowedSinc"
  if (!all_exists(all_fnames)) {
    if (native) {
      if (verbose > 0) {
        message("Resampling Image")
      }
      # resample to 1x1x1 for
      resampled = lapply(
        prenormalize$images,
        resample_image,
        parameters = c(1, 1, 1),
        parameter_type = "mm",
        interpolator = "windowedsinc")

      if (!is.null(brain_mask)) {
        resampled_brain_mask = resample_image(
          prenormalize$brain_mask,
          parameters = c(1, 1, 1),
          parameter_type = "mm",
          interpolator = "nearestneighbor")
        writenii(resampled_brain_mask, filename = brain_fname)

        resampled = lapply(
          resampled,
          mask_img,
          mask = resampled_brain_mask)
      }

      if (!is.null(brain_pct)) {
        resampled_brain_pct = resample_image(
          prenormalize$brain_pct,
          parameters = c(1, 1, 1),
          parameter_type = "mm",
          interpolator = "windowedsinc")
        writenii(resampled_brain_pct, filename = brain_pct_fname)
      }

      # this should be after remasking
      mapply(function(img, fname){
        writenii(img, filename = fname)
      }, resampled, fnames)

      if (!is.null(gold_standard)) {
        # resampled_gs = resample_image(
        #   prenormalize$GOLD_STANDARD,
        #   parameters = c(1, 1, 1),
        #   parameter_type = "mm",
        #   interpolator = "nearestneighbor")
        # using reample_to_target so that you can use
        # genericLabel interpolator
        resampled_gs = resample_to_target(
          prenormalize$GOLD_STANDARD,
          target = prenormalize$images[[1]],
          verbose = verbose > 1,
          interpolator = "genericLabel")
        writenii(resampled_gs, filename = gs_fname)
      }

      t1_reg = NULL
      template_fname = prenormalize$images[[1]]

    } else {

      template_fname = switch(
        template,
        Eve = EveTemplate::getEvePath(what = "Brain"),
        MNI = MNITemplate::getMNIPath(what = "Brain", res = "1mm")
      )
      # template_img = readnii(template_fname)

      if (verbose > 0) {
        message(paste0(
          "Registering T1 to ", template, " Template")
        )
      }
      t1_reg = registration(
        filename = prenormalize$images$T1,
        template.file = template_fname,
        skull_strip = FALSE,
        correct = FALSE,
        outprefix = nii.stub(prenormalize$images$T1),
        verbose = verbose > 1,
        typeofTransform = typeofTransform,
        interpolator = interpolator)

      resampled = lapply(
        prenormalize$images,
        function(r) {
          ants_apply_transforms(
            moving = r,
            fixed = template_fname,
            transformlist = t1_reg$fwdtransforms,
            interpolator = interpolator)
        })

      if (!is.null(brain_mask)) {
        resampled_brain_mask = ants_apply_transforms(
          moving = prenormalize$brain_mask,
          fixed = template_fname,
          transformlist = t1_reg$fwdtransforms,
          interpolator = "nearestNeighbor")

        writenii(resampled_brain_mask, filename = brain_fname)

        resampled = lapply(
          resampled,
          mask_img,
          mask = resampled_brain_mask)

      }

      if (!is.null(brain_pct)) {
        resampled_brain_pct = ants_apply_transforms(
          moving = prenormalize$brain_pct,
          fixed = template_fname,
          transformlist = t1_reg$fwdtransforms,
          interpolator = interpolator)
        writenii(resampled_brain_pct, filename = brain_pct_fname)
      }

      # this should be after remasking
      mapply(function(img, fname){
        writenii(img, filename = fname)
      }, resampled, fnames)

      if (!is.null(gold_standard)) {
        resampled_gs = ants_apply_transforms(
          moving = prenormalize$GOLD_STANDARD,
          fixed = template_fname,
          transformlist = t1_reg$fwdtransforms,
          interpolator = "genericLabel")
        writenii(resampled_gs, filename = gs_fname)
      }

    }
    t1_reg = t1_reg[c("fwdtransforms",
                      "invtransforms",
                      "typeofTransform",
                      "interpolator")]
  } else {
    if (native) {
      t1_reg = NULL
    } else {
      outprefix = nii.stub(prenormalize$images$T1)
      t1_reg = list(
        fwdtransforms = paste0(
          outprefix,
          c("1Warp.nii.gz",
          "0GenericAffine.mat")),
        invtransforms = paste0(
          outprefix,
          c("0GenericAffine.mat",
          "1InverseWarp.nii.gz")
        ),
        typeofTransform = typeofTransform,
        interpolator = interpolator)
    }
  }



  resampled = lapply(fnames, identity)
  resampled_brain_mask = brain_fname
  resampled_brain_pct = brain_pct_fname
  resampled_gs = gs_fname

  L = list(
    images = resampled,
    suffix = suffix,
    gs_suffix = paste0(prenormalize$gs_suffix, app),
    outdir = outdir,
    template = template
  )
  L$template_fname = template_fname
  L$GOLD_STANDARD = gold_standard
  L$brain_mask = resampled_brain_mask
  L$reg_to_template = t1_reg
  L$brain_pct = resampled_brain_pct

  return(L)
}