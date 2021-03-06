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
#' @param dis_interpolator Interpolation done, passed to
#' \code{\link{registration}} for discrete data
#' @param copy_origin Copy image origin from image
#' being registered, using \code{\link{antsCopyOrigin}}
#' @return List of images, suffix, brain mask, gold standard,
#' and registration if applicable
#' @export
#' @importFrom EveTemplate getEvePath
#' @importFrom MNITemplate getMNIPath
#' @importFrom extrantsr registration resample_to_target resample_image
#' @importFrom extrantsr antsCopyOrigin transformlist_from_outprefix
spatial_normalize = function(
  prenormalize,
  template = c("none", "Eve", "MNI"),
  verbose = TRUE,
  typeofTransform = "SyN",
  interpolator = "lanczosWindowedSinc",
  dis_interpolator = "genericLabel",
  copy_origin = TRUE
) {

  template = match.arg(template)
  native = template == "none"

  if (native) {
    app = "_resampled"
  } else {
    app = paste0("_regto", template)
  }

  if (!"outdir" %in% names(prenormalize)) {
    stop("prenormalize$outdir must not be NULL")
  }
  if (!"images" %in% names(prenormalize)) {
    stop("prenormalize$images must not be NULL")
  }
  imgs = prenormalize$images
  if (!"T1" %in% names(imgs)) {
    stop("prenormalize$images must have T1 element")
  }

  if (!"suffix" %in% names(prenormalize)) {
    warning(paste0("prenormalize$suffix is NULL",
                   ", likely to cause naming problems"))
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
    if (!"gs_suffix" %in% names(prenormalize)) {
      warning(paste0("prenormalize$gs_suffix is NULL",
                     ", likely to cause naming problems"))
    }
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
        if (verbose > 1) {
          message("Resampling Brain Mask")
        }
        resampled_brain_mask = resample_to_target(
          prenormalize$brain_mask,
          target = resampled[[1]],
          verbose = verbose > 1,
          interpolator = dis_interpolator,
          copy_origin = copy_origin)
        writenii(resampled_brain_mask, filename = brain_fname)

        resampled = lapply(
          resampled,
          mask_img,
          mask = resampled_brain_mask)
      }

      if (!is.null(brain_pct)) {
        if (verbose > 1) {
          message("Resampling Brain Percent")
        }
        resampled_brain_pct = resample_to_target(
          prenormalize$brain_pct,
          target = resampled[[1]],
          verbose = verbose > 1,
          interpolator = interpolator,
          copy_origin = copy_origin)
        writenii(resampled_brain_pct, filename = brain_pct_fname)
      }

      # this should be after remasking
      mapply(function(img, fname){
        writenii(img, filename = fname)
      }, resampled, fnames)

      if (!is.null(gold_standard)) {
        if (verbose > 1) {
          message("Resampling Gold Standard")
        }
        # resampled_gs = resample_image(
        #   prenormalize$GOLD_STANDARD,
        #   parameters = c(1, 1, 1),
        #   parameter_type = "mm",
        #   interpolator = "nearestneighbor")
        # using reample_to_target so that you can use
        # genericLabel interpolator
        resampled_gs = resample_to_target(
          prenormalize$GOLD_STANDARD,
          target = resampled[[1]],
          verbose = verbose > 1,
          interpolator = dis_interpolator,
          copy_origin = copy_origin)
        writenii(resampled_gs, filename = gs_fname)
      }

      t1_reg = NULL
      template_fname = fnames[1]

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
      t1 = prenormalize$images$T1
      if (is.null(t1)) {
        stop("T1 is not in prenormalize$images!")
      }
      t1_reg = registration(
        filename = prenormalize$images$T1,
        template.file = template_fname,
        skull_strip = FALSE,
        correct = FALSE,
        outprefix = nii.stub(prenormalize$images$T1),
        verbose = verbose > 1,
        typeofTransform = typeofTransform,
        remove.warp = FALSE,
        interpolator = interpolator,
        copy_origin = copy_origin)
      keeper = function(x) {
        x[ grep("Generic|Warp", x)]
      }
      t1_reg$fwdtransforms = keeper(t1_reg$fwdtransforms)
      t1_reg$invtransforms = keeper(t1_reg$invtransforms)

      t1 = check_ants(prenormalize$images$T1)

      resampled = lapply(
        prenormalize$images,
        function(r) {
          if (copy_origin) {
            r = check_ants(r)
            r = antsCopyOrigin(reference = t1, target = r)
          }
          ants_apply_transforms(
            moving = r,
            fixed = template_fname,
            transformlist = t1_reg$fwdtransforms,
            interpolator = interpolator)
        })

      if (!is.null(brain_mask)) {
        bm = prenormalize$brain_mask
        if (copy_origin) {
          bm = check_ants(bm)
          bm = antsCopyOrigin(reference = t1, target = bm)
        }
        resampled_brain_mask = ants_apply_transforms(
          moving = bm,
          fixed = template_fname,
          transformlist = t1_reg$fwdtransforms,
          interpolator = dis_interpolator)

        writenii(resampled_brain_mask, filename = brain_fname)

        resampled = lapply(
          resampled,
          mask_img,
          mask = resampled_brain_mask)

      }

      if (!is.null(brain_pct)) {
        bm = prenormalize$brain_pct
        if (copy_origin) {
          bm = check_ants(bm)
          bm = antsCopyOrigin(reference = t1, target = bm)
        }
        resampled_brain_pct = ants_apply_transforms(
          moving = bm,
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
        bm = prenormalize$GOLD_STANDARD
        if (copy_origin) {
          bm = check_ants(bm)
          bm = antsCopyOrigin(reference = t1, target = bm)
        }
        resampled_gs = ants_apply_transforms(
          moving = bm,
          fixed = template_fname,
          transformlist = t1_reg$fwdtransforms,
          interpolator = dis_interpolator)
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
      template_fname = fnames[1]

    } else {
      template_fname = switch(
        template,
        Eve = EveTemplate::getEvePath(what = "Brain"),
        MNI = MNITemplate::getMNIPath(what = "Brain", res = "1mm")
      )
      outprefix = nii.stub(prenormalize$images$T1)

      t1_reg = extrantsr::transformlist_from_outprefix(
        outprefix = outprefix,
        typeofTransform = typeofTransform)
      # have_warp = (typeofTransform %in%
      #      c("Translation", "Rigid", "Similarity", "TRSAA")) |
      #     grepl("Rigid", typeofTransform) |
      #     grepl("Affine", typeofTransform)
      # have_warp = !have_warp
      # endings = c("0GenericAffine.mat")
      # if (have_warp) {
      #   endings = c("1Warp.nii.gz", endings)
      # }
      # inv_endings = rev(endings)
      # inv_endings = sub("1Warp", "1InverseWarp", inv_endings)

      # t1_reg = list(
      #   fwdtransforms = paste0(
      #     outprefix,
      #     c("1Warp.nii.gz",
      #       "0GenericAffine.mat")),
      #   invtransforms = paste0(
      #     outprefix,
      #     c("0GenericAffine.mat",
      #       "1InverseWarp.nii.gz")
      #   ))
      t1_reg$typeofTransform = typeofTransform
      t1_reg$interpolator = interpolator
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
    template = template,
    typeofTransform = typeofTransform,
    interpolator = interpolator
  )
  L$spatial_suffix = app
  L$template_fname = checkimg(template_fname)
  # L$GOLD_STANDARD = gold_standard
  L$GOLD_STANDARD = resampled_gs
  L$brain_mask = resampled_brain_mask
  L$reg_to_template = t1_reg
  L$brain_pct = resampled_brain_pct

  return(L)
}