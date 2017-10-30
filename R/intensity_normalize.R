
#' Normalize using Min/Max
#' @param img Image to normalize
#' @param mask Mask same dimensions as image
#' @param remask Should the image be remasked after normalizing?
#'
#' @return Object of same class as img
#' @export
#' @importFrom neurobase finite_img mask_vals
minmax = function(img, mask = NULL, remask = FALSE) {
  img = check_nifti(img, allow.array = TRUE)

  dimg = dim(img)
  if (is.null(mask)) {
    mask = array(1, dim = dimg)
  }
  vals = mask_vals(img, mask = mask)
  r = range(vals)
  img = (img - r[1]) / (r[2] - r[1])
  img = finite_img(img)
  if (remask) {
    img = mask_img(img, mask)
  }
  return(img)
}

#' SMRI Normalize
#'
#' @param x List of images
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#' @param mask mask to normalize over
#' @param normalization Normalization Method
#' @param remask should mask_img be applied to the output images
#' after normalization
#' @param suffix what the append to end of the output filename
#' @param ... Additional arguments to MALF
#'
#' @return List of images
#' @export
#'
#' @importFrom neurobase zscore_img quantile_img
#' @importFrom WhiteStripe whitestripe_norm whitestripe whitestripe_hybrid
intensity_normalize = function(
  x,
  mask = NULL,
  normalization = c("z", "trimmed_z", "quantile", "minmax",
                    "whitestripe_T1", "whitestripe_T2",
                    "whitestripe_hybrid"),
  outdir = tempdir(),
  verbose = TRUE,
  remask = FALSE,
  suffix = "_z",
  ...
) {


  nii_names = names(x)
  if (length(nii_names) != length(x)) {
    stop("x must be a named vector or named list")
  }

  if (verbose > 0) {
    message("Intensity Normalization")
  }
  #################################
  # Bias correct the images using N4
  #################################
  fnames = file.path(
    outdir,
    paste0(nii_names,
           suffix,
           ".nii.gz"))
  names(fnames) = nii_names

  if (!all_exists(fnames)) {

    args = list(...)
    nargs = names(args)

    normalization = match.arg(normalization)
    if (normalization %in% c("z", "trimmed_z")) {
      if (normalization == "z") {
        args$trim = 0
      }
      args$mask = mask
      args$centrality = "trimmed_mean"
      args$variability = "trimmed_sd"
      norm = lapply(x, function(r) {
        args$img = r
        do.call("zscore_img", args = args)
      })
    }

    if (normalization %in% "quantile") {
      args$mask = mask
      norm = lapply(x, function(r) {
        args$img = r
        do.call("quantile_img", args = args)
      })
    }

    if (normalization %in% c("minmax")) {
      args$mask = mask
      norm = lapply(x, function(r) {
        args$img = r
        do.call("minmax", args = args)
      })
    }

    if (normalization %in%
        c("whitestripe_T1",
          "whitestripe_T2",
          "whitestripe_hybrid")) {
      x = check_nifti(x)
      stripped = TRUE

      if (normalization %in%
          c("whitestripe_T1",
            "whitestripe_T2")) {

        type = sub("whitestripe_", "", normalization)
        type = toupper(type)
        stopifnot(type %in% nii_names)

        img = x[[type]]
        if (!is.null(mask)) {
          # joint mask
          mask = check_nifti(mask)
          run_mask = mask > 0 & img > 0

          img = mask_vals(img, mask = run_mask)
          class(img) = "img_voi"
          stripped = FALSE
        }

        args$img = img
        args$type = type
        func = "whitestripe"
      }
      if (normalization %in% "whitestripe_hybrid") {
        stopifnot(all(c("T1", "T2") %in% nii_names))
        t1 = x$T1
        t2 = x$T2
        if (!is.null(mask)) {
          # joint mask
          mask = check_nifti(mask)
          run_mask = t1 > 0 & t2 > 0 & mask > 0

          t1 = mask_vals(t1, mask = run_mask)
          class(t1) = "img_voi"

          t2 = mask_vals(t2, mask = run_mask)
          class(t2) = "img_voi"

          stripped = FALSE
        }
        args$t1 = t1
        args$t2 = t2
        func = "whitestripe_hybrid"
      }

      # data is usually stripped
      if (!"stripped" %in% nargs) {
        args$stripped = stripped
      }
      args$verbose = verbose
      ws = do.call(func, args = args)
      indices = ws$whitestripe.ind
      ws_mask = ws$mask.img
      if (!is.null(mask)) {
        ind = which(run_mask > 0)
        indices = ind[indices]
        ws_mask = run_mask * 0
        ws_mask[indices] = 1
      }
      fname = file.path(
        outdir,
        paste0(normalization,
               "_mask.nii.gz"))
      writenii(mask, fname)

      norm = lapply(
        x,
        WhiteStripe::whitestripe_norm,
        indices = indices)

    }

    if (remask) {
      if (!is.null(mask)) {
        mask = check_nifti(mask)
        norm = lapply(norm, mask_img, mask = mask)
      }
    }

    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, norm, fnames)
  }

  norm = lapply(
    fnames,
    identity)

  return(norm)
}
