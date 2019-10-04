#' Bias Correct, Remove Neck, Reduce Image
#'
#' @param x List of images or character vector.  Must be named with
#' the imaging modalities.
#' @param remove_negative Force values less than zero to zero
#' @param gold_standard Gold Standard image/filename, if applicable
#' @param gs_space space the Gold Standard is located
#' @param probs passed to \code{\link{winsor}} for Winsorization
#' @param non_zero Should zeroes be excluded from the calculation
#' of quantiles? Passed to \code{\link{winsor}}.
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#' @param reduce If \code{TRUE}, the image size will be reduced using
#' \code{\link{reduce_img}}.
#' @param cleanup Argument to \code{\link{oMask}}.  If >0, morphological
#' operations will be applied to clean up the mask by
#' eroding away small or weakly-connected areas, and closing holes.
#' @return List of output filenames
#' @export
#' @importFrom neurobase applyEmptyImageDimensions
#' @importFrom neurobase writenii robust_window
bc_noneck_reduce = function(
  x,
  remove_negative = TRUE,
  gold_standard = NULL,
  gs_space = NULL,
  reduce = TRUE,
  cleanup = 1,
  probs = c(0, 0.999),
  non_zero = TRUE,
  outdir = tempdir(),
  verbose = TRUE){

  rpi_done = rpi_raw(
    x = x,
    gold_standard = gold_standard,
    verbose = verbose)

  gold_standard = rpi_done$GOLD_STANDARD
  rpi_done$GOLD_STANDARD = NULL

  rm_negative = function(x) {
    tmp_res = lapply(x, function(r) {
      img = RNifti::readNifti(r)
      img[is.na(img)] = 0
      if (min(img) < 0) {
        img[img < 0] = 0
        RNifti::writeNifti(img, r)
        rm(img); gc()
      }
      return(TRUE)
    })
    rm(tmp_res);
    return(x)
  }
  if (remove_negative) {
    tmp_res = rm_negative(rpi_done)
  }
  suffix = "_N4"
  n4 = n4_raw(
    x = rpi_done,
    verbose = verbose,
    mask = NULL,
    outdir = outdir,
    suffix = suffix
  )

  if (remove_negative) {
    tmp_res = rm_negative(n4)
  }

  suffix = paste0(suffix, "_noneck")
  noneck = noneck(
    x = n4,
    verbose = verbose,
    outdir = outdir,
    suffix = suffix
  )

  if (reduce) {
    suffix = paste0(suffix, "_reduced")
    rm_neck = reduce_img(
      x = noneck,
      verbose = verbose,
      outdir = outdir,
      suffix = suffix,
      cleanup = cleanup
    )
  } else {
    rm_neck = noneck
  }


  #################################
  # Removing dims if gold standard exists
  #################################
  if (!is.null(gold_standard)) {
    les_fname = file.path(
      outdir,
      paste0("GOLD_STANDARD", suffix, ".nii.gz")
    )
    if (!all_exists(les_fname)) {
      nn = noneck[[gs_space]]
      dd = mask_reduce(nn,
                       cleanup = cleanup)
      les_rm_neck = applyEmptyImageDimensions(
        img = gold_standard,
        inds = dd$inds
      )
      writenii(les_rm_neck,
               filename = les_fname)
    }
  } else {
    les_fname = NULL
  }

  suffix = paste0(suffix, "_winsor")
  win_imgs = winsor(
    x = rm_neck,
    verbose = verbose,
    outdir = outdir,
    probs = probs,
    non_zero = non_zero,
    suffix = suffix)

  win_imgs$GOLD_STANDARD = les_fname

  L = list(
    images = win_imgs,
    gs_space = gs_space,
    outdir = outdir,
    suffix = suffix
  )
  return(L)
}




