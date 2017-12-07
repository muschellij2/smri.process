#' Bias Correct, Remove Neck, Reduce Image
#'
#' @param x List of images or character vector.  Must be named with
#' the imaging modalities. T1 and FLAIR must be included
#' @param gold_standard Gold Standard image/filename, if applicable
#' @param gs_space space the Gold Standard is located
#' @param probs passed to \code{\link{winsor}} for Winsorization
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#' @return List of output filenames
#' @export
#' @importFrom neurobase applyEmptyImageDimensions
#' @importFrom neurobase writenii robust_window
bc_noneck_reduce = function(
  x,
  gold_standard = NULL,
  gs_space = NULL,
  probs = c(0, 0.999),
  outdir = tempdir(),
  verbose = TRUE){

  rpi_done = rpi_raw(
    x = x,
    gold_standard = gold_standard,
    verbose = verbose)

  gold_standard = rpi_done$GOLD_STANDARD
  rpi_done$GOLD_STANDARD = NULL

  suffix = "_N4"
  n4 = n4_raw(
    x = rpi_done,
    verbose = verbose,
    mask = NULL,
    outdir = outdir,
    suffix = suffix
    )

  suffix = paste0(suffix, "_noneck")
  noneck = noneck(
    x = n4,
    verbose = verbose,
    outdir = outdir,
    suffix = suffix
  )

  suffix = paste0(suffix, "_reduced")
  rm_neck = reduce_img(
    x = noneck,
    verbose = verbose,
    outdir = outdir,
    suffix = suffix
  )


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
      dd = mask_reduce(nn)
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




