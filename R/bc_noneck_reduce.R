#' @title Bias Correct, Remove Neck, Reduce Image
#'
#' @param x List of images or character vector.  Must be named with
#' the imaging modalities. T1 and FLAIR must be included
#' @param gold_standard Gold Standard image/filename, if applicable
#' @param gs_space space the Gold Standard is located
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#' @return List of output filenames
#' @export
#' @importFrom neurobase applyEmptyImageDimensions
#' @importFrom neurobase writenii
bc_noneck_reduce = function(
  x,
  gold_standard = NULL,
  gs_space = NULL,
  outdir = tempdir(),
  verbose = TRUE){

  rpi_done = rpi_raw(
    x = x,
    gold_standard = gold_standard,
    verbose = verbose)

  gold_standard = rpi_done$GOLD_STANDARD
  rpi_done$GOLD_STANDARD = NULL

  n4 = n4_raw(
    x = rpi_done,
    verbose = verbose,
    outdir = outdir)

  noneck = noneck(
    x = n4,
    verbose = verbose,
    outdir = outdir)

  rm_neck = reduce_img(
    x = noneck,
    verbose = verbose,
    outdir = outdir)

  nn = noneck[[gs_space]]

  #################################
  # Removing dims if gold standard exists
  #################################
  if (!is.null(gold_standard)) {
    les_fname = file.path(
      outdir,
      "GOLD_STANDARD_reduced.nii.gz")
    if (!all_exists(les_fname)) {
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

  rm_neck$GOLD_STANDARD = les_fname

  L = list(
    images = rm_neck,
    gs_space = gs_space,
    outdir = outdir
  )
  return(L)
}




