#' Register Images to T1 image
#'
#' @param x List of images or character vector.  Must be named with
#' the imaging modalities. T1 and FLAIR must be included
#' @param gs_space space the Gold Standard is located (if applicable)
#' @param interpolator interpolation done in \code{\link{registration}}
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#' @param gs_interpolator interpolation done in \code{\link{registration}}
#' for gold standard
#' @param suffix Name to append to the image filename
#'
#' @return List of Images
#' @export
#' @importFrom neurobase nii.stub
#' @importFrom extrantsr ants_apply_transforms
reg_to_t1 = function(
  x,
  gs_space = NULL,
  interpolator = "Linear",
  gs_interpolator = "NearestNeighbor",
  outdir = tempdir(),
  verbose = TRUE,
  suffix = "_regtoT1"
) {

  nii_names = names(x)
  nii_names = toupper(nii_names)
  if (length(nii_names) != length(x)) {
    stop("x must be a named vector or named list")
  }
  if (!("T1" %in% nii_names)) {
    stop("T1 must be in the images vector")
  }
  x = checkimg(x)
  x = lapply(x, identity)
  names(x) = nii_names

  gold_standard = x$GOLD_STANDARD
  nii_names = nii_names[ !(nii_names %in% "GOLD_STANDARD") ]
  x = x[nii_names]
  if (!is.null(gold_standard)) {
    if (is.null(gs_space)) {
      stop("gs_space must be specified if gold_standard specified")
    }
    gs_space = match.arg(gs_space, choices = names(x))
    gold_standard = checkimg(gold_standard)
  }

  ####################################################
  # Making fname stubs
  ####################################################
  fnames = file.path(
    outdir,
    paste0(nii_names,
           suffix,
           ".nii.gz"))
  names(fnames) = nii_names

  les_fname = file.path(
    outdir,
    paste0("GOLD_STANDARD", suffix, ".nii.gz")
  )
  if (!is.null(gold_standard)) {
    all_fnames = c(fnames, GOLD_STANDARD = les_fname)
  } else {
    all_fnames = fnames
  }

  if (verbose > 0) {
    message("Registering Images to T1")
  }

  if (!all_exists(all_fnames)) {
    # t1 doesn't need to register to itself
    t1 = x$T1
    t1_fname = fnames["T1"]
    file.copy(t1, t1_fname, overwrite = TRUE)

    # removing T1 from the image list
    x = x[ !(names(x) %in% "T1")]
    fnames = fnames[ !(names(fnames) %in% "T1")]

    # run registration
    reg = mapply(function(run_img, run_name, fname) {
      if (verbose > 0) {
        message(paste0("Registering ", run_name, " to T1"))
      }
      outprefix = nii.stub(fname)
      res = registration(
        filename = run_img,
        template.file = t1,
        skull_strip = FALSE,
        correct = FALSE,
        outfile = fname,
        outprefix = outprefix,
        typeofTransform = "Rigid",
        interpolator = interpolator,
        verbose = verbose > 1)
    }, x, names(x), fnames, SIMPLIFY = FALSE)

    # gold standard applying
    if (!is.null(gold_standard)) {
      # if T1 space - then just copy over
      # no registration needed
      if (gs_space != "T1") {
        transformlist = reg[[gs_space]]$fwdtransforms
        gs_img = ants_apply_transforms(
          fixed = t1,
          moving = gold_standard,
          transformlist = transformlist,
          interpolator = gs_interpolator)
        writenii(gs_img, filename = les_fname)
      } else {
        file.copy(gold_standard, les_fname, overwrite = TRUE)
      }
    }
  }

  res_imgs = lapply(
    all_fnames,
    identity)
  return(res_imgs)

}