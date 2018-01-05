#' MALF Tissue and Structure Segmentation
#'
#' @param t1 T1-weighted image/filename
#' @param outdir Output directory
#' @param num_templates Number of templates to use for MALF
#' @param interpolator interpolation done in \code{\link{registration}}
#' @param typeofTransform Transformation to align the templates to the T1
#' @param verbose Print diagnostic messages
#' @param ... additional options to pass to \code{\link{malf}}
#'
#' @return Object of class nifti
#' @export
#' @importFrom extrantsr malf
#' @importFrom neurobase niftiarr
#' @importFrom neurobase check_outfile
#' @importFrom extrantsr reapply_malf
t1_segment = function(
  t1,
  outdir = tempdir(),
  num_templates = 35,
  interpolator = "genericLabel",
  typeofTransform = "SyNAggro",
  verbose = TRUE,
  ...){

  ind = seq(num_templates)
  templates = malf.templates::malf_images()
  brains = templates$brains[ind]
  tissues = templates$tissues[ind]
  labels = templates$labels[ind]

  #######################################
  # Try MALF for Tissues with MASS Templates
  #######################################
  tissue_fname = file.path(
    outdir,
    "TISSUES.nii.gz")
  tissue_list_fnames = file.path(
    outdir,
    paste0(c("CSF", "GM", "WM"), ".nii.gz")
  )
  label_fname = file.path(
    outdir,
    "STRUCTURES.nii.gz")

  fnames = c(TISSUES = tissue_fname,
             STRUCTURES = label_fname,
             CSF = tissue_list_fnames[1],
             GM = tissue_list_fnames[2],
             WM = tissue_list_fnames[3])

  if (!all_exists(fnames)) {

    args = list(...)
    inverted = args$inverted

    regs = malf(
      infile = t1,
      template.images = brains,
      template.structs = tissues,
      keep_images = TRUE,
      retimg = TRUE,
      func = "pct",
      keep_regs = TRUE,
      interpolator = "Linear",
      other_interpolator = interpolator,
      invert_interpolator = interpolator,
      verbose = verbose,
      typeofTransform = typeofTransform,
      ...
    )

    tissue_img = seg_list_to_hard(
      regs$outimg,
      ties.method = "first") - 1
    writenii(tissue_img, filename = tissue_fname)

    tissue_list = regs$outimg[2:4]

    mapply(function(img, filename) {
      writenii(img, filename = filename)
    }, tissue_list, tissue_list_fnames)

    if (is.null(inverted)) {
      inverted = FALSE
    }
    label_img = reapply_malf(
      infile = t1,
      regs = regs$regs,
      template.structs = labels,
      keep_images = TRUE,
      retimg = FALSE,
      outfile = label_fname,
      func = "mode",
      interpolator = interpolator,
      verbose = verbose,
      inverted = inverted
    )
    label_img = NULL
  }
  outfiles = lapply(fnames,
                    identity)
  return(outfiles)
}


seg_list_to_hard = function(
  img,
  ties.method = c("first", "last", "random")
){
  stopifnot(inherits(img, "list"))
  xmax = sapply(img, c)
  ties.method = match.arg(
    ties.method,
    c("first", "last", "random")
  )

  maxs = max.col(xmax,
                 ties.method = ties.method)
  res = neurobase::niftiarr(img[[1]], maxs)
  res
}