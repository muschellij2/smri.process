
#' Multiple FAST for PD/T1/T2
#'
#' @param x List of processed images with T1/T2/PD
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#' @param bias_correct should bias correction be done?  passed to
#'  \code{\link{fast}}
#' @param ... additional arguments to \code{\link{fast}}
#'
#' @return List of filenames
#' @export
#' @importFrom fslr fast
multi_fast = function(
  x,
  outdir = tempdir(),
  bias_correct = FALSE,
  verbose = TRUE, ...) {

  nii_names = names(x)
  nii_names = toupper(nii_names)
  if (length(nii_names) != length(x)) {
    stop("x must be a named vector or named list")
  }

  mod = "T1"
  mods = c("T1", "T2", "PD") #, "FLAIR")
  out = vector(mode = "list", length = 3)
  names(out) = mods

  for (mod in mods) {
    if (mod %in% nii_names) {
      img = x[[mod]]
      if (verbose > 0) {
        message(
          paste0(
            "Running FAST on ",
            mod))
      }
      fast_outfile = file.path(
        outdir,
        paste0(mod,
               "_FAST"))
      out_type = c("seg", "mixeltype", "pve_0",
                   "pve_1", "pve_2", "pveseg")
      fnames = paste0(
        fast_outfile,
        "_",
        out_type, ".nii.gz")
      if (!all_exists(fnames)) {
        image_type = mod
        if (mod == "FLAIR") {
          image_type = "T2" # took out FLAIR because unsure of htis
        }
        res = fast(
          img,
          bias_correct = FALSE,
          outfile = fast_outfile,
          verbose = verbose > 1,
          retimg = FALSE,
          type = image_type,
          ...)
      } else {
        res = fnames
      }
      out[[mod]] = res
    } else {
      out[[mod]] = NULL
    }
  }
  nulls = sapply(out, is.null)
  out = out[!nulls]
  if (length(out) == 0) {
    out = NULL
  }
  return(out)
}