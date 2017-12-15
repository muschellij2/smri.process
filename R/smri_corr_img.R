#' SMRI Local Correlation Images
#'
#' Calculates the local correlation between sets of images
#' @param x List of images
#' @param outdir Output directory
#' @param verbose print diagnostic messages
#' @param mask mask to normalize over
#' @param suffix what the append to end of the output filename
#' @param radius radius of neighborhood.  See \code{\link{corr_img}}
#' @param method Method for correlation. See \code{\link{corr_img}}
#'
#' @return List of filenames
#' @export
#'
#' @importFrom utils combn setTxtProgressBar txtProgressBar
smri_corr_img = function(
  x,
  mask = NULL,
  outdir = tempdir(),
  verbose = TRUE,
  suffix = "_corr",
  radius = rep(1, 3),
  method = c("pearson", "spearman")
){


  nii_names = names(x)
  if (length(nii_names) != length(x)) {
    stop("x must be a named vector or named list")
  }

  if (verbose > 0) {
    message("Correlating Images")
  }
  eg = t(combn(nii_names, m = 2))
  eg = data.frame(eg, stringsAsFactors = FALSE)
  colnames(eg) = c("img1", "img2")
  eg$prefix = paste0(eg$img1, "_", eg$img2)
  eg$filename = paste0(eg$prefix, suffix, ".nii.gz")
  eg$filename = file.path(outdir, eg$filename)


  method = match.arg(method)

  #################################
  # Bias correct the images using N4
  #################################
  fnames = eg$filename
  names(fnames) = eg$prefix

  if (!all_exists(fnames)) {

    mask = checkimg(mask)
    N = nrow(eg)
    if (verbose) {
      pb = txtProgressBar(max = N, initial = 0, style = 3)
    }
    for (irow in seq(N)) {
      ieg = eg[irow,]
      fname = ieg$filename[1]

      img1 = x[[ieg$img1[1]]]
      img2 = x[[ieg$img2[1]]]
      if (verbose > 1) {
        msg = paste0("\nCorrelating the ", ieg$img1[1], " and ", ieg$img2[1])
        message(msg)
      }
      myres = extrantsr::corr_img(
        img1,
        img2,
        mask = mask,
        radius = radius,
        method = method,
        verbose = verbose > 1)
      if (verbose > 1) {
        msg = paste0("Writing ", fname)
        message(msg)
      }
      writenii(myres, filename = fname)
      gc(); gc();
      if (verbose) {
        setTxtProgressBar(pb, value = irow)
      }
    }
    if (verbose) {
      close(pb)
    }

  }
  norm = lapply(
    fnames,
    identity)

  return(norm)
}
