#' Skull Strip image
#'
#' @param file Filename of T1 image
#' @param template Template full image
#' @param template_mask Template brain mask
#' @param verbose print diagnostic messages
#' @param n_iter Number of iterations undergone with \code{\link{abpN4}}
#' and \code{\link{abpBrainExtraction}}
#' @param pad should zero padding be done
#'
#' @return A list of output images, brain and corrected
#' @export
#' @importFrom ANTsR abpBrainExtraction abpN4
#' @importFrom ANTsRCore is.antsImage antsImageClone
#' @importFrom methods formalArgs
#' @importFrom stats median
#' @importFrom penn115 penn115_image_fname penn115_brain_fname
#' @importFrom penn115 penn115_brain_mask_fname
#' @importFrom extrantsr tempants oro2ants
#' @importFrom neurobase zero_pad check_nifti
n4_skull_strip = function(
  file,
  template = penn115_image_fname(),
  template_mask = penn115_brain_mask_fname(),
  verbose = TRUE,
  pad = TRUE,
  n_iter = 2) {

  # load the file
  if (verbose) {
    message('Loading file:')
  }
  reader = function(x) {
    if (!is.antsImage(x)) {
      y = antsImageRead(x)
    } else {
      y = antsImageClone(x)
    }
    return(y)
  }
  if (pad) {
    img = check_nifti(file)
    img = zero_pad(img, kdim = c(3, 3, 3))
    file = checkimg(img)
    rm(img)
  }

  subimg = reader(file)
  submask = subimg * 0 + 1
  # load template files

  if (verbose) {
    message("Loading template...")
  }


  temp = reader(template)
  tempmask =  reader(template_mask)

  # two rounds of N4-BrainExtract to skull strip
  if (verbose) {
    message("Skull stripping...")
  }
  for (i in 1:n_iter) {
    if (verbose) {
      message(paste0("Running iteration ", i))
    }
    args = list(
      img = subimg,
      mask = submask,
      verbose = verbose > 1)
    if (!"..." %in% formalArgs(abpN4)) {
      args$verbose = NULL
    }
    subimg = do.call(abpN4, args = args)

    args = list(
      img = subimg,
      tem = temp,
      temmask = tempmask,
      verbose = verbose > 1)
    if (!"verbose" %in% formalArgs(abpBrainExtraction)) {
      args$verbose = NULL
    }
    bextract = do.call(abpBrainExtraction, args = args)
    rm(submask)
    submask = bextract$bmask
    if (verbose) {
      message(
        paste0("SubMask number of voxels ", sum(submask))
      )
    }
  }
  if (pad) {
    img = check_nifti(submask)
    img = zero_pad(img, kdim = c(3, 3, 3), invert = TRUE)
    submask = oro2ants(img)
  }

  # simg = n4 * submask

  # L = list(n4 = n4,
  #          brain_mask = submask,
  #          n4_brain = simg)
  # return(L)
  return(submask)
}
