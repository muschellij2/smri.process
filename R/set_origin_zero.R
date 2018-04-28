#' Set Origin to Zero for the Images
#'
#' @param fnames character vector of filenames
#'
#' @return Logical indicating if the setting was successful
#'
#' @export
#' @importFrom ANTsRCore antsImageRead antsSetOrigin antsGetOrigin
#' @importFrom ANTsRCore antsImageWrite
set_origin_zero = function(fnames) {
  stopifnot(is.character(fnames))
  res = sapply(fnames, function(fname) {
    img = antsImageRead(fname)
    orig_origin = antsGetOrigin(img)
    if (all(orig_origin == 0)) {
      return(TRUE)
    }
    antsSetOrigin(img, origin = c(0, 0, 0))
    antsImageWrite(img, filename = fname)
    img_again = antsImageRead(fname)
    all(antsGetOrigin(img_again) == c(0, 0, 0)   )
  })
  return(res)
}