#' Threshold Quantile Image
#'
#' @param img Character vector, or object of class \code{nifti}
#' @param mask Mask to determine cumulative distribution function (cdf) from
#' @param probs Probability to threshold at
#'
#' @return Object of class \code{\link{nifti}}
#' @export
quantile_mask = function(img, mask, probs = 0.85) {
  qimg = quantile_img(img = img, mask = mask)
  mask = qimg >= probs
}