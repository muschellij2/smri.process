#' @title All Files exist
#' @description Wrapper for all files existing
#' @param ... files to be passed to \code{\link{file.exists}}
#' @return Logical length of \code{...}
#' @export
all_exists = function(...){
  all(file.exists(...))
}