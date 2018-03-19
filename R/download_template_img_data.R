#' @title Download Eve-Template Normed data
#'
#' @description Download Eve-Template Normed data
#' @return Logical indicator if the files were downloaded.
#'
#' @param ... Additionl arguments passed to
#' \code{\link{template_img_data}}
#' @param modality Modality or modalities to download
#' @param statistic Statistic image(s) to download
#' Passed to \code{\link{template_img_data}}
#' @param lib.loc a character vector with path names of R libraries to
#' download the data.
#'
#' @examples
#' tfile = tempfile()
#' dir.create(tfile)
#' eg = template_img_data(
#' modality = "T1_Pre", statistic = "Mean", lib.loc = tfile)
#' dl = download_template_img_data(
#' modality = "T1_Pre", statistic = "Mean", lib.loc = tfile)
#' @export
#' @importFrom utils download.file
download_template_img_data = function(
  ...,
  lib.loc = NULL){

  eg = template_img_data(..., lib.loc = lib.loc)
  img_files = eg$file

  if (!is.null(lib.loc)) {
    desc = system.file("DESCRIPTION", package = "smri.process")
    lib.dir = file.path(lib.loc, "smri.process")
    if (!dir.exists(lib.dir)) {
      dir.create(lib.dir, recursive = TRUE)
    }
    out_desc = file.path(lib.dir, "DESCRIPTION")
    if (!file.exists(out_desc)) {
      if (file.exists(desc)) {
        file.copy(desc, out_desc)
      }
    }
  }

  if (!all(file.exists(img_files))) {
    stubs = eg$fname
    lib_dir = system.file(package = "smri.process",
                          lib.loc = lib.loc)
    for (istub in stubs) {
      urlfile <- file.path(lib_dir, istub)
      url = paste0("http://muschellij2.github.io/smri.process/", istub)
      download.file(url, urlfile, quiet = TRUE)
    }
  }
  eg = template_img_data(..., lib.loc = lib.loc)
  img_files = eg$file

  return(all(file.exists(img_files)))
}

#' @export
#' @rdname download_template_img_data
template_img_data = function(
  modality = c("T1_Pre", "T1_Post", "FLAIR", "T2", "PD"),
  statistic = c("Mean", "Median", "SD"),
  lib.loc = NULL){

  modality = match.arg(modality, several.ok = TRUE)
  statistic = match.arg(statistic, several.ok = TRUE)
  eg = expand.grid(modality = modality,
                   statistic = statistic, stringsAsFactors = FALSE)
  eg$fname = paste0("Template_", eg$modality, "_", eg$statistic, ".nii.gz")
  eg$file = sapply(eg$fname, system.file, package = "smri.process",
                   lib.loc = lib.loc)
  return(eg)
}
