% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_template_img_data.R
\name{download_template_img_data}
\alias{download_template_img_data}
\alias{template_img_data}
\title{Download Eve-Template Normed data}
\usage{
download_template_img_data(..., lib.loc = NULL)

template_img_data(modality = c("T1_Pre", "T1_Post", "FLAIR", "T2", "PD"),
  statistic = c("Mean", "Median", "SD"), lib.loc = NULL)
}
\arguments{
\item{...}{Additionl arguments passed to
\code{\link{template_img_data}}}

\item{lib.loc}{a character vector with path names of R libraries to
download the data.}

\item{modality}{Modality or modalities to download}

\item{statistic}{Statistic image(s) to download
Passed to \code{\link{template_img_data}}}
}
\value{
Logical indicator if the files were downloaded.
}
\description{
Download Eve-Template Normed data
}
\examples{
tfile = tempfile()
dir.create(tfile)
eg = template_img_data(
modality = "T1_Pre", statistic = "Mean", lib.loc = tfile)
dl = download_template_img_data(
modality = "T1_Pre", statistic = "Mean", lib.loc = tfile)
}
