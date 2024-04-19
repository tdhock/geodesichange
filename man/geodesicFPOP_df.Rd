\name{geodesicFPOP_df}
\alias{geodesicFPOP_df}
\alias{geodesicFPOP_vec}
\alias{geodesicFPOP_dir}
\title{
geodesic Functional Pruning Optimal Partitioning
}
\description{
Optimal change-point detection for 1d angular data.
}
\usage{
geodesicFPOP_vec <- function(angle.vec, pen.num)
geodesicFPOP_dir(problem.dir, penalty.param, db.file=NULL)
geodesicFPOP_df(radians.df, pen.num, base.dir = tempdir())
}
\arguments{
  \item{angle.vec}{numeric vector of data to segment.}
  \item{radians.df}{data frame of data to segment, with columns start, end, radians.}
  \item{pen.num}{numeric penalty for each change-point.}
  \item{pen.param}{numeric penalty for each change-point.}
  \item{base.dir}{where to save output files.}
  \item{problem.dir}{directory with data to segment.}
  \item{db.file}{where to write temporary cost function database, can be
    useful to put this on a drive with faster I/O.}
}
\value{
  List of segmentation results.
}
\author{
  Toby Dylan Hocking
}
\examples{

angles <- c(0.1,1,6)
median(angles)
fit <- geodesichange::geodesicFPOP_vec(angles, pen.num=Inf)
fit$segments$param
angle.vec <- c(0.1, 6, 6.1, 3, 3.1, 2.9)
result <- geodesichange::geodesicFPOP_vec(angle.vec, pen.num=1)
result$segments
result$loss

}
