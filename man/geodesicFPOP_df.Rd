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
geodesicFPOP_vec(angle.vec, penalty, db.file=NULL, container.str="list")
geodesicFPOP_dir(problem.dir, penalty, db.file=NULL, container.str="list")
geodesicFPOP_df(
  radians.df, penalty, base.dir = tempdir(), db.file=NULL, container.str="list")
}
\arguments{
  \item{angle.vec}{numeric vector of data to segment.}
  \item{radians.df}{data frame of data to segment, with columns start, end, radians.}
  \item{penalty}{numeric penalty for each change-point.}
  \item{base.dir}{where to save output files.}
  \item{problem.dir}{directory with data to segment.}
  \item{db.file}{where to write temporary cost function database, can be
    useful to put this on a drive with faster I/O.}
  \item{container.str}{string: container to store cost functions, either
    \code{"list"} or \code{"map"}}
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
fit <- geodesichange::geodesicFPOP_vec(angles, penalty=Inf)
fit$segments$param
angle.vec <- c(0.1, 6, 6.1, 3, 3.1, 2.9)
result <- geodesichange::geodesicFPOP_vec(angle.vec, penalty=1)
result$segments
result$loss

}
