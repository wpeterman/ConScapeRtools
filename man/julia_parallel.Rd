% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/julia_parallel.R
\name{julia_parallel}
\alias{julia_parallel}
\title{Setup Julia for parallel processing}
\usage{
julia_parallel(workers, jl_home, attempts = 9)
}
\arguments{
\item{workers}{Number of parallel workers to create}

\item{attempts}{Number of attempts to try and process files}

\item{julia_home}{Path to the Julia 'bin' directory}
}
\value{
NULL.
}
\description{
Set up parallel processing to run ConScape with parallel.
}
\author{
Bill Peterman
}
\keyword{internal}
