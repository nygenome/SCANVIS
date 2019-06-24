\name{IR2Mat}
\alias{IR2Mat}
\title{IRanges to Matrix}
\description{converts IRanges interval object to matrix}
\usage{IR2Mat(I)}
\arguments{
  \item{I}{IRanges interval object}
}
\details{This function is called upon by SCANVISscan}
\value{a matrix with start and end coordinates for the intervals in I}
\examples{
library(IRanges)
I=IRanges(1:10,21:30)
m=IR2Mat(I)
}
\keyword{SCANVISscan}
