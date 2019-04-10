\name{SCANVIS.read_STAR}
\alias{SCANVIS.read_STAR}
\title{upload SJ.tab STAR file in SCANVIS use}
\description{This function is a little wrapper for reading in splice junction details from the SJ.tab file output by the STAR alignment tool.}
\usage{
SCANVIS.read_STAR(starSJ)
}
\arguments{
  \item{starSJ}{url to SJ file output by STAR aligner}
}

\value{SJ data in matrix format as required for SCANVIS functions}

%\author{Phaedra Agius <pagius@nygenome.org>}

%\examples{
%sj=SCANVIS.read_STAR(??)
%}
\keyword{STAR}
