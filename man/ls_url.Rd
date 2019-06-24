\name{ls_url}
\alias{ls_url}
\title{list files available at annotation/gencode url}
\description{Function called upon by SCANVISannotation for GTF file pulldown from url supplied}
\usage{
ls_url(url)
}
\arguments{
  \item{url}{url to GTF files}
}
\details{calls upon functions in Rcurl and rtracklayer}

\value{a list of files for download at url}

\examples{
ftpfiles=ls_url('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/')
}
\keyword{gencode}
\keyword{annotation}
\keyword{GTF}
