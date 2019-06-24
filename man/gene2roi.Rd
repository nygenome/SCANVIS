\name{gene2roi}
\alias{gene2roi}
\title{gene name/s to region of interest}
\description{Converts gene name/s to genomic coordinates using gene annotation file from SCANVISannotation}
\usage{
gene2roi(g,gen)
}
\arguments{
  \item{g}{vector of one or more gene names or gene ids in the same chromosome}
  \item{gen}{gene annotation object as output by SCANIVS.annotation}
}
\details{This function is called upon by SCANVISlinkvar and SCANVISvisual}

\value{chr, start and end of the union of genomic intervals that overlap the genes in g}

\examples{
data(SCANVISexamples)
g=c('TDRD6','PLA2G7')
roi=gene2roi(g,gen19)
}
\keyword{interval}
\keyword{gene}
\keyword{SCANVISannotation}
\keyword{SCANVISlinkvar}
\keyword{SCANVISvisual}

