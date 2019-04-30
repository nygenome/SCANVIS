\name{gene2roi}
\alias{gene2roi}
\title{gene name/s to region of interest}
\description{Converts gene name/s to genomic coordinates using gene annotation file from SCANVIS.annotation}
\usage{
gene2roi(g,gen)
}
\arguments{
  \item{g}{vector of one or more gene names or gene ids in the same chromosome}
  \item{gen}{gene annotation object as output by SCANIVS.annotation}
}
\details{This function is called upon by SCANVIS.linkvar and SCANVIS.visual}

\value{chr, start and end of the union of genomic intervals that overlap the genes in g}

\examples{
data(SCANVIS_examples)
g=c('TDRD6','PLA2G7')
roi=gene2roi(g,gen19)
}
\keyword{interval}
\keyword{gene}
\keyword{SCANVIS.annotation}
\keyword{SCANVIS.linkvar}
\keyword{SCANVIS.visual}

