\name{SCANVIS.annotation}
\alias{SCANVIS.annotation}
\title{assembles annotation from gtf file into SCANVIS-readable format}
\description{This function ftps to the supplied gtf url, downloads gtf
to current directory and assembles annotation details into a SCANVIS-readable object}
\usage{
SCANVIS.annotation(ftp.url)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ftp.url}{web url to gtf file}
}
%\details{Requires Rcurl and rtracklayer}

\value{a gencode object compatible (and required) for use with most SCANVIS functions}

%\references{
%% ~put references to the literature/web site here ~
%}
%\author{Phaedra Agius <pagius@nygenome.org>}

\note{Web access required. If variants are available and intended for use with SCANVIS.linkvar, the gencode reference genome must be the same as that used for the variant calls.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
\dontrun{gen19=SCANVIS.annotation('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/')}
}
\keyword{gencode}
\keyword{annotation}
