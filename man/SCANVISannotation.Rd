\name{SCANVISannotation}
\alias{SCANVISannotation}
\title{assembles annotation from gtf file into SCANVISreadable format}
\description{This function ftps to the supplied gtf url, downloads gtf
to current directory and assembles annotation details into a SCANVISreadable object}
\usage{
SCANVISannotation(ftp.url)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ftp.url}{
}
}
%\details{Requires Rcurl and rtracklayer}

\value{a gencode object compatible (and required) for use with most SCANVISfunctions}

%\references{
%% ~put references to the literature/web site here ~
%}
%\author{Phaedra Agius <pagius@nygenome.org>}

\note{Web access required. If variants are available and intended for use with SCANVISlinkvar, the gencode reference genome must be the same as that used for the variant calls.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
\dontrun{gen19=SCANVISannotation('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/')}
}
\keyword{gencode}
\keyword{annotation}
