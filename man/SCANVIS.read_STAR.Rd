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

\examples{
tmp=cbind(rep('chr1',10),seq(100,1000,100),seq(100,1000,100)+99,rep(2,10),rep(2,10),rep(0,10),c(rep(500,5),rep(8,5)),rep(0,10),rep(50,10))
write.table(tmp,'tmp',sep='\t',quote=F,row.names=F,col.names=F)
sj=SCANVIS.read_STAR('tmp')
#sj is now suitable as input for SCANVIS.scan
}
\keyword{STAR}
