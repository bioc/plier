\name{justPlier}
\alias{justPlier}

\title{ Implements the PLIER algorithm }
\description{
 Provides a wrapper around Affymetrix's API to provide an implementation of the PLIER alogrimth
}
\usage{
justPlier(eset=ReadAffy(),replicate=1:length(eset),get.affinities=FALSE,normalize=F,norm.type="together",augmentation=0.1,defaultaffinity=1.0,defaultconcentration=1.0,attenuation=0.005,seaconvergence=0.000001,seaiteration=3000,gmcutoff=0.15,probepenalty=0.001,concpenalty=0.000001,usemm=TRUE,usemodel=FALSE,fitaffinity=T,plierconvergence=0.000001,plieriteration=3000,dropmax=3.0,lambdalimit=0.01,optimization=0)
}
\arguments{
  \item{eset}{ An AffyBatch object containing the raw data }
  \item{replicate}{ A factor containing the replicate structure to use for grouping samples }
  \item{get.affinities}{ If TRUE, then return affinities in the preprocessing slot of the exprSet object }
  \item{normalize}{ If TRUE then apply quantile normalization to the probes before generating expression calls }
  \item{norm.type}{ Can be 'separate', 'pmonly', 'mmonly' or 'together' }
  \item{augmentation}{ Model parameter }
  \item{defaultaffinity}{ Model parameter }
  \item{defaultconcentration}{ Model parameter }
  \item{attenuation}{ Model parameter }
  \item{seaconvergence}{ Model parameter }
  \item{seaiteration}{ Model parameter }
  \item{gmcutoff}{ Model parameter }
  \item{probepenalty}{ Model parameter }
  \item{concpenalty}{ Model parameter }
  \item{usemm}{ Model parameter }
  \item{usemodel}{ Model parameter }
  \item{fitaffinity}{ Model parameter }
  \item{plierconvergence}{ Model parameter }
  \item{plieriteration}{ Model parameter }
  \item{dropmax}{ Model parameter }
  \item{lambdalimit}{ Model parameter }
  \item{optimization}{ Model parameter }
}
\details{
}

\value{
  An Expression set containing PLIER generated expression calls 
}
\references{ \url{bioinf.picr.man.ac.uk} \url{www.affymetrix.com} }
\author{ Crispin J Miller (wrapper), Earl Hubbell (algorithm) }


\seealso{ normalize.AffyBatch.quantiles }
\examples{

}
\keyword{ misc }

