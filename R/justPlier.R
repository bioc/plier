"justPlier" <-
function(eset=ReadAffy(),replicate=1:length(eset),get.affinities=FALSE,normalize=FALSE,norm.type="together",augmentation=0.1,defaultaffinity=1.0,defaultconcentration=1.0,attenuation=0.005,seaconvergence=0.000001,seaiteration=3000,gmcutoff=0.15,probepenalty=0.001,concpenalty=0.000001,usemm=TRUE,usemodel=FALSE,fitaffinity=TRUE,plierconvergence=0.000001,plieriteration=3000,dropmax=3.0,lambdalimit=0.01,optimization=0) {
  if(normalize) {
    cat("Quantile normalizing...");
    eset <- normalize.AffyBatch.quantiles(eset,norm.type);
    cat("Done.\n");
  }
  else {
    warning("##### normalise=FALSE: No normalization will be performed in plier");
  }
  pns       <- probeNames(eset)
  num_exp   <- length(sampleNames(eset))
  o <- order(pns)
  pns <- pns[o]
  if(missing(replicate)) { replicate <- 1:num_exp }
  else { if(!is.numeric(replicate)) { stop("Parameter 'replicate' must be a vector of integers") } }
  if(length(replicate)!=num_exp) { stop("Parameter 'replicate' should be the same length as the number of samples in 'eset'") }
  pms <- as.double(pm(eset)[o,])
  mms <- mm(eset)

  if(!usemm) { 
    cat("usemm=FALSE: setting mms to 1.0\n")
    mms <- rep(1.0,length(pms))
  }

  else { if(dim(mms)[1] == length(o)) {
      mms <- as.double(mms[o,])
    }	
  
    else {
      cat("Warning: mms only contained dim",dim(mms),"mm probes!\n")
      cat("Creating dummy mms..\n");
      mms <- rep(1.0,length(pms))
    }
  }
  num_probe <- length(pns)

  r <- .C("an_experiment",as.logical(TRUE),as.double(augmentation), as.double(gmcutoff), as.double(probepenalty), as.double(concpenalty), as.double(defaultaffinity), as.double(defaultconcentration), as.double(attenuation), as.double(seaconvergence), as.integer(seaiteration), as.double(plierconvergence), as.integer(plieriteration), as.logical(usemm), as.logical(usemodel), as.logical(fitaffinity), as.double(dropmax), as.double(lambdalimit), as.integer(optimization), as.integer(num_exp), as.integer(num_probe), as.integer(replicate), pms, mms, as.character(pns), concentration=double(num_exp * length(unique(pns))), affinity=double(num_probe), error.code=integer(1),PACKAGE="plier")

  x <- log2(t(matrix(r$concentration,nrow=num_exp)))
  colnames(x) <- sampleNames(eset)
  rownames(x) <- featureNames(eset)
  res <- new("ExpressionSet", 	
             exprs       = x,
             phenoData   = phenoData(eset),
             annotation  = annotation(eset), 
             experimentData = experimentData(eset));
  if(get.affinities) {
    a <-  (r$affinity)
    names(a) <- probeNames(eset)
    res@description@preprocessing$affinity=a;
  }
  return(res)
}

.testoneprobeset <- function() {
# here are the default test parameters

augmentation         <- 0.1
defaultaffinity      <- 1.0
defaultconcentration <- 1.0
attenuation          <- 0.005
seaconvergence       <- 0.000001
seaiteration         <- 2000
gmcutoff             <- 0.15
probepenalty         <- 0.001
concpenalty          <- 0.000001
usemm                <- TRUE
usemodel             <- FALSE
fitaffinity          <- FALSE
plierconvergence     <- 0.000001
plieriteration       <- 3000
dropmax              <- 3.0
lambdalimit          <- 0.01
optimization         <- 0


pm <- c(4071.0,3742.0,3517.0,4231.0,4037.0,3615.0,6374.0,5431.0,5102.0)
mm <- c(503.0,377.0,321.0,354.0,353.0,362.0,1693.0,1436.0,1640.0)
num_exp   <- 3
num_probe <- 3
replicate <- 1:num_exp
 r <- .C("one_probeset",as.logical(TRUE),as.double(augmentation), as.double(gmcutoff), as.double(probepenalty), as.double(concpenalty), as.double(defaultaffinity), as.double(defaultconcentration), as.double(attenuation), as.double(seaconvergence), as.integer(seaiteration), as.double(plierconvergence), as.integer(plieriteration), as.logical(usemm), as.logical(usemodel), as.logical(fitaffinity), as.double(dropmax), as.double(lambdalimit), as.integer(optimization), as.integer(num_exp), as.integer(num_probe), as.integer(replicate), as.double(pm), as.double(mm), concentration=double(num_exp), affinity=double(num_probe), error.code=integer(1),PACKAGE="plier")
print(r)
}
