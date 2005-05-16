.First.lib <-  function(lib,pkg,where) {
  library.dynam("plier",pkg,lib);
  require(affy,quietly=TRUE);
  require(methods,quietly=TRUE);
  where <- match(paste("package:", pkg, sep=""), search());

  cacheMetaData(as.environment(where));
  cat("PLIER V 1.0                                                  \n");
  cat("      See: http://www.affymetrix.com for the algorithm       \n");
  cat("      For details on the R port:                             \n");
  cat("      http://bioinformatics.picr.man.ac.uk                   \n");
  cat("      mailto: microarray@picr.man.ac.uk                      \n");
}
