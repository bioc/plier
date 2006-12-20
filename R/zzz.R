.First.lib <-  function(lib,pkg,where) {
  library.dynam("plier",pkg,lib);
  require(affy,quietly=TRUE);
  require(methods,quietly=TRUE);
  where <- match(paste("package:", pkg, sep=""), search());

  cacheMetaData(as.environment(where));
}
