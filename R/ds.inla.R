

ds.inla <- function(formula = NULL, data = NULL, type='split',
                     do.checks = TRUE,
                     datasources=NULL){
  
  # look for DS connections
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }
  
  # ensure datasources is a list of DSConnection-class
  if(!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {methods::is(d,"DSConnection")}))))){
    stop("The 'datasources' were expected to be a list of DSConnection-class objects", call.=FALSE)
  }
  
  if(is.null(formula)){
    stop("Please provide the input formula!", call.=FALSE)
  }
  
  
  
  if(do.checks){
    
    # check if the input objects are defined in all the studies
    #isDefined(datasources, x)
    isDefined(datasources, all.vars(formula)[1])
    
    # call the internal function that checks the input object is of the same class in all studies.
    #typ1 <- checkClass(datasources, x)
    typ2 <- checkClass(datasources, all.vars(formula)[1])
  }
  
  browser()
  calltext <- paste0("inlaDS(", deparse(formula), ", ",  data, ")", collapse = "")
  
  output <- DSI::datashield.aggregate(datasources, as.symbol(calltext))
  
  ## Only type split supported now 
  if(type == "split"){
    
    result <- lapply(output, function(x) unlist(x))
    return(result)
  }
}


