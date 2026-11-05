#' Calculate Local Moran's I
#'
#' Calculates Local Moran's I for a numeric variable using a spatial weights
#' object across one or more DataSHIELD studies.
#'
#' @param x Character string specifying the name of the numeric vector on the
#'   server to be analysed.
#' @param listw Character string specifying the name of the spatial weights
#'   object of class \code{listw} on the server.
#' @param type Character string specifying the output format. Currently only
#'   \code{"split"} is supported.
#' @param vars Optional variable specification. Currently unused.
#' @param do.checks Logical. If \code{TRUE}, checks that the input objects are
#'   defined and have consistent classes across the studies.
#' @param nsim Numeric value specifying the number of simulations used for
#'   calculating the Monte Carlo p-value. Default is 199.
#' @param datasources A list of \code{DSConnection} objects. If \code{NULL},
#'   active DataSHIELD connections are used.
#'
#' @return A named list with one element per study. Each study contains a
#'   named numeric vector with:
#'   \itemize{
#'     \item \code{moran_I}: Moran's I statistic.
#'     \item \code{z}: standardised Moran's I statistic.
#'   }
#'
#' @examples
#' \dontrun{
#' DSOpal::Opal()
#' builder <- DSI::newDSLoginBuilder()
#' builder$append(
#'   server = "mersey-demo",
#'   url = "http://localhost:8880",
#'   user = "administrator",
#'   password = "password",
#'   driver = "OpalDriver",
#'   profile = "ds-geospatial"
#' )
#' 
#' logindata <- builder$build()
#' conns <- DSI::datashield.login(logins = logindata)
#' DSI::datashield.assign.resource(
#' conns,
#' symbol = "cells",
#' resource = paste0(project, ".", "cells")
#' )
#' DSI::datashield.assign.resource(
#' conns,
#' symbol = "W",
#' resource = paste0(project, ".", "W_D_5km")
#' )
#' ## 2.8 Materialise arbitrary R objects ----
#' DSI::datashield.assign.resource(
#'   conns,
#'   symbol = "cells_resource",
#'   resource = paste0(project, ".cells")
#' )
#' 
#' DSI::datashield.assign.expr(
#'   conns,
#'   symbol = "cells_object",
#'   expr = quote(as.resource.object(cells_resource))
#' )
#' ## 2.9 Map specific columns ----
#' DSI::datashield.assign.expr(
#'   conns,
#'   symbol = "obesity_qof_2024_25",
#'   expr = quote(cells_object$obesity_qof_2024_25)
#' )
#' DSI::datashield.assign.expr(
#'   conns,
#'   symbol = "W_object",
#'   expr = quote(as.resource.object(W))
#' )
#' DSI::datashield.assign.expr(
#'   conns,
#'   symbol = "W_band800",
#'   expr = quote(W_object$band800)
#' )
#' 
#' ds.moransI(
#'   x = "obesity_qof_2024_25",
#'   listw = "W_band800",
#'   nsim = 199,
#'   datasources = conns
#' )
#' 
#'
#' @export

ds.localmoran <- function(x = NULL, listw = NULL, type='split',
                       vars = NULL, do.checks = TRUE,
                       nsim = 199, datasources=NULL){
  
  # look for DS connections
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }
  
  # ensure datasources is a list of DSConnection-class
  if(!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {methods::is(d,"DSConnection")}))))){
    stop("The 'datasources' were expected to be a list of DSConnection-class objects", call.=FALSE)
  }
  
  if(is.null(x)){
    stop("Please provide the input vector 'x', a numeric!", call.=FALSE)
  }
  
  if(is.null(listw)){
    stop("Please provide the input vector 'listw'!, a list object", call.=FALSE)
  }
  
  if(do.checks){
    
    # check if the input objects are defined in all the studies
    isDefined(datasources, x)
    isDefined(datasources, listw)
    
    # call the internal function that checks the input object is of the same class in all studies.
    typ1 <- checkClass(datasources, x)
    typ2 <- checkClass(datasources, listw)
  }
  
  if(!is.numeric(nsim) & length(nsim != 1)){
    stop("nsim should be numeric of length 1")
  }
  calltext <- paste0("Localmoran(", x, ", ", listw, ")")

  output <- DSI::datashield.aggregate(datasources, as.symbol(calltext))

  ## Only type split supported now 
  if(type == "split"){
    return(output)
  }
  
}


