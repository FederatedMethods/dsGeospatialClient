#' @title Fit a Spatial Lag of X (SLX) Regression Model
#'
#' @description
#' Fits a Spatial Lag of X (SLX) regression model on DataSHIELD servers
#' using the server-side function \code{lmSLXDS()}. The function
#' estimates the association between an outcome and both local and
#' neighbouring explanatory variables without transferring individual-level
#' data to the client.
#'
#' @details
#' The SLX model incorporates spatially lagged explanatory variables
#' and is defined as:
#'
#' \deqn{y = X\beta + WX\theta + \epsilon}
#'
#' where \eqn{y} is the dependent variable, \eqn{X} is the matrix of
#' explanatory variables, \eqn{W} is the spatial weights matrix,
#' \eqn{\beta} represents the coefficients of local explanatory variables,
#' and \eqn{\theta} represents the coefficients of spatially lagged
#' explanatory variables.
#'
#' The function calls \code{lmSLXDS()} on each DataSHIELD server,
#' where disclosure checks are applied before fitting the model using
#' \code{spatialreg::lmSLX()}.
#'
#' Currently, only \code{type = "split"} is supported, returning
#' results separately for each server.
#'
#' @param formula An R formula specifying the dependent variable and
#'   explanatory variables.
#' @param data Character string specifying the name of the data frame
#'   on the DataSHIELD servers.
#' @param listw Character string specifying the name of a server-side
#'   spatial weights object of class \code{listw}.
#' @param type Character string specifying how results are returned.
#'   Currently, only \code{"split"} is supported. Defaults to
#'   \code{"split"}.
#' @param do.checks Logical. If \code{TRUE}, checks whether the spatial
#'   weights object exists on the DataSHIELD servers and retrieves its
#'   class. Defaults to \code{TRUE}.
#' @param datasources A list of \code{DSConnection} objects representing
#'   the DataSHIELD connections. If \code{NULL}, existing connections
#'   are retrieved using \code{datashield.connections_find()}.
#'
#' @return
#' A named list containing results from each DataSHIELD server.
#' When \code{type = "split"}, each server's output is converted
#' to a vector using \code{unlist()}.
#'
#' If the server-side disclosure checks pass, the results contain
#' estimated regression coefficients, including coefficients for
#' spatially lagged explanatory variables.
#'
#' If the disclosure checks fail, the server returns information
#' about the failed checks instead of fitting the model.
#'
#' @seealso
#' \code{\link[spatialreg]{lmSLX}}
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
#' ds.lmSLX(obesity_qof_2024_25 ~ dist_food_fastfood + obesity_qof_2023_24,
#' data = "cells_object", listw = "W_band800", datasources = conns)
#' 
#' @export

ds.gwrbasic <- function(formula = NULL, data = NULL, dMat = NULL, type='split',
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

  dMat_arg <- ifelse((is.null(dMat)), "NULL", dMat)
  calltext <- paste0("gwrbasicDS(", deparse(formula), ", ", data, ", ", dMat_arg,")", collapse = "")

  output <- DSI::datashield.aggregate(datasources, as.symbol(calltext))
  
  ## Only type split supported now 
  if(type == "split"){
    result <- output
    return(result)
  }
}


