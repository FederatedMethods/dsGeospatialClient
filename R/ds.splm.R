#' @title Spatial Linear Regression in DataSHIELD
#'
#' @description
#' Performs spatial linear regression on data stored across one or more
#' DataSHIELD servers using the server-side \code{splmDS()} function,
#' which relies on \code{spmodel::splm()}.
#'
#' The function estimates regression relationships while accounting
#' for spatial dependence between observations through a spatial
#' covariance structure.
#'
#' @details
#' Spatial linear regression extends ordinary linear regression by
#' explicitly modelling spatial correlation in the regression errors.
#' Unlike geographically weighted regression (GWR), it estimates
#' global regression coefficients rather than location-specific
#' coefficients.
#'
#' The spatial linear model is expressed as:
#'
#' \deqn{y_i = \beta_0 + \sum_{k=1}^{p}\beta_k x_{ik} + \epsilon_i}
#'
#' where \eqn{y_i} is the dependent variable, \eqn{x_{ik}} represents
#' the explanatory variables, and \eqn{\epsilon_i} is a spatially
#' correlated error term.
#'
#' The underlying server-side function applies DataSHIELD disclosure
#' checks before model estimation, including:
#'
#' \itemize{
#'   \item A minimum number of complete observations, determined by
#'   \code{nfilter.tab}.
#'   \item A maximum ratio of model parameters to observations,
#'   determined by \code{nfilter.glm}.
#' }
#'
#' If \code{do.checks = TRUE}, the client function checks whether
#' the dependent variable specified in \code{formula} is defined
#' on the connected DataSHIELD servers and retrieves its class
#' information.
#'
#' The regression is performed independently on each datasource.
#' No pooled or federated model is estimated.
#'
#' Currently, only \code{type = "split"} is supported.
#'
#' @param formula An R formula specifying the dependent variable
#'   and one or more independent variables for the spatial
#'   linear regression model.
#'
#' @param data A character string specifying the name of the
#'   server-side spatial data object, typically an \code{sf}
#'   object, containing the variables specified in \code{formula}
#'   and the spatial geometry required for model estimation.
#'
#' @param type A character string specifying the output format.
#'   Currently, only \code{"split"} is supported, returning
#'   results separately for each datasource.
#'   Default is \code{"split"}.
#'
#' @param do.checks Logical. If \code{TRUE}, performs preliminary
#'   checks on the dependent variable across the connected
#'   DataSHIELD servers.
#'   Default is \code{TRUE}.
#'
#' @param datasources A list of \code{DSConnection} objects
#'   representing the DataSHIELD servers.
#'   If \code{NULL}, connections are obtained using
#'   \code{datashield.connections_find()}.
#'
#' @return
#' A named list containing the spatial linear regression results
#' from each DataSHIELD server.
#'
#' When the server-side disclosure checks are satisfied, each
#' datasource returns a list containing:
#'
#' \describe{
#'   \item{\code{fixed}}{
#'     Estimated fixed-effect regression coefficients.
#'   }
#'   \item{\code{spcov}}{
#'     Estimated spatial covariance parameters describing the
#'     spatial dependence structure of the regression errors.
#'   }
#'   \item{\code{randcov}}{
#'     Estimated random-effect covariance parameters, where
#'     applicable. May be \code{NULL} if no random effects
#'     are specified.
#'   }
#' }
#'
#' If the server-side disclosure checks fail, the corresponding
#' datasource returns a list containing:
#'
#' \describe{
#'   \item{\code{Table_valid}}{
#'     Logical value \code{FALSE}, indicating that the model
#'     failed the disclosure checks.
#'   }
#'   \item{\code{Nvalid}}{
#'     Number of complete observations.
#'   }
#'   \item{\code{Nmissing}}{
#'     Number of observations containing missing values.
#'   }
#'   \item{\code{Ntotal}}{
#'     Total number of observations.
#'   }
#'   \item{\code{Warning}}{
#'     Message explaining why the disclosure checks failed.
#'   }
#' }
#'
#' Results are returned separately for each datasource without
#' flattening the underlying list structure.
#'
#' @seealso
#' \code{\link[spmodel]{splm}},
#' \code{\link[DSI]{datashield.aggregate}}
#'
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
#' ds.splm(
#'   formula = obesity_qof_2024_25 ~ dist_food_fastfood,
#'   data = "cells_object",
#'   datasources = conns
#' )
#' @export

ds.splm <- function(formula = NULL, data = NULL,  
                         type='split', do.checks = TRUE, datasources=NULL){
  
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
  
  calltext <- paste0("splmDS(", deparse(formula), ", ", data,  ")", collapse = "")
  
  output <- DSI::datashield.aggregate(datasources, as.symbol(calltext))
  
  ## Only type split supported now 
  if(type == "split"){
    
    result <- output
    return(result)
  }
}


