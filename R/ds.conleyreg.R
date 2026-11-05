
#' @title Linear Regression with Conley Spatial Standard Errors
#'
#' @description
#' Performs linear regression with spatial heteroskedasticity and
#' autocorrelation consistent (Conley) standard errors on data stored
#' across one or more DataSHIELD servers.
#'
#' The function calls the server-side \code{conleyregDS()} function,
#' which uses \code{conleyreg::conleyreg()} to estimate ordinary least
#' squares (OLS) regression coefficients with spatially adjusted
#' standard errors.
#'
#' @details
#' Conley regression accounts for spatial dependence in regression
#' errors when calculating standard errors, t-statistics, and p-values.
#' The regression coefficients themselves are estimated using OLS.
#'
#' Spatial dependence is considered within the distance specified by
#' \code{dist_cutoff}. Observations beyond this distance do not
#' contribute to the spatial covariance adjustment.
#'
#' The server-side function applies DataSHIELD disclosure checks
#' using \code{nfilter.tab} and \code{nfilter.glm} before fitting
#' the regression model.
#'
#' If \code{do.checks = TRUE}, the client function checks whether
#' the dependent variable specified in \code{formula} exists on
#' the connected DataSHIELD servers and retrieves its class information.
#'
#' Results are returned separately for each datasource.
#' Currently, only \code{type = "split"} is supported.
#'
#' @param formula An R formula specifying the dependent variable
#'   and one or more independent variables for the regression model.
#'
#' @param data A character string specifying the name of the
#'   server-side spatial data object containing the variables
#'   specified in \code{formula}.
#'
#' @param dist_cutoff A positive numeric value specifying the
#'   maximum distance over which spatial correlation is accounted
#'   for when estimating standard errors.
#'   The distance units depend on the spatial distance configuration
#'   used by \code{conleyreg::conleyreg()}.
#'   Default is \code{1000}.
#'
#' @param type A character string specifying the output format.
#'   Currently, only \code{"split"} is supported, returning
#'   results separately for each datasource.
#'   Default is \code{"split"}.
#'
#' @param do.checks Logical. If \code{TRUE}, checks whether the
#'   dependent variable exists on each datasource and retrieves
#'   its class information.
#'   Default is \code{TRUE}.
#'
#' @param datasources A list of \code{DSConnection} objects
#'   representing the DataSHIELD servers.
#'   If \code{NULL}, connections are obtained using
#'   \code{datashield.connections_find()}.
#'
#' @return
#' A named list containing the results from each DataSHIELD server.
#'
#' When the model passes the server-side disclosure checks,
#' each server returns:
#'
#' \describe{
#'   \item{\code{coef}}{
#'     Regression coefficient estimates, Conley-adjusted standard
#'     errors, t-statistics, and p-values.
#'   }
#' }
#'
#' If the disclosure checks fail, the server may return:
#'
#' \describe{
#'   \item{\code{Table_valid}}{
#'     Logical value \code{FALSE}.
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
#'     Message explaining why the disclosure check failed.
#'   }
#' }
#'
#' With \code{type = "split"}, each server's output is flattened
#' into a named vector using \code{unlist()}.
#'
#' @seealso
#' \code{\link[conleyreg]{conleyreg}},
#' \code{\link[DSI]{datashield.aggregate}}
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
#' ds.conleyreg(
#'   formula = obesity_qof_2024_25 ~ dist_food_fastfood,
#'   data = "cells_object",
#'   dist_cutoff = 1000,
#'   datasources = conns)
#' @export

ds.conleyreg <- function(formula = NULL, data = NULL,  dist_cutoff = 1000,
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
  

  calltext <- paste0("conleyregDS(", deparse(formula), ", ", data, ", ", dist_cutoff, ")", collapse = "")
  
  output <- DSI::datashield.aggregate(datasources, as.symbol(calltext))
  
  ## Only type split supported now 
  if(type == "split"){
    
    result <- output
    return(result)
  }
}


