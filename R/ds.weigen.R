#' @title Calculate Spatial Eigenvectors in DataSHIELD
#'
#' @description
#' Calculates spatial eigenvectors from geographical data using
#' \code{spmoran::weigen()} through the server-side
#' \code{weigenDS()} function.
#'
#' The function enables spatial eigenvector analysis across one or
#' more DataSHIELD servers without transferring individual-level
#' spatial data to the client.
#'
#' @details
#' Spatial eigenvectors represent spatial patterns among geographical
#' observations and can be used to account for spatial autocorrelation
#' in statistical models.
#'
#' The function calls \code{weigenDS()} independently on each
#' DataSHIELD server. The server-side function extracts centroid
#' coordinates from the supplied spatial dataset and calculates
#' spatial eigenvectors using \code{spmoran::weigen()}.
#'
#' The spatial neighbourhood structure is determined by the
#' \code{weights} argument. For example, \code{"knn"} constructs
#' a k-nearest-neighbour structure, while \code{"tri"} uses
#' Delaunay triangulation.
#'
#' The \code{k} argument specifies the number of nearest neighbours
#' when \code{weights = "knn"}. The \code{threshold} and \code{enum}
#' arguments control the selection and maximum number of spatial
#' eigenvectors.
#'
#' When \code{do.checks = TRUE}, the function checks whether the
#' input object exists on the specified servers and examines its
#' class before executing the analysis.
#'
#' The server-side function applies DataSHIELD disclosure checks
#' based on the minimum observation threshold (\code{nfilter.tab}).
#' Only the server-generated summaries are returned to the client.
#' The observation-level eigenvector matrix is not returned.
#'
#' Currently, only \code{type = "split"} is supported. Results
#' are returned separately for each DataSHIELD server, without
#' pooling estimates across studies.
#'
#' @param x A character string specifying the name of a server-side
#'   \code{sf} object containing the spatial observations.
#'   The object must be available on each selected DataSHIELD server.
#'
#' @param weights A character string specifying the method used to
#'   construct the spatial neighbourhood structure. Supported values
#'   include \code{"knn"} for k-nearest neighbours and \code{"tri"}
#'   for Delaunay triangulation. Default is \code{"knn"}.
#'
#' @param k An integer specifying the number of nearest neighbours
#'   when \code{weights = "knn"}. Default is \code{4}.
#'
#' @param threshold A numeric value specifying the relative eigenvalue
#'   threshold used to select spatial eigenvectors.
#'   Default is \code{0.25}.
#'
#' @param enum An integer specifying the maximum number of spatial
#'   eigenvectors to extract. Default is \code{100}.
#'
#' @param type A character string specifying the output format.
#'   Currently, only \code{"split"} is supported, returning
#'   results separately for each DataSHIELD server.
#'   Default is \code{"split"}.
#'
#' @param do.checks A logical value indicating whether to check
#'   that the input object exists and examine its class across
#'   the selected DataSHIELD servers. Default is \code{TRUE}.
#'
#' @param datasources A list of \code{DSConnection} objects specifying
#'   the DataSHIELD servers on which the analysis should be performed.
#'   If \code{NULL}, the currently active DataSHIELD connections
#'   are used.
#'
#' @return
#' A named list containing the results returned by each DataSHIELD
#' server when \code{type = "split"}.
#'
#' For successful analyses, each server returns:
#' \describe{
#'   \item{\code{n_eigenvectors}}{
#'     Number of retained spatial eigenvectors.
#'   }
#'   \item{\code{eigenvalue_summary}}{
#'     Named numeric vector containing the first quartile
#'     (\code{Q1}), median (\code{Median}), mean (\code{Mean}),
#'     and third quartile (\code{Q3}) of the retained eigenvalues.
#'   }
#'   \item{\code{spectrum}}{
#'     Numeric vector containing the retained eigenvalues,
#'     if included in the server-side output.
#'   }
#' }
#'
#' If the server-side disclosure check fails, the returned list
#' contains \code{Table_valid = FALSE}, together with the
#' observation counts and a warning message, as provided by
#' \code{weigenDS()}.
#'
#' @note
#' The analysis is performed independently on each DataSHIELD server.
#' No individual spatial observations or eigenvector matrices are
#' returned to the client.
#'
#' The server-side minimum observation-count check does not, by
#' itself, guarantee that all returned summaries are disclosure-safe.
#' In particular, releasing the complete eigenvalue spectrum or
#' exact counts may require additional disclosure assessment.
#'
#' @seealso
#' \code{\link[spmoran]{weigen}},
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
#'ds.weigen(
#'   x = "cells_object",
#'   weights = "knn",
#'   k = 4,
#'   threshold = 0.25,
#'   enum = 100,
#'   datasources = conns
#' )
#'
#' @export

ds.weigen <- function(x = NULL, weights = "knn", 
                      k = 4, threshold = 0.25, 
                      enum = 100, type='split',
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
  
  if(is.null(x)){
    stop("Please provide the input vector 'x', an sf object!", call.=FALSE)
  }
  
  
  
  if(do.checks){
    
    # check if the input objects are defined in all the studies
    isDefined(datasources, x)

    # call the internal function that checks the input object is of the same class in all studies.
    typ1 <- checkClass(datasources, x)
  }
  

  calltext <- paste0("weigenDS(", x, ", ", shQuote(weights), ", ", k, ", ", threshold, ",", enum, ")")
  
  output <- DSI::datashield.aggregate(datasources, as.symbol(calltext))
  
  ## Only type split supported now 
  if(type == "split"){
    result <- output
    return(result)
  }
}


