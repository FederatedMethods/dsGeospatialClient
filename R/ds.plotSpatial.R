#' Plot Spatial Distribution of Cases and Negatives
#'
#' @description
#' Creates a spatial map of cases and negative observations using
#' precomputed display coordinates and a basemap. The function calls
#' the corresponding server-side function, \code{plotSpatialDS()},
#' through DataSHIELD and uses the returned data to construct a
#' \code{ggplot2} map.
#'
#' @details
#' The function requests spatial data from the server and generates
#' a map showing cases and negative observations within a specified
#' geographical window.
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Calls the server-side \code{plotSpatialDS()} function
#'   using \code{DSI::datashield.aggregate()}.
#'   \item Retrieves the case observations, negative observations,
#'   and basemap returned by the server.
#'   \item Defines a geographical display window with a 1 km buffer.
#'   \item Overlays the observations on the basemap, displaying
#'   negative observations as blue points and cases as orange triangles.
#'   \item Adds a dashed boundary indicating the specified window
#'   and a 1 km scale bar.
#' }
#'
#' Cases are selected by the server-side function according to the
#' minimum number of neighbouring cases specified by \code{k}.
#' The neighbourhood is defined by the supplied network distance
#' matrix \code{D}.
#'
#' @param D Character string specifying the name of the network
#'   distance matrix assigned on the DataSHIELD server.
#'
#' @param study Character string specifying the name of the study
#'   dataset assigned on the DataSHIELD server. The dataset must
#'   contain the binary outcome variable specified by \code{case_var}.
#'
#' @param mapping Character string specifying the name of the
#'   server-side mapping table linking study observations to
#'   distance matrix row identifiers.
#'
#' @param pos Character string specifying the name of the
#'   server-side object containing precomputed display coordinates.
#'
#' @param basemap Character string specifying the name of the
#'   server-side basemap object containing the raster image
#'   and its spatial extent.
#'
#' @param area Character string specifying the geographical area
#'   to display. Defaults to \code{"east"}.
#'
#' @param case_var Character string specifying the binary outcome
#'   variable, where 1 represents a case and 0 represents a
#'   negative observation. Defaults to \code{"case"}.
#'
#' @param k Integer specifying the minimum number of other cases
#'   required within the network neighbourhood for a case to
#'   be displayed. Defaults to 5.
#'
#' @param xmin Numeric value specifying the minimum easting
#'   of the geographical window. Defaults to 334000.
#'
#' @param ymin Numeric value specifying the minimum northing
#'   of the geographical window. Defaults to 388500.
#'
#' @param xmax Numeric value specifying the maximum easting
#'   of the geographical window. Defaults to 339000.
#'
#' @param ymax Numeric value specifying the maximum northing
#'   of the geographical window. Defaults to 393500.
#'
#' @return
#' A \code{ggplot} object displaying the spatial distribution
#' of cases and negative observations over the basemap.
#' The plot includes a dashed geographical window boundary
#' and a 1 km scale bar.
#'
#' @seealso
#' \code{plotSpatialDS()}
#'ds.plotSpatial(
#' D = "D_object",
#' study = "study_object",
#' mapping = "mapping_object",
#' pos = "pos_object",
#' basemap = "basemap_object",
#' area = "east",
#' case_var = "case",
#' k = 5,
#' xmin = 334000, 
#' ymin = 388500, 
#' xmax = 339000, 
#' ymax = 393500, 
#' datasources = conns)
#' @export
#' 
ds.plotSpatial <- function(
    D,
    study,
    mapping,
    pos,
    basemap,
    area = "east",
    case_var = "case",
    k = 5,
    xmin = 334000, 
    ymin = 388500, 
    xmax = 339000, 
    ymax = 393500, 
    datasources
) {
  
  # Select the requested area
  
  window <- c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)
  view <- window + c(-1000, -1000, 1000, 1000)
  
  
  calltext <- paste0(
    "plotSpatialDS(", D, ",", study, ",", mapping, ",", pos, ",",
    basemap, ",",
    encodeString(area, quote = '"'), ",",
    encodeString(case_var, quote = '"'), ",",
    k, ")"
  )

  output <- DSI::datashield.aggregate(datasources, as.symbol(calltext))
  
  basemap <- output[[1]]$basemap
  cases <- output[[1]]$cases
  negatives <- output[[1]]$negatives
  
  # 5. Plot on basemap
 
  bar_x <- view[["xmin"]] + 300
  bar_y <- view[["ymin"]] + 300
  
  map_base <- ggplot2::ggplot() +
    annotation_raster(basemap$img, basemap$xmin, basemap$xmax, basemap$ymin, basemap$ymax) +
    annotate("rect", xmin = window[["xmin"]], xmax = window[["xmax"]],
             ymin = window[["ymin"]], ymax = window[["ymax"]],
             fill = NA, colour = "grey55", linetype = "dashed", linewidth = 0.3) +
    annotate("segment", x = bar_x, xend = bar_x + 1000, y = bar_y, yend = bar_y, linewidth = 0.8) +
    annotate("text", x = bar_x + 500, y = bar_y + 110, label = "1 km", size = 3) +
    coord_equal(xlim = view[c("xmin", "xmax")], ylim = view[c("ymin", "ymax")], expand = FALSE) +
    theme_void()
  
  p <- map_base +
    ggplot2::geom_point(
      data = negatives,
      ggplot2::aes(x = x, y = y),
      colour = "#2a78d6",
      size = 1.4
    ) +
    ggplot2::geom_point(
      data = cases,
      ggplot2::aes(x = x, y = y),
      colour = "#eb6834",
      shape = 17,
      size = 2.2
    ) +
    ggplot2::theme_void()
  
  return(p)
}


