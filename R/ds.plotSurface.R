#' Plot Spatial Relative Risk Surface
#'
#' @description
#' Creates a spatial relative risk surface using kernel density
#' estimates of cases and negative observations. The function calls
#' the corresponding server-side function, \code{plotSpatialDS()},
#' through DataSHIELD and uses the returned data to construct a
#' \code{ggplot2} map displaying the log relative risk surface.
#'
#' @details
#' The function requests spatial data from the server and generates
#' a smoothed relative risk surface within a specified geographical
#' window.
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Calls the server-side \code{plotSpatialDS()} function
#'   using \code{DSI::datashield.aggregate()}.
#'   \item Retrieves the case observations, negative observations,
#'   and basemap returned by the server.
#'   \item Defines a geographical display window with a 1 km buffer
#'   and generates a regular spatial grid at 50 m resolution.
#'   \item Calculates Gaussian kernel density estimates for cases
#'   and negative observations using a bandwidth of 800 m.
#'   \item Computes the log relative risk as the logarithm of the
#'   ratio between the case and negative kernel density estimates.
#'   \item Suppresses grid cells containing fewer than five negative
#'   observations within an 800 m radius.
#'   \item Overlays the relative risk surface on the basemap using
#'   a diverging colour scale, with a dashed geographical window
#'   boundary and a 1 km scale bar.
#' }
#'
#' Cases are selected by the server-side function according to the
#' minimum number of neighbouring cases specified by \code{k}.
#' The neighbourhood is defined by the supplied network distance
#' matrix \code{D}.
#'
#' The relative risk surface is calculated using normalised kernel
#' density estimates. Positive log relative risk values indicate
#' higher relative case density, while negative values indicate
#' higher relative negative-observation density.
#'
#' The displayed colour scale is limited to log relative risk values
#' between \code{-log(4)} and \code{log(4)}, with values outside
#' this range assigned the corresponding endpoint colour.
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
#'   be included in the relative risk calculation. Defaults to 5.
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
#' @param datasources DataSHIELD connection object specifying
#'   the server(s) on which the analysis is performed.
#'
#' @return
#' A \code{ggplot} object displaying a smoothed spatial log
#' relative risk surface over the basemap. The plot includes
#' a diverging colour scale, a dashed geographical window
#' boundary, and a 1 km scale bar.
#'
#' @seealso
#' \code{plotSpatialDS()}, \code{ds.plotSpatial()}
#'
#' @examples
#' \dontrun{
#' ds.plotSurface(
#'   D = "D_object",
#'   study = "study_object",
#'   mapping = "mapping_object",
#'   pos = "pos_object",
#'   basemap = "basemap_object",
#'   area = "east",
#'   case_var = "case",
#'   k = 5,
#'   xmin = 334000,
#'   ymin = 388500,
#'   xmax = 339000,
#'   ymax = 393500,
#'   datasources = conns
#' )
#' }
#'
#' @export




ds.plotSurface <- function(
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
  

  grid <- CJ(
    x = seq(view["xmin"] + 25, view["xmax"] - 25, by = 50),
    y = seq(view["ymin"] + 25, view["ymax"] - 25, by = 50)
  )
  
  grid_sum <- function(px, py, w, h, kernel = c("gaussian", "count")) {
    kernel <- match.arg(kernel)
    out <- numeric(nrow(grid))
    chunks <- split(seq_len(nrow(grid)), ceiling(seq_len(nrow(grid)) / 2000))
    for (grid_rows in chunks) {
      sq_dist <- outer(grid$x[grid_rows], px, "-")^2 + outer(grid$y[grid_rows], py, "-")^2
      if (kernel == "gaussian") {
        kernel_weights <- exp(-sq_dist / (2 * h^2))
      } else {
        kernel_weights <- sq_dist <= h^2
      }
      out[grid_rows] <- drop(kernel_weights %*% w)
    }
    out
  }
  
  h <- 800
  case_density <- grid_sum(cases$x, cases$y, rep(1, nrow(cases)), h) / nrow(cases)
  
  
  negative_density <- grid_sum(negatives$x, negatives$y, rep(1, nrow(negatives)), h) / nrow(negatives)
  negatives_near <- grid_sum(negatives$x, negatives$y, rep(1, nrow(negatives)), h, "count")
  
  risk <- grid[, .(x, y)]
  risk[, value := log(pmax(case_density, 1e-12) / pmax(negative_density, 1e-12))]
  risk <- risk[negatives_near >= 5]
  
  diverging <- c("#104281", "#3987e5", "#f0efec", "#e34948", "#8f1d1d")
  
  map_base <- ggplot2::ggplot() +
    annotation_raster(basemap$img, basemap$xmin, basemap$xmax, basemap$ymin, basemap$ymax) +
    annotate("rect", xmin = window[["xmin"]], xmax = window[["xmax"]],
             ymin = window[["ymin"]], ymax = window[["ymax"]],
             fill = NA, colour = "grey55", linetype = "dashed", linewidth = 0.3) +
    annotate("segment", x = bar_x, xend = bar_x + 1000, y = bar_y, yend = bar_y, linewidth = 0.8) +
    annotate("text", x = bar_x + 500, y = bar_y + 110, label = "1 km", size = 3) +
    coord_equal(xlim = view[c("xmin", "xmax")], ylim = view[c("ymin", "ymax")], expand = FALSE) +
    theme_void()
  
  
  p <- map_base() +
    geom_raster(data = risk, aes(x, y, fill = value), alpha = 0.5) +
    scale_fill_gradientn(colours = diverging, limits = c(-log(4), log(4)),
                         oob = scales::squish, name = "log relative risk")
  
  
  return(p)
}


