#' Calculate Grouped Summary Statistics Across DataSHIELD Studies
#'
#' @description
#' Computes summary statistics for a numerical variable grouped by a
#' categorical or geographical variable across one or more DataSHIELD
#' studies. Results can be combined across studies or returned separately
#' for each study. Optionally, summary statistics can be visualised as
#' choropleth maps using LSOA geographical boundaries.
#'
#' @details
#' The function calls the server-side function \code{MeanSdGpDS()} to
#' calculate grouped summary statistics without transferring individual-level
#' data to the client.
#'
#' Two output modes are supported:
#' \itemize{
#'   \item \code{"combine"}: Combines statistics across studies using
#'   group sample sizes as weights.
#'   \item \code{"split"}: Returns grouped statistics separately for
#'   each study.
#' }
#'
#' Available statistics include the mean, standard deviation (SD),
#' standard error of the mean (SEM), sample size and the 5th, 10th,
#' 25th, 50th, 75th, 90th and 95th quantiles.
#'
#' In combined mode, means and quantiles are calculated as sample-size
#' weighted averages of the study-specific estimates. The combined
#' standard deviation is calculated from the sample-size weighted
#' study-specific variances. These calculations do not reconstruct
#' pooled individual-level quantiles or the exact pooled variance.
#'
#' Combined results require identical LSOA group names and ordering
#' across studies. If these differ, the function warns and attempts
#' to switch to split mode.
#'
#' When geographical plotting is enabled, 2011 LSOA boundaries are
#' retrieved using \code{boundr} and joined to the summary statistics.
#' The region name \code{"Cheshire Merseyside"} is automatically expanded
#' into its constituent local authority districts.
#'
#' @param x Character string specifying the numerical variable to
#'   summarise, in the format \code{"data$variable"}.
#' @param y Character string specifying the grouping variable, in the
#'   format \code{"data$group"}. For geographical summaries, this should
#'   identify LSOA groups.
#' @param type Character string specifying whether to combine statistics
#'   across studies (\code{"combine"}, default) or return study-specific
#'   statistics (\code{"split"}).
#' @param do.checks Logical. If \code{TRUE}, checks that \code{x} and
#'   \code{y} are defined in the connected studies and invokes class
#'   checks. Default is \code{FALSE}.
#' @param draw.plot Logical. If \code{TRUE}, returns a geographical
#'   visualisation of the selected statistic. If \code{FALSE} (default),
#'   returns summary statistics as a list.
#' @param names_region A list of character vectors specifying local
#'   authority district names for geographical boundary retrieval.
#'   For \code{type = "combine"}, a list of length one is expected.
#'   For \code{type = "split"}, the list should contain one element
#'   per datasource.
#' @param metric Character string specifying the statistic to visualise.
#'   Supported values are \code{"Mean"}, \code{"SD"}, \code{"SEM"},
#'   \code{"5th_Quantile"}, \code{"10th_Quantile"},
#'   \code{"25th_Quantile"}, \code{"50th_Quantile"},
#'   \code{"75th_Quantile"}, \code{"90th_Quantile"} and
#'   \code{"95th_Quantile"}.
#' @param datasources A list of \code{DSConnection} objects representing
#'   the connected DataSHIELD studies. If \code{NULL}, connections are
#'   obtained using \code{datashield.connections_find()}.
#'
#' @return
#' If \code{draw.plot = FALSE}, returns a named list containing:
#' \describe{
#'   \item{Mean_gp / Mean_gp_study}{Group means, combined or by study.}
#'   \item{StDev_gp / StDev_gp_study}{Group standard deviations.}
#'   \item{Nvalid_gp / Nvalid_gp_study}{Valid sample counts by group.}
#'   \item{SEM_gp / SEM_gp_study}{Standard errors of group means.}
#'   \item{Q_gp}{Seven quantiles for each group. A data frame in
#'     combined mode or a list of matrices in split mode.}
#'   \item{Total_Nvalid}{Total number of valid observations.}
#'   \item{Total_Nmissing}{Total number of missing observations.}
#'   \item{Total_Ntotal}{Total number of observations.}
#'   \item{LSOAnames}{LSOA identifiers corresponding to the groups.}
#' }
#'
#' If \code{draw.plot = TRUE} and valid geographical parameters are
#' supplied, returns a \code{ggplot} object showing the selected
#' statistic as an LSOA choropleth map, faceted by study or combined
#' result.
#'
#' If the server reports an invalid disclosure-control cell count,
#' returns the character message
#' \code{"At least 1 cell count is 1-nfilter, please regroup"}.
#'
#' @examples
#' \dontrun{
#' library(DSOpal)
#' library(dsBase)
#' library(dsBaseClient)
#'
#'
#' builder <- DSI::newDSLoginBuilder()
#'
#' builder$append(server = "study1",
#'                url = "http://localhost:8880",
#'                user = "administrator",
#'                password = "password",
#'                table = "UPRN.greenspace_distance_metrics",
#'                profile = "grove")
#'
#' builder$append(server = "study2",
#'                url = "http://localhost:8880",
#'                user = "administrator",
#'                password = "password",
#'                table = "UPRN.greenspace_distance_metrics",
#'                profile = "grove")
#'
#' logindata <- builder$build()
#'
#' conns <- DSI::datashield.login(logins = logindata, assign = TRUE, symbol = "D")
#'
#' # Combined grouped summary statistics
#' result <- ds.summarystat(
#'   x = "D$has_asthma", y = "D$lsoa11cd",
#'   type = "combine"
#' )
#'
#' # Study-specific summary statistics
#' result_split <- ds.summarystat(
#'   x = "D$has_asthma", y = "D$lsoa11cd",
#'   type = "split"
#' )
#'
#' # Plot mean values across LSOAs
#' map <- ds.summarystat(
#'   x = "D$has_asthma", y = "D$lsoa11cd",
#'   type = "combine",
#'   draw.plot = TRUE,
#'   names_region = list(c("Liverpool", "Sefton")),
#'   metric = "Mean"
#' )
#'
#' }
#'
#' @seealso
#' \code{\link[DSI]{datashield.aggregate}}
#'
#' @export
#'
ds.summarystat <- function(
  x = NULL,
  y = NULL,
  type = 'combine',
  do.checks = FALSE,
  draw.plot = FALSE,
  names_region = NULL,
  metric = NULL,
  datasources = NULL
) {
  # look for DS connections
  if (is.null(datasources)) {
    datasources <- datashield.connections_find()
  }

  # ensure datasources is a list of DSConnection-class
  if (
    !(is.list(datasources) &&
      all(unlist(lapply(datasources, function(d) {
        methods::is(d, "DSConnection")
      }))))
  ) {
    stop(
      "The 'datasources' were expected to be a list of DSConnection-class objects",
      call. = FALSE
    )
  }

  if (is.null(x)) {
    stop("Please provide the name of the input vector!", call. = FALSE)
  }

  if (is.null(y)) {
    stop("Please provide the name of the input vector!", call. = FALSE)
  }

  if (do.checks) {
    # check if the input objects are defined in all the studies
    isDefined(datasources, x)
    isDefined(datasources, y)

    # call the internal function that checks the input object is of the same class in all studies.
    typ1 <- checkClass(datasources, x)
    typ2 <- checkClass(datasources, y)
  }

  if (!is.logical(draw.plot)) {
    stop("draw.plot should be logical")
  }
  # names of the studies
  stdnames <- names(datasources)

  # variable names

  xvarname <- strsplit(x, "\\$", perl = TRUE)[[1]][2]
  yvarname <- strsplit(y, "\\$", perl = TRUE)[[1]][2]

  # call the server side function that calculates mean and standard deviation
  # by group in each study
  calltext <- paste0("MeanSdGpDS(", x, ",", y, ")")
  output <- DSI::datashield.aggregate(datasources, as.symbol(calltext))
  numsources <- length(output)

  lsoa_names <- lapply(output, function(x) names(x$Mean_gp))

  all.tables.valid <- 1
  for (i in 1:numsources) {
    if (unlist(output[[i]]$Table_valid == FALSE)) {
      all.tables.valid <- all.tables.valid - 1
    }
  }

  if (all.tables.valid < 1) {
    warning.message <- "At least 1 cell count is 1-nfilter, please regroup"
    return(warning.message)
  }

  ## Plot matrix

  if (!type %in% c("split", "combine")) {
    warning("Invalide 'type' argument. Defaulting to type 'combine'")
    type <- "combine"
  }

  plot_matrix <- list()
  split <- type == 'split'

  ## names_region check
  # check if list and all elements are character

  if (
    is.list(names_region) &&
      length(names_region) > 0 &&
      all(is.character(unlist(names_region)))
  ) {
    names_region.valid = TRUE

    names_region <- lapply(names_region, function(i) {
      if ("Cheshire Merseyside" %in% i) {
        indx <- which(unlist(i) == "Cheshire Merseyside")
        i <- i[-indx]
        i <- c(
          i,
          "Cheshire East",
          "Cheshire West and Chester",
          "Halton",
          "Knowsley",
          "Liverpool",
          "Sefton",
          "St. Helens",
          "Warrington",
          "Wirral"
        )
      } else {
        i <- i
      }
    })

    # names_region is valid check type combination

    if (split) {
      if (numsources == 1) {
        warning(
          "type 'split', expects length of datasources to be more than 1,
              defaulting to type 'combine'"
        )
        split = FALSE
      } else {
        # check list length = numsources
        # check is all entries match when length is >1
        if (length(names_region) != numsources) {
          warning(
            paste0(
              "length of names_region should match length of datasources,
                       using only first elements in names_region to return"
            )
          )
          names_region <- lapply(1:numsources, function(i) {
            names_region[[i]] <- names_region[[1]]
          })
          names_region <- lapply(names_region, unique)
        } else {
          names_region <- lapply(names_region, unique)
        }
      }
    } else if (!split) {
      # when combine, list must have length 1
      # ifnot repeat list[[1]] for all numsources
      if (length(names_region) != 1) {
        warning(
          "length of names_region must be 1 for type 'combine',
                 using only first elements in names_region to return"
        )
        names_region = list(names_region[[1]])
        names_region <- lapply(names_region, unique)
      } else {
        names_region = list(names_region[[1]])
        names_region <- lapply(names_region, unique)
      }
    }

    shape_list <- lapply(names_region, function(region) {
      boundr::bounds(
        "lsoa",
        within_level = "lad",
        within_names = region,
        lookup_year = 2011,
        opts = boundr::boundr_options(resolution = "BFC")
      ) |>
        dplyr::select(lsoa11cd, geometry)
    })
    q.names <- c(
      "5th_Quantile",
      "10th_Quantile",
      "25th_Quantile",
      "50th_Quantile",
      "75th_Quantile",
      "90th_Quantile",
      "95th_Quantile"
    )

    if (is.null(metric) || !metric %in% c("Mean", "SD", "SEM", q.names)) {
      warning("Invalid metric type, returning available metric table")
      draw.plot = FALSE
      names_region.valid = FALSE
    } else {
      if (!split) {
        metric_map <- list("Mean" = "mean.gp", "SD" = "SD.gp", "SEM" = "SEM.gp")
      } else {
        metric_map <- list(
          "Mean" = "mean.gp.study",
          "SD" = "SD.gp.study",
          "SEM" = "SEM.gp.study"
        )
      }
    }
  } else {
    warning("invalid region names/format, returning the whole table")
    draw.plot = FALSE
    names_region.valid = FALSE
  }

  ##########

  if (all.tables.valid == 1) {
    # COMBINE
    if (type != "split") {
      numsources <- length(output)
      mean.matrix <- NULL
      sd.matrix <- NULL
      qq.matrix <- NULL
      n.matrix <- NULL
      Nvalid <- 0
      Nmissing <- 0
      Ntotal <- 0

      for (j in 1:numsources) {
        mean.matrix <- rbind(mean.matrix, as.numeric(unlist(output[[j]][2])))
        sd.matrix <- rbind(sd.matrix, as.numeric(unlist(output[[j]][3])))
        qq.matrix <- rbind(qq.matrix, as.numeric(unlist(output[[j]][4])))
        n.matrix <- rbind(n.matrix, as.numeric(unlist(output[[j]][5])))
        Nvalid <- Nvalid + as.numeric(unlist(output[[j]][6]))
        Nmissing <- Nmissing + as.numeric(unlist(output[[j]][7]))
        Ntotal <- Ntotal + as.numeric(unlist(output[[j]][8]))
      }

      var.matrix <- sd.matrix^2
      nsum.vector <- rep(1, numsources)
      # Calculate weighted means across studies in each group
      mean.gp <- (diag(t(mean.matrix) %*% n.matrix)) /
        (t(n.matrix) %*% nsum.vector)

      # Calculate weighted quantiles across studies in each group
      # Repeat each LSOA sample size for its 7 quantiles
      qq.n.matrix <- t(
        apply(n.matrix, 1, function(x) rep(x, each = 7))
      )

      # Weighted quantiles across studies
      qq.gp <- (diag(t(qq.matrix) %*% qq.n.matrix)) /
        (t(qq.n.matrix) %*% nsum.vector)

      qq.gp <- matrix(
        qq.gp,
        ncol = 7,
        byrow = TRUE
      )

      # Add names
      rownames(qq.gp) <- names(output[[1]][[4]])

      colnames(qq.gp) <- c(
        "5%_gp",
        "10%_gp",
        "25%_gp",
        "50%_gp",
        "75%_gp",
        "90%_gp",
        "95%_gp"
      )

      qq.gp <- as.data.frame(qq.gp)

      # Calculate weighted SDs across studies in each group
      var.gp <- (diag(t(var.matrix) %*% n.matrix)) /
        (t(n.matrix) %*% nsum.vector)
      SD.gp <- sqrt(var.gp)
      N.gp <- (t(n.matrix) %*% nsum.vector)
      SEM.gp <- SD.gp / sqrt(N.gp)

      # create names
      names.gp <- rep(NA, length(mean.gp))

      for (k in 1:length(mean.gp)) {
        names.gp[k] <- paste0(yvarname, "_", k)
      }

      dimnames(mean.gp) <- c(list(names.gp), list("Mean_gp"))
      dimnames(SD.gp) <- c(list(names.gp), list("SD_gp"))
      dimnames(N.gp) <- c(list(names.gp), list("Nvalid_gp"))
      dimnames(SEM.gp) <- c(list(names.gp), list("SEM_gp"))
    }

    # SPLIT
    if (type != "combine") {
      numsources <- length(output)
      mean.matrix <- NULL
      sd.matrix <- NULL
      qq.matrix <- list()
      n.matrix <- NULL
      Nvalid <- 0
      Nmissing <- 0
      Ntotal <- 0

      for (j in 1:numsources) {
        mean.matrix <- rbind(mean.matrix, as.numeric(unlist(output[[j]][2])))
        sd.matrix <- rbind(sd.matrix, as.numeric(unlist(output[[j]][3])))
        qq.matrix[[j]] <- as.numeric(unlist(output[[j]][4]))
        n.matrix <- rbind(n.matrix, as.numeric(unlist(output[[j]][5])))
        Nvalid <- Nvalid + as.numeric(unlist(output[[j]][6]))
        Nmissing <- Nmissing + as.numeric(unlist(output[[j]][7]))
        Ntotal <- Ntotal + as.numeric(unlist(output[[j]][8]))
      }
      var.matrix <- sd.matrix^2
      mean.gp.study <- t(mean.matrix)
      SD.gp.study <- t(sd.matrix)
      N.gp.study <- t(n.matrix)
      SEM.gp.study <- SD.gp.study / sqrt(N.gp.study)

      qq.gp.study <- lapply(qq.matrix, function(i) {
        matrix(i, ncol = 7, byrow = TRUE)
      })

      qq.gp.study <- lapply(qq.gp.study, function(i) {
        colnames(i) <- c(
          "5%_gp",
          "10%_gp",
          "25%_gp",
          "50%_gp",
          "75%_gp",
          "90%_gp",
          "95%_gp"
        )
        i
      })

      # create names
      names.gp <- rep(NA, dim(mean.gp.study)[1])

      for (k in 1:dim(mean.gp.study)[1]) {
        names.gp[k] <- paste0(yvarname, "_", k)
      }

      names.study <- names(datasources)

      dimnames(mean.gp.study) <- c(list(names.gp), list(names.study))
      dimnames(SD.gp.study) <- c(list(names.gp), list(names.study))
      dimnames(N.gp.study) <- c(list(names.gp), list(names.study))
      dimnames(SEM.gp.study) <- c(list(names.gp), list(names.study))
    }

    if (type == "combine") {
      # if type is combine check whether lsoa levels are same in the two studies
      names.identical <- all(vapply(
        lsoa_names,
        identical,
        logical(1),
        y = lsoa_names[[1]]
      ))

      if (!names.identical) {
        warning(
          "LSOA levels across studies are not the same. ",
          "Cannot proceed with type = 'combine'. ",
          "Defaulting to type = 'split'."
        )
        type <- "split"
      } else {
        lsoanames <- lsoa_names[[1]]
        if (draw.plot) {
          ## If plot is TRUE please select which metric to plot else table is returned

          if (metric %in% q.names) {
            sel.metric.data <- data.frame(lsoa11cd = rownames(qq.gp))
            sel.metric.data$value <- qq.gp[, which(metric == q.names)]
            sel.metric.data$server <- 'combine'
          } else {
            sel.metric <- get(metric_map[[metric]])

            sel.metric.data <- data.frame(
              lsoa11cd = lsoanames,
              value = as.numeric(sel.metric),
              server = 'combine'
            )
          }
          shape_sf <- shape_list[[1]] |>
            dplyr::left_join(sel.metric.data, by = "lsoa11cd") |>
            sf::st_as_sf()
        } else {
          result <- list(
            mean.gp,
            SD.gp,
            N.gp,
            SEM.gp,
            qq.gp,
            Nvalid,
            Nmissing,
            Ntotal,
            lsoa_names[[1]]
          )
          names(result) <- list(
            "Mean_gp",
            "StDev_gp",
            "Nvalid_gp",
            "SEM_gp",
            "Q_gp",
            "Total_Nvalid",
            "Total_Nmissing",
            "Total_Ntotal",
            "LSOAnames"
          )
          return(result)
        }
      }
    }

    if (type == "split") {
      if (draw.plot) {
        if (metric %in% q.names) {
          sel.metric.data <-
            lapply(1:length(qq.gp.study), function(i) {
              data.frame(
                server = names(datasources)[i],
                value = unlist(qq.gp.study[[i]][, which(metric == q.names)]),
                lsoa11cd = lsoa_names[[i]],
                row.names = NULL
              )
            })
        } else {
          sel.metric <- get(metric_map[[metric]])

          sel.metric.data <-
            lapply(seq_len(ncol(sel.metric)), function(i) {
              data.frame(
                server = colnames(sel.metric)[i],
                value = sel.metric[, i],
                lsoa11cd = lsoa_names[[i]],
                row.names = NULL
              )
            })
        }

        shape_sf <- do.call(
          rbind,
          lapply(1:numsources, function(i) {
            shape_list[[i]] |>
              dplyr::left_join(sel.metric.data[[i]], by = "lsoa11cd") |>
              sf::st_as_sf()
          })
        )
      } else {
        lsoanames <- matrix(unlist(lsoa_names), ncol = numsources)
        dimnames(lsoanames) <- dimnames(mean.gp.study)
        result <- list(
          mean.gp.study,
          SD.gp.study,
          N.gp.study,
          SEM.gp.study,
          qq.gp.study,
          Nvalid,
          Nmissing,
          Ntotal,
          lsoanames
        )
        names(result) <- list(
          "Mean_gp_study",
          "StDev_gp_study",
          "Nvalid_gp_study",
          "SEM_gp_study",
          "Q_gp",
          "Total_Nvalid",
          "Total_Nmissing",
          "Total_Ntotal",
          "LSOAnames"
        )
        return(result)
      }
    }
  }

  if (draw.plot) {
    plotresult <- ggplot2::ggplot(shape_sf) +
      ggplot2::geom_sf(
        ggplot2::aes(fill = value)
      ) +
      ggplot2::scale_fill_gradientn(
        colours = grDevices::colorRampPalette(
          c(
            "#440154",
            "#414487",
            "#2A788E",
            "#22A884",
            "#7AD151",
            "#FDE725"
          )
        )(100),
        na.value = "grey90"
      ) +
      ggplot2::labs(fill = metric) +
      ggplot2::facet_wrap(~server) +
      ggplot2::theme_minimal()

    return(plotresult)
  }
}
#ds.meanSdGp
