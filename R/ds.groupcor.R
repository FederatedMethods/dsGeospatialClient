#' Calculate Grouped Pearson Correlations Across DataSHIELD Studies
#'
#' @description
#' Calculates Pearson correlation coefficients between two numerical
#' variables within geographical groups across one or more DataSHIELD
#' studies. Correlations can be calculated separately for each study
#' or combined across studies using aggregated summary statistics.
#' Results can be returned as a table or visualised as choropleth maps
#' using Lower Layer Super Output Area (LSOA) boundaries.
#'
#' @details
#' The function calls the server-side function \code{groupcorDS()} to
#' obtain aggregated statistics required for calculating correlations
#' without transferring individual-level data to the client.
#'
#' Pearson correlation coefficients are calculated from covariance
#' matrices constructed using sums of products, sums of variables
#' and complete-case observation counts. The covariance matrices
#' are converted to correlation matrices using \code{stats::cov2cor()}.
#'
#' Two analysis modes are supported:
#' \itemize{
#'   \item \code{"split"}: Calculates correlation coefficients
#'   independently for each study and geographical group.
#'   \item \code{"combine"}: Combines aggregated sufficient statistics
#'   across studies for each common geographical group before
#'   calculating the correlation coefficient.
#' }
#'
#' In combined mode, the function pools sums of products, variable
#' sums and complete-case counts across studies. This produces a
#' correlation based on the combined observations rather than an
#' average of study-specific correlation coefficients.
#'
#' The function is intended for numerical variables, including
#' continuous, discrete and binary numerical variables. Character
#' variables and unordered categorical variables are not suitable
#' for direct Pearson correlation analysis. The grouping variable
#' can represent categorical geographical identifiers such as LSOAs.
#'
#' For combined analysis, only LSOAs present in all connected studies
#' are included. The number of common LSOAs must satisfy the
#' disclosure-control threshold obtained from the server.
#'
#' When geographical plotting is enabled, 2011 LSOA boundaries are
#' retrieved using \code{boundr}. Correlation coefficients are joined
#' to the corresponding geographical boundaries and displayed using
#' \code{ggplot2}.
#'
#' The region name \code{"Cheshire Merseyside"} is automatically
#' expanded into its constituent local authority districts.
#'
#' @param x Character string specifying the name of the first
#'   numerical variable on the DataSHIELD server.
#'
#' @param y Character string specifying the name of the second
#'   numerical variable on the DataSHIELD server.
#'
#' @param index Character string specifying the grouping variable,
#'   typically containing LSOA geographical identifiers.
#'
#' @param type Character string specifying the analysis mode.
#'   Either \code{"split"} (default), which returns study-specific
#'   correlations, or \code{"combine"}, which calculates pooled
#'   correlations across studies.
#'
#' @param datasources A list of \code{DSConnection} objects
#'   representing the connected DataSHIELD studies. If \code{NULL},
#'   connections are obtained using
#'   \code{datashield.connections_find()}.
#'
#' @param do.checks Logical. If \code{TRUE} (default), checks whether
#'   the input variables are defined in the connected studies and
#'   invokes class checks for \code{x} and \code{y}.
#'
#' @param names_region A list of character vectors specifying local
#'   authority district names for geographical boundary retrieval.
#'   For \code{type = "combine"}, a list of length one is expected.
#'   For \code{type = "split"}, one list element per datasource is
#'   expected. The function signature defaults to \code{"Liverpool"},
#'   but a list must be supplied to pass input validation.
#'
#' @param draw.plot Logical. If \code{TRUE} (default), returns a
#'   choropleth map of correlation coefficients. If \code{FALSE},
#'   returns a data frame containing the calculated correlations.
#'
#' @return
#' If \code{draw.plot = FALSE}, returns a data frame containing:
#' \describe{
#'   \item{lsoa11cd}{LSOA geographical identifier.}
#'   \item{complete.counts}{Number of complete observation pairs
#'     used to calculate the correlation.}
#'   \item{correlation.coef}{Pearson correlation coefficient between
#'     the two numerical variables.}
#'   \item{server}{Name of the contributing DataSHIELD study, or
#'     \code{"combine"} for pooled results.}
#' }
#'
#' In split mode, the returned data frame also includes geographical
#' boundary information as an \code{sf} object.
#'
#' If \code{draw.plot = TRUE}, returns a \code{ggplot} object
#' displaying correlation coefficients across LSOAs, faceted
#' by study or combined result.
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
#' # Calculate study-specific correlations
#' result <- ds.groupcor(
#'   x = "D$has_asthma", y = "D$distance_local_greenspace",
#'   index = "D$lsoa11cd",
#'   type = "split",
#'   names_region = list("Liverpool", "Sefton"),
#'   draw.plot = FALSE
#' )
#'
#' # Calculate combined correlations across studies
#' combined <- ds.groupcor(
#'   x = "D$has_asthma", y = "D$distance_local_greenspace",
#'   index = "D$lsoa11cd",
#'   type = "combine",
#'   names_region = list(c("Liverpool", "Sefton")),
#'   draw.plot = FALSE
#' )
#'
#' # Visualise combined correlations
#' map <- ds.groupcor(
#'   x = "D$has_asthma", y = "D$distance_local_greenspace",
#'   index = "D$lsoa11cd",
#'   type = "combine",
#'   names_region = list("Liverpool"),
#'   draw.plot = TRUE
#' )
#' print(map)
#' }
#'
#' @seealso
#' \code{\link[stats]{cor}},
#' \code{\link[stats]{cov2cor}},
#' \code{\link[DSI]{datashield.aggregate}}
#'
#' @export

ds.groupcor <- function(
  x = NULL,
  y = NULL,
  index = NULL,
  type = "split",
  datasources = NULL,
  do.checks = TRUE,
  names_region = list("Liverpool"),
  draw.plot = FALSE
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
    stop(
      "x=NULL. Please provide the name of a matrix or dataframe or the names of two numeric vectors!",
      call. = FALSE
    )
  }

  if (is.null(y)) {
    stop("Please provide the name of the input vector!", call. = FALSE)
  }

  if (is.null(index)) {
    stop("Please provide the name of the input vector!", call. = FALSE)
  }

  # check the type of the input objects
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

  calltext <- paste0("groupcorDS(", x, ",", y, ",", index, ")")

  output <- DSI::datashield.aggregate(datasources, as.symbol(calltext))
  nfilter.tab <- output[[1]][["Nfilter.tab"]]
  output <- lapply(output, function(x) {
    x[-length(x)]
  })

  # name of the studies to be used in the output
  stdnames <- names(datasources)

  if (!type %in% c("split", "combine")) {
    warning("Invalide 'type' argument. Defaulting to type 'combine'")
    type <- "combine"
  }
  if (!is.list(names_region)) {
    stop("provide names_region as a list")
  }
  plot_matrix <- list()
  split <- type == 'split'
  numsources <- length(datasources)

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
    } else {
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
  } else {
    warning("invalid region names/format, returning the whole table")
    draw.plot = FALSE
    names_region.valid = FALSE
  }

  ## Loading map

  # call the server side function
  if (type == "split") {
    names(shape_list) <- stdnames
    covariance <- list()
    sqrt.diag <- list()
    correlation <- list()
    splitresults <- list()
    results <- list()
    for (i in 1:length(output)) {
      for (j in 1:length(output[[i]])) {
        # -1 to exclude [[nfiltertab]]
        covariance[[j]] <- matrix(
          0,
          ncol = dim(output[[i]][[j]][[1]])[2],
          nrow = dim(output[[i]][[j]][[1]])[1]
        )
        correlation[[j]] <- matrix(
          0,
          ncol = dim(output[[i]][[j]][[1]])[2],
          nrow = dim(output[[i]][[j]][[1]])[1]
        )
        colnames(correlation[[j]]) <- colnames(output[[i]][[j]][[1]])
        rownames(correlation[[j]]) <- colnames(output[[i]][[j]][[1]])
        for (m in 1:dim(output[[i]][[j]][[1]])[1]) {
          for (n in 1:dim(output[[i]][[j]][[1]])[2]) {
            covariance[[j]][m, n] <- (1 / (output[[i]][[j]][[3]][m, n] - 1)) *
              (output[[i]][[j]][[1]][m, n]) -
              (1 /
                (output[[i]][[j]][[3]][m, n] *
                  (output[[i]][[j]][[3]][m, n] - 1))) *
                output[[i]][[j]][[2]][m, n] *
                output[[i]][[j]][[2]][n, m]
          }
        }
        correlation[[j]] <- stats::cov2cor(covariance[[j]])
        splitresults[[j]] <- list(
          output[[i]][[j]][[4]][[1]],
          output[[i]][[j]][[4]][[2]],
          correlation[[j]],
          output[[i]][[j]][[3]]
        )
        n1 <- "Number of missing values in each variable"
        n2 <- "Number of missing values casewise"
        n3 <- "Correlation Matrix"
        n4 <- "Number of complete cases used"
        names(splitresults[[j]]) <- c(n1, n2, n3, n4)
      }
      names(splitresults) <- names(output[[i]])
      results[[i]] <- splitresults
    }
    names(results) <- stdnames

    results_data <- do.call(
      rbind,
      lapply(names(results), function(server_name) {
        server_results <- results[[server_name]]

        correlation_data <- data.frame(
          lsoa11cd = names(server_results),

          complete.counts = vapply(
            server_results,
            function(x) {
              x[["Number of complete cases used"]][1, 2]
            },
            numeric(1)
          ),

          correlation.coef = vapply(
            server_results,
            function(x) {
              x[["Correlation Matrix"]][1, 2]
            },
            numeric(1)
          ),

          server = server_name,
          row.names = NULL
        )

        if (!draw.plot) {
          dplyr::left_join(
            shape_list[[server_name]],
            correlation_data,
            by = "lsoa11cd"
          )
        } else {
          correlation_data
        }
      })
    )

    if (draw.plot) {
      results_data <- sf::st_as_sf(results_data)
      rownames(results_data) <- NULL
    } else {
      rownames(results_data) <- NULL
      return(results_data)
    }
  } else if (type == "combine") {
    common_lsoas <- sort(Reduce(intersect, lapply(output, names)))

    if (
      is.null(common_lsoas) ||
        length(common_lsoas) < nfilter.tab
    ) {
      stop("Common lsoa in servers is less than the safe limit")
    }
    # lsoa in all servers should match or warning that common lsoa is only analysed
    # in that case check if the matrix size is not disclosive
    combined.results <- list()
    for (lsoa in common_lsoas) {
      reference <- output[[1]][[lsoa]]

      combined.sums.of.products <- matrix(
        0,
        nrow = nrow(reference[[1]]),
        ncol = ncol(reference[[1]]),
        dimnames = dimnames(reference[[1]])
      )

      combined.sums <- matrix(
        0,
        nrow = nrow(reference[[2]]),
        ncol = ncol(reference[[2]]),
        dimnames = dimnames(reference[[2]])
      )

      combined.complete.cases <- matrix(
        0,
        nrow = nrow(reference[[3]]),
        ncol = ncol(reference[[3]]),
        dimnames = dimnames(reference[[3]])
      )

      combined.missing.cases.vector <- matrix(
        0,
        nrow = nrow(reference[[4]][[1]]),
        ncol = ncol(reference[[4]][[1]]),
        dimnames = dimnames(reference[[4]][[1]])
      )

      combined.missing.cases.matrix <- matrix(
        0,
        nrow = nrow(reference[[4]][[2]]),
        ncol = ncol(reference[[4]][[2]]),
        dimnames = dimnames(reference[[4]][[2]])
      )

      combined.sums.of.squares <- matrix(
        0,
        nrow = nrow(reference[[5]]),
        ncol = ncol(reference[[5]]),
        dimnames = dimnames(reference[[5]])
      )

      # Combine sufficient statistics across relevant servers
      for (i in stdnames) {
        lsoa_output <- output[[i]][[lsoa]]

        combined.sums.of.products <-
          combined.sums.of.products + lsoa_output[[1]]

        combined.sums <-
          combined.sums + lsoa_output[[2]]

        combined.complete.cases <-
          combined.complete.cases + lsoa_output[[3]]

        combined.missing.cases.vector <-
          combined.missing.cases.vector + lsoa_output[[4]][[1]]

        combined.missing.cases.matrix <-
          combined.missing.cases.matrix + lsoa_output[[4]][[2]]

        combined.sums.of.squares <-
          combined.sums.of.squares + lsoa_output[[5]]
      }

      # Calculate covariance matrix for this LSOA
      combined.covariance <- matrix(
        NA_real_,
        nrow = nrow(combined.sums.of.products),
        ncol = ncol(combined.sums.of.products),
        dimnames = dimnames(combined.sums.of.products)
      )

      for (m in seq_len(nrow(combined.sums.of.products))) {
        for (n in seq_len(ncol(combined.sums.of.products))) {
          n_complete <- combined.complete.cases[m, n]

          if (!is.na(n_complete) && n_complete > 1) {
            combined.covariance[m, n] <-
              combined.sums.of.products[m, n] /
              (n_complete - 1) -
              (combined.sums[m, n] *
                combined.sums[n, m]) /
                (n_complete * (n_complete - 1))
          }
        }
      }

      # Convert covariance to correlation
      combined.correlation <- tryCatch(
        stats::cov2cor(combined.covariance),
        error = function(e) {
          matrix(
            NA_real_,
            nrow = nrow(combined.covariance),
            ncol = ncol(combined.covariance),
            dimnames = dimnames(combined.covariance)
          )
        }
      )

      combined.results[[lsoa]] <- list(
        "Number of missing values in each variable" = combined.missing.cases.vector,

        "Number of missing values casewise" = combined.missing.cases.matrix,

        "Number of complete cases used" = combined.complete.cases,

        "Covariance Matrix" = combined.covariance,

        "Correlation Matrix" = combined.correlation
      )
    }

    results <- combined.results

    correlation_data <- data.frame(
      lsoa11cd = names(results),
      complete.counts = sapply(results, function(x) {
        x$`Number of complete cases used`[1, 2]
      }),
      correlation.coef = sapply(results, function(x) {
        x$`Correlation Matrix`[1, 2]
      }),
      server = if (type == "split") {
        names(results)
      } else {
        rep("combine", length(results))
      },
      row.names = NULL
    )
    if (draw.plot) {
      results_data <- dplyr::left_join(
        shape_list[[1]],
        correlation_data,
        by = "lsoa11cd"
      ) |>
        dplyr::mutate(server = "combine")
    } else {
      results_data <- correlation_data
      return(results_data)
    }
  } else {
    stop(
      'Function argument "type" has to be either "combine" or "split"',
      call. = FALSE
    )
  }

  # Plotting
  ##### Load shapefile
  corplot <- ggplot2::ggplot(results_data) +
    ggplot2::geom_sf(ggplot2::aes(fill = correlation.coef)) +
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
    ggplot2::labs(fill = "Correlation Coefficient") +
    ggplot2::facet_wrap(~server) +
    ggplot2::theme_minimal()

  return(corplot)
}
