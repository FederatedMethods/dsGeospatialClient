#'
#' @title  Computes the mean and standard deviation across groups defined by one factor
#' @description This function calculates the mean and SD of a continuous variable for each class of
#' a single factor.
#' @details This function calculates the mean, standard deviation (SD), N (number of observations)
#' and the standard error of the mean (SEM) of a continuous variable broken down into subgroups
#' defined by a single factor. 
#' 
#' There are important differences between \code{ds.meanSdGp} function compared to 
#' the function \code{ds.meanByClass}:
#' 
#' (A) \code{ds.meanSdGp} does not actually subset the data it simply calculates the required statistics
#' and reports them. This means you cannot use this function if you wish to physically break the
#' data into subsets. On the other hand, it makes the function very much faster than \code{ds.meanByClass}
#' if you do not need to create physical subsets. \cr
#' (B) \code{ds.meanByClass} allows you to specify up to
#' three categorising factors, but \code{ds.meanSdGp} only allows one. However, this is not a serious
#' problem. If you have two factors (e.g. sex with two levels \code{[0,1]} and \code{BMI.categorical} with
#' three levels \code{[1,2,3]}) you simply need to create a new factor that combines the two together in a
#' way that gives each combination of levels a different value in the new factor. So, in the
#' example given, the calculation \code{newfactor = (3*sex) + BMI} gives you six values: \cr
#' (1) \code{sex = 0} and \code{BMI = 1} -> \code{newfactor = 1} \cr
#' (2) \code{sex = 0} and \code{BMI = 2} -> \code{newfactor = 2} \cr
#' (3) \code{sex = 0} and \code{BMI = 3} -> \code{newfactor = 3} \cr
#' (4) \code{sex = 1} and \code{BMI = 1} -> \code{newfactor = 4} \cr
#' (5) \code{sex = 1} and \code{BMI = 2} -> \code{newfactor = 5} \cr
#' (6) \code{sex = 1} and \code{BMI = 3} -> \code{newfactor = 6} \cr
#' 
#' (C) At present, \code{ds.meanByClass} calculates the sample size in each group to mean the 
#' total sample size (i.e. it
#' includes all observations in each group regardless of whether or not they include missing values
#' for the continuous variable or the factor). The calculation of sample size in each group by
#' \code{ds.meanSdGp} always reports the number of observations that are non-missing both for the
#' continuous variable and the factor. This makes sense - in the case of \code{ds.meanByClass},
#' the total size of the physical subsets was important, 
#' but when it comes down only to \code{ds.meanSdGp} which
#' undertakes analysis without physical subsetting,  it is only the observations with non-missing
#' values in both variables that contribute to the calculation of means and SDs within each group
#' and so it is logical to consider those counts as primary. The only reference \code{ds.meanSdGp} makes
#' to missing counts is in the reporting of \code{Ntotal} and \code{Nmissing} overall (ie not broken down by
#' group). 
#' 
#' For the future, we plan to extend \code{ds.meanByClass} to report both total and non-missing
#' counts in subgroups.
#' 
#' Depending on the variable \code{type} can be carried out different analysis:\cr
#' (1) \code{"combine"}: a pooled table of results is generated. \cr
#' (2) \code{"split"} a table of results is generated for each study. \cr
#' 
#' Server function called: \code{meanSdGpDS}
#' @param x a character string specifying the name of a numeric continuous
#' variable.
#' @param y  a character string specifying the name of a categorical
#' variable of class factor.
#' @param type a character string that represents the type of analysis to carry out.
#' This can be set as: \code{"combine"} or \code{"split"}. 
#' Default \code{"combine"}. 
#' For more information see \strong{Details}.  
#' @param do.checks logical. If TRUE the administrative checks
#' are undertaken to ensure that the input objects are defined in all studies and that the
#' variables are of equivalent class in each study. 
#' Default is FALSE to save time.  
#' @param datasources a list of \code{\link[DSI]{DSConnection-class}} 
#' objects obtained after login. If the \code{datasources} argument is not specified
#' the default set of connections will be used: see \code{\link[DSI]{datashield.connections_default}}.
#' @return \code{ds.meanSdGp} returns to the client-side the mean, SD, Nvalid and SEM combined
#' across studies and/or separately for each study, depending on the argument \code{type}. 
#' 
#' @author DataSHIELD Development Team
#' @seealso \code{\link{ds.subsetByClass}} to subset by the classes of factor vector(s).
#' @seealso \code{\link{ds.subset}} to subset by complete cases (i.e. removing missing values), threshold,
#' columns and rows.
#' @export
#' @examples
#' \dontrun{
#' 
#'  ## Version 6, for version 5 see the Wiki
#'   
#'   # connecting to the Opal servers
#' 
#'   require('DSI')
#'   require('DSOpal')
#'   require('dsBaseClient')
#'
#'   builder <- DSI::newDSLoginBuilder()
#'   builder$append(server = "study1", 
#'                  url = "http://192.168.56.100:8080/", 
#'                  user = "administrator", password = "datashield_test&", 
#'                  table = "SURVIVAL.EXPAND_NO_MISSING1", driver = "OpalDriver")
#'   builder$append(server = "study2", 
#'                  url = "http://192.168.56.100:8080/", 
#'                  user = "administrator", password = "datashield_test&", 
#'                  table = "SURVIVAL.EXPAND_NO_MISSING2", driver = "OpalDriver")
#'   builder$append(server = "study3",
#'                  url = "http://192.168.56.100:8080/", 
#'                  user = "administrator", password = "datashield_test&", 
#'                  table = "SURVIVAL.EXPAND_NO_MISSING3", driver = "OpalDriver")
#'   logindata <- builder$build()
#'   
#'   connections <- DSI::datashield.login(logins = logindata, assign = TRUE, symbol = "D") 
#'
#'
#'   #Example 1: Calculate the mean, SD, Nvalid and SEM of the continuous variable age.60 (age in
#'   #years centralised at 60), broken down by time.id (a six level factor relating to survival time)
#'   #and report the pooled results combined across studies.
#'  
#'   ds.meanSdGp(x = "D$age.60",
#'               y = "D$time.id",
#'               type = "combine",
#'               do.checks = FALSE,
#'               datasources = connections)
#'               
#'   #Example 2: Calculate the mean, SD, Nvalid and SEM of the continuous variable age.60 (age in
#'   #years centralised at 60), broken down by time.id (a six level factor relating to survival time)
#'   #and report both study-specific results and the pooled results combined across studies.
#'   #Save the returned output to msg.b.
#'   
#'   ds.meanSdGp(x = "D$age.60",
#'               y = "D$time.id",
#'               type = "both",
#'               do.checks = FALSE,
#'               datasources = connections)  
#'              
#'   # clear the Datashield R sessions and logout
#'   datashield.logout(connections)
#' }
#'
ds.summarystat <- function(x=NULL, y=NULL, type='combine', do.checks=FALSE, 
                           draw.plot = FALSE, names_region = NULL, 
                           metric = NULL, datasources=NULL){
  
  # look for DS connections
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }
  
  # ensure datasources is a list of DSConnection-class
  if(!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {methods::is(d,"DSConnection")}))))){
    stop("The 'datasources' were expected to be a list of DSConnection-class objects", call.=FALSE)
  }
  
  if(is.null(x)){
    stop("Please provide the name of the input vector!", call.=FALSE)
  }
  
  if(is.null(y)){
    stop("Please provide the name of the input vector!", call.=FALSE)
  }
  
  if(do.checks){
    
    # check if the input objects are defined in all the studies
    isDefined(datasources, x)
    isDefined(datasources, y)
    
    # call the internal function that checks the input object is of the same class in all studies.
    typ1 <- checkClass(datasources, x)
    typ2 <- checkClass(datasources, y)
  }
  
  if(!is.logical(draw.plot)){
    stop("draw.plot should be logical")
  }
  # names of the studies
  stdnames <- names(datasources)
  
  # variable names

  xvarname <- strsplit(x, "\\$", perl=TRUE)[[1]][2]
  yvarname <- strsplit(y, "\\$", perl=TRUE)[[1]][2]
  
  # call the server side function that calculates mean and standard deviation
  # by group in each study
  calltext <- paste0("MeanSdGpDS(", x, ",", y, ")")
  output <- DSI::datashield.aggregate(datasources, as.symbol(calltext))
  numsources <- length(output)

  lsoa_names <- lapply(output, function(x) names(x$Mean_gp))
  
  all.tables.valid<-1
  for(i in 1:numsources){
    if(unlist(output[[i]]$Table_valid==FALSE)) all.tables.valid <- all.tables.valid-1
  }
  
  if(all.tables.valid<1){
    warning.message<-"At least 1 cell count is 1-nfilter, please regroup"
    return(warning.message)
  }
  
## Plot matrix

  if(!type %in% c("split", "combine"))
  {
    warning("Invalide 'type' argument. Defaulting to type 'combine'")
    type <- "combine"
  }

  plot_matrix <- list()
  split <- type == 'split'
  
  ## names_region check 
  # check if list and all elements are character 

  if (is.list(names_region) &&
      length(names_region) > 0 &&
      all(is.character(unlist(names_region)))) {
    names_region.valid = TRUE
    
    
    names_region <- lapply(names_region, function(i){
      if("CheshireMercyside" %in% i){
        indx <- which(unlist(i) == "CheshireMercyside")
        i <- i[-indx]
        i <- c(i, "Cheshire East","Cheshire West and Chester",
               "Halton","Knowsley","Liverpool","Sefton", 
               "St. Helens","Warrington","Wirral")
      }else{
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
      } else{
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
        } else{
          names_region <- lapply(names_region, unique)
        }
      }
      } else if(!split) {
      # when combine, list must have length 1
      # ifnot repeat list[[1]] for all numsources
      if (length(names_region) != 1) {
        warning(
          "length of names_region must be 1 for type 'combine',
                 using only first elements in names_region to return"
        )
        names_region = list(names_region[[1]])
        names_region <- lapply(names_region, unique)
        
      } else{
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
    q.names <- c("5th_Quantile", "10th_Quantile", "25th_Quantile", "50th_Quantile",
                 "75th_Quantile", "90th_Quantile", "95th_Quantile")
    
    if(is.null(metric) || !metric %in% c("Mean", "SD","SEM", q.names)){
      warning("Invalid metric type, returning available metric table")
      draw.plot = FALSE
      names_region.valid =FALSE
    } else {
      
      if(!split){
      metric_map <- list("Mean" = "mean.gp", 
                         "SD" = "SD.gp", 
                         "SEM" = "SEM.gp")
      
      } else {
      metric_map <- list("Mean" = "mean.gp.study", 
                         "SD" = "SD.gp.study", 
                         "SEM" = "SEM.gp.study")
      }
    }
    
  }else{
    warning("invalid region names/format, returning the whole table")
    draw.plot = FALSE
    names_region.valid = FALSE
  }
 

##########
  
  if(all.tables.valid==1){
    # COMBINE 
    if(type!="split"){
      numsources <- length(output)
      mean.matrix <- NULL
      sd.matrix <- NULL
      qq.matrix <- NULL
      n.matrix <- NULL
      Nvalid <- 0
      Nmissing <- 0
      Ntotal <- 0
      
      for(j in 1:numsources){
        mean.matrix <- rbind(mean.matrix,as.numeric(unlist(output[[j]][2])))
        sd.matrix <- rbind(sd.matrix,as.numeric(unlist(output[[j]][3])))
        qq.matrix <- rbind(qq.matrix,as.numeric(unlist(output[[j]][4])))
        n.matrix <- rbind(n.matrix,as.numeric(unlist(output[[j]][5])))
        Nvalid <- Nvalid+as.numeric(unlist(output[[j]][6]))
        Nmissing <- Nmissing+as.numeric(unlist(output[[j]][7]))
        Ntotal <- Ntotal+as.numeric(unlist(output[[j]][8]))
        

      }
     
      var.matrix <- sd.matrix^2
      nsum.vector <- rep(1,numsources)
      # Calculate weighted means across studies in each group
      mean.gp <- (diag(t(mean.matrix)%*%n.matrix))/(t(n.matrix)%*%nsum.vector)
      
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
      
      colnames(qq.gp) <- c("5%_gp", "10%_gp", "25%_gp", "50%_gp", "75%_gp", "90%_gp", "95%_gp")
      
      qq.gp <- as.data.frame(qq.gp)
      
      # Calculate weighted SDs across studies in each group
      var.gp <- (diag(t(var.matrix)%*%n.matrix))/(t(n.matrix)%*%nsum.vector)
      SD.gp <- sqrt(var.gp)
      N.gp <- (t(n.matrix) %*% nsum.vector)
      SEM.gp <- SD.gp/sqrt(N.gp)
      
      # create names
      names.gp <- rep(NA,length(mean.gp))
      
      for(k in 1:length(mean.gp)){
        names.gp[k] <- paste0(yvarname,"_",k)
      }
      
      dimnames(mean.gp) <- c(list(names.gp),list("Mean_gp"))
      dimnames(SD.gp) <- c(list(names.gp),list("SD_gp"))
      dimnames(N.gp) <- c(list(names.gp),list("Nvalid_gp"))
      dimnames(SEM.gp) <- c(list(names.gp),list("SEM_gp"))
      
    }
    
    # SPLIT 
    if(type!="combine"){
      numsources <- length(output)
      mean.matrix <- NULL
      sd.matrix <- NULL
      qq.matrix <- list()
      n.matrix <- NULL
      Nvalid <- 0
      Nmissing <- 0
      Ntotal <- 0
      
      for(j in 1:numsources){
        mean.matrix <- rbind(mean.matrix,as.numeric(unlist(output[[j]][2])))
        sd.matrix <- rbind(sd.matrix,as.numeric(unlist(output[[j]][3])))
        qq.matrix[[j]] <- rbind(qq.matrix,as.numeric(unlist(output[[j]][4])))
        n.matrix <- rbind(n.matrix,as.numeric(unlist(output[[j]][5])))
        Nvalid <- Nvalid+as.numeric(unlist(output[[j]][6]))
        Nmissing <- Nmissing+as.numeric(unlist(output[[j]][7]))
        Ntotal <- Ntotal+as.numeric(unlist(output[[j]][8]))
        

      }
      var.matrix <- sd.matrix^2
      mean.gp.study <- t(mean.matrix)
      SD.gp.study <- t(sd.matrix)
      N.gp.study <- t(n.matrix)
      SEM.gp.study <- SD.gp.study/sqrt(N.gp.study)
      
      qq.gp.study <- lapply(qq.matrix, function(i){
                              matrix(i, ncol = 7, 
                                     byrow = TRUE)})
      
      qq.gp.study <- lapply(qq.gp.study, function(i) {
        colnames(i) <- c("5%_gp", "10%_gp", "25%_gp", "50%_gp",
                         "75%_gp", "90%_gp", "95%_gp")
        i
      })
      
      
     # create names
      names.gp <- rep(NA,dim(mean.gp.study)[1])
      
      for(k in 1:dim(mean.gp.study)[1]){
        names.gp[k] <- paste0(yvarname,"_",k)
      }
      
      names.study <- names(datasources)
      
      dimnames(mean.gp.study) <- c(list(names.gp),list(names.study))
      dimnames(SD.gp.study) <- c(list(names.gp),list(names.study))
      dimnames(N.gp.study) <- c(list(names.gp),list(names.study))
      dimnames(SEM.gp.study) <- c(list(names.gp),list(names.study))
      
    }
    
    if(type=="combine"){
      
      # if type is combine check whether lsoa levels are same in the two studies
      names.identical <- all(vapply(lsoa_names, identical, logical(1),y = lsoa_names[[1]]))
      
      if (!names.identical) {
        warning(
          "LSOA levels across studies are not the same. ",
          "Cannot proceed with type = 'combine'. ",
          "Defaulting to type = 'split'."
        )
        type <- "split"
        
      } else {
        
      lsoanames <- lsoa_names[[1]]
      if(draw.plot){
        
      ## If plot is TRUE please select which metric to plot else table is returned
      
      if(metric %in% q.names){
        
        sel.metric.data <- data.frame(lsoa11cd = rownames(qq.gp))
        sel.metric.data$value <- qq.gp[ , which(metric == q.names)]
        sel.metric.data$server <- 'combine'
       
      } else {
        sel.metric <- get(metric_map[[metric]])
        
        sel.metric.data <- data.frame(lsoa11cd = lsoanames,
                              value = as.numeric(sel.metric),
                              server = 'combine')
      }
      shape_sf <- shape_list[[1]] |>
        dplyr::left_join(sel.metric.data, by = "lsoa11cd") |>
        sf::st_as_sf()
      

      } else {
        result <- list(mean.gp,SD.gp,N.gp,SEM.gp,qq.gp,Nvalid,Nmissing,Ntotal, lsoa_names[[1]])
        names(result) <- list("Mean_gp","StDev_gp","Nvalid_gp","SEM_gp", "Q_gp","Total_Nvalid",
                              "Total_Nmissing","Total_Ntotal", "LSOAnames")
        return(result)
      }
     }
    }
    
    if(type=="split"){
    
    if(draw.plot){
      browser()
      if(metric %in% q.names){
        
        sel.metric.data <- 
          lapply(1:length(qq.gp.study), function(i) {
            data.frame(
              server = names(datasources)[i],
              value = qq.gp.study[[i]][ ,which(metric == q.names)],
              lsoa11cd = lsoa_names[[i]],
              row.names = NULL
            )
          })
        
      }else{
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
      lapply(1:numsources, function(i){
      shape_list[[i]] |>
      dplyr::left_join(sel.metric.data[[i]], by = "lsoa11cd") |>
      sf::st_as_sf()
    })
    )
    

    }  else {
      
      lsoanames <- matrix(unlist(lsoa_names), ncol = numsources)
      dimnames(lsoanames) <- dimnames(mean.gp.study) 
      result <- list(mean.gp.study,SD.gp.study,N.gp.study,SEM.gp.study,qq.gp.study,Nvalid,
                     Nmissing,Ntotal, lsoanames)
      names(result) <- list("Mean_gp_study","StDev_gp_study","Nvalid_gp_study",
                            "SEM_gp_study","Q_gp","Total_Nvalid","Total_Nmissing",
                            "Total_Ntotal", "LSOAnames")
      return(result)
    }
    }
    
  }
  
  if(draw.plot){
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