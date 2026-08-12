#' Calculate chi-squared statistics for grouped contingency tables
#'
#' @description
#' Calculates Pearson's chi-squared tests for two categorical variables
#' separately across the levels of a third grouping variable. The analysis
#' can be performed separately for each DataSHIELD datasource or after
#' combining contingency-table counts across datasources.
#'
#' Optionally, when the grouping variable contains 2011 LSOA codes, the
#' chi-squared statistics can be joined to LSOA boundary data and returned
#' as spatial data or plotted as maps.
#'
#' @details
#' `ds.chisqr()` constructs three-dimensional contingency tables from
#' `rvar`, `cvar`, and `stvar`. The `rvar` argument defines the rows of
#' each contingency table, `cvar` defines the columns, and `stvar`
#' identifies the groups for which separate two-dimensional contingency
#' tables and chi-squared tests are produced.
#'
#' For example, if `rvar` represents sex, `cvar` represents asthma status,
#' and `stvar` represents LSOA, a separate contingency table and Pearson's
#' chi-squared test are calculated for each LSOA.
#'
#' Variables supplied to `rvar`, `cvar`, and `stvar` may be numeric,
#' integer, or factor variables. Their levels are identified on the
#' server side before constructing tables with a consistent structure
#' across the participating datasources.
#'
#' When `split = TRUE`, chi-squared tests are calculated independently
#' for each datasource. When `split = FALSE`, cell counts are first
#' combined across the valid datasources and the chi-squared tests are
#' calculated from the combined contingency tables.
#'
#' The returned chi-squared statistics are calculated using
#' [stats::chisq.test()]. For every level of `stvar`, the output includes
#' the chi-squared statistic, degrees of freedom, p-value, and datasource.
#'
#' If `names_region` is supplied, LSOA boundaries are obtained using
#' `boundr::bounds()` and matched against the levels of `stvar`. Only
#' LSOAs occurring within the requested region or regions are included.
#'
#' When `draw.plot = TRUE` and valid region names have been supplied,
#' the chi-squared statistics are displayed geographically using
#' [ggplot2::geom_sf()]. Separate maps are faceted by datasource when
#' `split = TRUE`.
#'
#' As with the underlying DataSHIELD table functionality, contingency
#' tables remain subject to the server-side disclosure-control settings.
#' A datasource that fails the relevant table checks is excluded from
#' the analysis and reported in the validity information.
#'
#' @param rvar Character string specifying the name of the variable defining
#'   the rows of each two-dimensional contingency table.
#'
#' @param cvar Character string specifying the name of the variable defining
#'   the columns of each two-dimensional contingency table.
#'
#' @param stvar Character string specifying the grouping variable defining
#'   the separate two-dimensional contingency tables. For spatial output,
#'   this should contain 2011 LSOA codes corresponding to the selected
#'   boundary data.
#'
#' @param exclude Levels to exclude when constructing the contingency
#'   tables. This argument is passed to the server-side table operation.
#'   See [base::table()] for details of the interaction between `exclude`
#'   and `useNA`.
#'
#' @param useNA Character string specifying whether missing values should
#'   be included in the contingency tables. Accepted values are `"no"`
#'   and `"always"`. Default is `"always"`.
#'
#' @param suppress.chisq.warnings Logical. If `TRUE`, warnings generated
#'   by [stats::chisq.test()], for example warnings concerning small
#'   expected cell counts, are suppressed. Default is `FALSE`.
#'
#' @param datasources A list of `DSConnection-class` objects obtained after
#'   login. If `NULL`, the function searches for available DataSHIELD
#'   connections using `datashield.connections_find()`.
#'
#' @param force.nfilter Optional character string specifying an integer
#'   value used to increase the server-side `nfilter.tab` disclosure
#'   threshold. The value can only increase, and not decrease, the
#'   disclosure-control threshold configured by the data custodian.
#'
#' @param names_region Optional character vector containing the names of
#'   regions used to obtain 2011 LSOA boundaries. When `split = TRUE`,
#'   region names correspond to the participating datasources. When
#'   `split = FALSE`, only one region is used.
#'
#' @param split Logical. If `TRUE`, chi-squared tests are calculated
#'   separately for each datasource. If `FALSE`, counts are combined
#'   across datasources before calculating the chi-squared tests.
#'   Default is `TRUE`.
#'
#' @param draw.plot Logical. If `TRUE` and valid `names_region` have been
#'   supplied, returns a faceted spatial plot of the chi-squared
#'   statistics. If `FALSE` and valid regions are supplied, the spatially
#'   matched results are returned as a data frame without geometry.
#'   Default is `TRUE`.
#'
#' @return
#' The return value depends on `split`, `names_region`, and `draw.plot`.
#'
#' When no valid region information is supplied, a data frame is returned
#' containing:
#'
#' \itemize{
#'   \item `lsoa11cd` - level of the grouping variable (`stvar`);
#'   \item `X.squared` - Pearson's chi-squared statistic;
#'   \item `df` - degrees of freedom;
#'   \item `p.value` - p-value of the chi-squared test;
#'   \item `datasource` - datasource name, or `"COMBINED"` when
#'     `split = FALSE`.
#' }
#'
#' If valid `names_region` are supplied and `draw.plot = FALSE`, the same
#' results are returned after matching them to the selected LSOA
#' boundaries.
#'
#' If valid `names_region` are supplied and `draw.plot = TRUE`, a
#' `ggplot` object containing maps of `X.squared` is returned, faceted
#' by datasource.
#'
#' If all participating studies fail the server-side table checks, a list
#' containing `validity.message` and `error.messages` is returned.
#'
#' @examples
#' \dontrun{
#' # Calculate chi-squared statistics separately for each datasource
#' ds.chisqr(
#'   rvar = "D$sex",
#'   cvar = "D$has_asthma",
#'   stvar = "D$lsoa11cd",
#'   split = TRUE,
#'   draw.plot = FALSE,
#'   datasources = conns
#' )
#'
#' # Calculate combined chi-squared statistics
#' ds.chisqr(
#'   rvar = "D$sex",
#'   cvar = "D$has_asthma",
#'   stvar = "D$lsoa11cd",
#'   split = FALSE,
#'   draw.plot = FALSE,
#'   datasources = conns
#' )
#'
#' # Map chi-squared statistics by LSOA
#' ds.chisqr(
#'   rvar = "D$sex",
#'   cvar = "D$has_asthma",
#'   stvar = "D$lsoa11cd",
#'   names_region = "Liverpool",
#'   split = FALSE,
#'   draw.plot = TRUE,
#'   datasources = conns
#' )
#' }
#'
#' @author DataSHIELD Development Team
#'
#' @export
#' 
ds.chisqr <- function(rvar=NULL, cvar=NULL, stvar=NULL, exclude=NULL,	
                      useNA ="always", suppress.chisq.warnings=FALSE, 
                      datasources=NULL, force.nfilter=NULL, names_region = NULL, 
                      split = TRUE, draw.plot = TRUE){
  
  # if no connection login details are provided look for 'connection' objects in the environment
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }

  # ensure datasources is a list of DSConnection-class
  if(!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {methods::is(d,"DSConnection")}))))){
    stop("The 'datasources' were expected to be a list of DSConnection-class objects", call.=FALSE)
  }
  
  # check if a value has been provided for rvar
  if(is.null(rvar)){
    return("Error: rvar must have a value which is a character string naming the row variable for the table")
  }
  
  # check if the input object is defined in all the studies
  isDefined(datasources, rvar)
  
  if(is.null(cvar) || !is.character(cvar)){
    return("Error: cvar must have a value which is a character string naming the column variable for the table")
  }
  
  if(!is.null(cvar)){
    isDefined(datasources, cvar)
    cvar.transmit<-cvar
  }
  
  if(is.null(stvar) || !is.character(stvar)){
    return("Error: if stvar must have a value which is a character string naming the variable coding separate tables for the table")
  }
  
  if(!is.null(stvar)){
    isDefined(datasources, stvar)
    stvar.transmit<-stvar
  }
  
  if(useNA!="no" && useNA!="always"){
    stop("useNA must be either 'no' or 'always'.")
  }
  
  if(!is.null(force.nfilter)&&!is.character(force.nfilter)){
    return("Error: if force.nfilter is not null, it must have a value which is a character string specifying an integer for the forced value of the nfilter")
  }
  
  
  #All arguments should be directly transmittable
  rvar.transmit<-rvar

  if(is.null(exclude))
  {
    exclude.transmit<-NULL
  }
  else
  {
    exclude.transmit<-paste0(as.character(exclude),collapse=",")
  }
  
  useNA.transmit<-useNA
  
  if(is.null(force.nfilter))
  {
    force.nfilter.transmit<-NULL
  }
  else
  {
    force.nfilter.transmit<-force.nfilter
  }
  #CALL THE asFactorDS3 SERVER SIDE FUNCTION (AN AGGREGATE FUNCTION)
  # FOR rvar, cvar AND stvar  
  #TO DETERMINE ALL OF THE LEVELS REQUIRED
  
  rvar.asfactor.calltext <- call("asFactorDS3", rvar)
  rvar.all.levels <- DSI::datashield.aggregate(datasources, rvar.asfactor.calltext)
  
  numstudies <- length(datasources)
  
  rvar.all.levels.all.studies <- NULL
  
  for(j in 1:numstudies){
    rvar.all.levels.all.studies <- c(rvar.all.levels.all.studies,rvar.all.levels[[j]])
  }
  
  if (length(rvar.all.levels.all.studies) == 0) {
    stop(paste0("Unable to obtain factors for rvar: '", rvar , "'"), call. = FALSE)
  }
  
  rvar.all.unique.levels <- as.character(unique(rvar.all.levels.all.studies))
  
  rvar.all.unique.levels.transmit <- paste0(rvar.all.unique.levels, collapse=",")
  
  ########################################################
  
  
    cvar.asfactor.calltext <- call("asFactorDS3", cvar)
    cvar.all.levels <- DSI::datashield.aggregate(datasources, cvar.asfactor.calltext)
    
    numstudies <- length(datasources)
    
    cvar.all.levels.all.studies <- NULL
    
    for(j in 1:numstudies){
      cvar.all.levels.all.studies <- c(cvar.all.levels.all.studies,cvar.all.levels[[j]])
    }
    
    if (length(cvar.all.levels.all.studies) == 0) {
      stop(paste0("Unable to obtain factors for cvar: '", cvar , "'"), call. = FALSE)
    }
    
    cvar.all.unique.levels <- as.character(unique(cvar.all.levels.all.studies))
    
    cvar.all.unique.levels.transmit <- paste0(cvar.all.unique.levels, collapse=",")
  
  ########################################################
    plot_matrix <- list()
    
    if(!is.null(names_region) && all(is.character(names_region))){
      names_region.valid = TRUE
      
      if(!split && length(names_region) > 1){
        warning(paste0("more than one region selected for split = ", split,
                       ", using only first names_region to return"))
        names_region = names_region[1]
      }
      
      if(split && length(names_region) != numstudies){
        names_region <- rep(names_region, times = numstudies)
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
      
    }else{
      warning("invalid region names, returning the whole table")
      draw.plot = FALSE
      names_region.valid = FALSE
    }

    stvar.asfactor.calltext <- call("asFactorDS3", stvar)
    stvar.all.levels <- DSI::datashield.aggregate(datasources, stvar.asfactor.calltext)
  
    if (exists("shape_list")){
    matched_lsoas <- Map(
      function(stvar_levels, shape) {
        stvar_levels[stvar_levels %in% shape$lsoa11cd]
      },
      stvar.all.levels,
      shape_list
    )
    stvar.all.levels <- matched_lsoas
    }
    
    numstudies <- length(datasources)
    
    stvar.all.levels.all.studies <- NULL
    
    for(j in 1:numstudies){
      stvar.all.levels.all.studies <- c(stvar.all.levels.all.studies,stvar.all.levels[[j]])
    }
    
    if (length(stvar.all.levels.all.studies) == 0) {
      stop(paste0("Unable to obtain factors for stvar: '", stvar , "'"), call. = FALSE)
    }
    
    stvar.all.unique.levels <- as.character(unique(stvar.all.levels.all.studies))
    
    stvar.all.unique.levels.transmit <- paste0(stvar.all.unique.levels, collapse=",")
  
  # CALL THE MAIN SERVER SIDE AGGREGATE FUNCTION
  
  calltext <- call("tableDS", rvar.transmit=rvar.transmit, cvar.transmit=cvar.transmit,
                   stvar.transmit=stvar.transmit,rvar.all.unique.levels.transmit=rvar.all.unique.levels.transmit,
                   cvar.all.unique.levels.transmit=cvar.all.unique.levels.transmit,
                   stvar.all.unique.levels.transmit=stvar.all.unique.levels.transmit,
                   exclude.transmit=exclude.transmit, useNA.transmit=useNA.transmit,
                   force.nfilter.transmit=force.nfilter.transmit)
  
  
  table.out<-DSI::datashield.aggregate(datasources, calltext)
  
  #END OF MAIN FUNCTION LEADING UP TO CALL
  ##############################################################################
  ## plotting 
  #Take serverside output and set up arrays with the correct dimensions
  #so all values of a given dimension in any study are included in the
  #tables produced for each individual study and for the combined values
  #across studies 
  
  #How many studies have returned tables
  numsources.orig<-length(table.out)
  table.out.orig<-table.out
  
  
  #Check whether return is a table or an error message
  valid.output <- rep(1, numsources.orig)
  error.messages <- table.out
  sum.valid <- 0
  study.names.valid <- NULL
  list.temp <- NULL
  
  
  for(ns in 1:numsources.orig)
  {
    if("character" %in% class(table.out[[ns]]))
    {
      valid.output[ns]<-0
      error.messages[[ns]]<-table.out[[ns]]
    }	
    else
    {
      error.messages[[ns]]<-"No errors reported from this study"
    }
  }
  
  num.valid.studies <- sum(valid.output)
  
  
  if(num.valid.studies==0)
  {
      validity.message<-"All studies failed for reasons identified below"
      message("\n",validity.message,"\n\n")
      for(ns in 1:numsources.orig)
      {
        message("\nStudy",ns,": ",error.messages[[ns]],"\n")
      }
      
      return(list(validity.message=validity.message,error.messages=error.messages))
    
  }
  
  
  for(ns in 1:numsources.orig)
  {
    if(valid.output[ns])
    {
      sum.valid<-sum.valid+1
      list.temp<-paste0(list.temp,"table.out.orig[[",ns,"]]")
      if(sum.valid < num.valid.studies)
      {
        list.temp<-paste0(list.temp,",")
      }
      study.names.valid<-c(study.names.valid,names(datasources)[ns])
    }
  }
  
  
  table.out.valid <- FALSE
  list.text<-paste0("table.out.valid<-list(",list.temp,")")
  
  eval(parse(text=list.text))
  
  table.out<-table.out.valid
  numsources<-length(table.out)
  
  if(num.valid.studies > 0 && num.valid.studies<numsources.orig)
  {
    validity.message<-"At least one study failed for reasons identified by 'error.messages':"
    
    for(ns in 1:numsources.orig)
    {
      message.add<-paste0("Study",ns,": ",error.messages[[ns]])
      validity.message<-c(validity.message,message.add)
    }
    
    #	message("\n",validity.message,"\n")
    #	for(ns in 1:numsources.orig)
    #		{
    #		message("\nStudy",ns,": ",error.messages[[ns]])
    #		}
    #		message("\n\n")
    
    #table.out<-table.out.valid
    #numsources<-length(table.out)
  }
  
  if(num.valid.studies==numsources.orig)
  {
    validity.message<-"Data in all studies were valid"

      message("\n",validity.message,"\n")
      for(ns in 1:numsources.orig)
      {
        message("\nStudy",ns,": ",error.messages[[ns]])
      }
      message("\n\n")

  }
  
  #check all tables from all sources have the same number of dimensions
  
  table.dimensions<-rep(0,numsources)
  for(ns in 1:numsources)
  {
    table.dimensions[ns]<-length(dim(table.out[[ns]]))
  }
  #table.dimensions
  all.dims.same<-TRUE
  if(numsources>1)
  {
    for(ns in 1:(numsources-1))
    {
      if(table.dimensions[ns]!=table.dimensions[ns+1])
      {
        all.dims.same<-FALSE
        return.message<-"Warning: tables in different sources have different numbers of dimensions. Please analyse and combine yourself from study.specific tables above" 
        return(return.message)   
      }
    }
 
  }
  
  ######################################################################
  #Work first with three dimensional tables

    #identify all possible values of each dimension
    
    rvar.dimnames<-NULL
    cvar.dimnames<-NULL
    stvar.dimnames<-NULL
    
    for(ns in 1:numsources)
    {
      rvar.dimnames<-c(rvar.dimnames,dimnames(table.out[[ns]])$rvar)
      rvar.dimnames[is.na(rvar.dimnames)]<-"NA"
      cvar.dimnames<-c(cvar.dimnames,dimnames(table.out[[ns]])$cvar)
      cvar.dimnames[is.na(cvar.dimnames)]<-"NA"
      stvar.dimnames<-c(stvar.dimnames,dimnames(table.out[[ns]])$stvar)
      stvar.dimnames[is.na(stvar.dimnames)]<-"NA"
    }
    
    rvar.dimnames.unique<-unique(rvar.dimnames)
    cvar.dimnames.unique<-unique(cvar.dimnames)
    stvar.dimnames.unique<-unique(stvar.dimnames)
    
    numcells.all.sources<-length(rvar.dimnames.unique)*length(cvar.dimnames.unique)*length(stvar.dimnames.unique)

    empty.table.all.sources.col.1<-rep(rvar.dimnames.unique,times=(length(cvar.dimnames.unique)*length(stvar.dimnames.unique)))
    empty.table.all.sources.col.1
    
    empty.table.all.sources.col.2<-rep(rep(cvar.dimnames.unique,each=length(rvar.dimnames.unique)),times=length(stvar.dimnames.unique))
    empty.table.all.sources.col.2
    
    empty.table.all.sources.col.3<-rep(rep(stvar.dimnames.unique,each=(length(rvar.dimnames.unique)*length(cvar.dimnames.unique))))
    empty.table.all.sources.col.3
    
    empty.table.all.sources.col.4<-rep(0,times=(length(rvar.dimnames.unique)*length(cvar.dimnames.unique)*length(stvar.dimnames.unique)))
    empty.table.all.sources.col.4
    
    empty.table.all.sources<-cbind(empty.table.all.sources.col.1,empty.table.all.sources.col.2,
                                   empty.table.all.sources.col.3,empty.table.all.sources.col.4)
    
    empty.table.all.sources[is.na(empty.table.all.sources)]<-"NA"
   
    dim.vector.all.sources<-c(length(rvar.dimnames.unique),length(cvar.dimnames.unique),
                              length(stvar.dimnames.unique),numsources)
    
    array.all.sources<-array(data=NA,dim=dim.vector.all.sources,
                             dimnames=list(rvar.dimnames.unique,cvar.dimnames.unique,stvar.dimnames.unique,NULL))
    
    names(dimnames(array.all.sources))<-c(rvar,cvar,stvar,"study")
    #print(array.all.sources) #length 108 all NAs so the study specific template for dimnames
    #has been correctly expanded to
    #include all studies
    
    #KEY LOOP
    for(ns in 1:numsources)
    {
      #start with study 1 then 2 etc etc
      
      
      
      numcells<-length(table.out[[ns]])
      #print(numcells)
      study.specific.dim.vect<-dim(table.out[[ns]])
      study.specific.dim.vect
      
      count.in.cell<-rep(NA,numcells)
      stvar.mark<-rep("",numcells)
      cvar.mark<-rep("",numcells)
      rvar.mark<-rep("",numcells)
      cells.so.far<-0
      
      for(ss in 1:study.specific.dim.vect[3])
      {
        for(cc in 1:study.specific.dim.vect[2])
        {
          for(rr in 1:study.specific.dim.vect[1])
          {
            cells.so.far<-cells.so.far+1
            count.in.cell[cells.so.far]<-table.out[[ns]][cells.so.far]
            rvar.mark[cells.so.far]<-rvar.dimnames[rr]
            cvar.mark[cells.so.far]<-cvar.dimnames[cc]
            stvar.mark[cells.so.far]<-stvar.dimnames[ss]
          }
        }
      }
      
      
      table.current.study<-cbind(rvar.mark,cvar.mark,stvar.mark,count.in.cell)
      table.current.study[is.na(table.current.study)]<-"NA"
      #message("current study =",ns)
      #print(table.current.study)
      
      
      array.current.study<-array(data=table.current.study[,4],dim=dim.vector.all.sources[1:3],
                                 dimnames=list(rvar.dimnames.unique,cvar.dimnames.unique,stvar.dimnames.unique))
      names(dimnames(array.current.study))<-c(rvar,cvar,stvar)
      #array.current.study
      
      #IS TABLE FOR CURRENT STUDY IDENTICAL IN STRUCTURE TO EMPTY TABLE OVERALL?
      
      etas<-as.vector(empty.table.all.sources[,1:3])
      tss<-as.vector(table.current.study[,1:3])
      
      tables.identical<-FALSE
      
      if((sum(etas==tss))==length(etas))tables.identical<-TRUE
      
      #if the table structure is identical to the structure of the empty table from all sources
      #then simply write the counts from the study specific table to the empty table from all sources
      #and then that becomes the study.specific.table for the given study
      
      #if structure not identical match rows in study specific table to rows in array.all.sources using
      #the dimnames.x.num to index the equivalent rows then the dimnames.x to check
      
      if(tables.identical)
      {
        #array.all.sources[,,,ns]<-array.current.study
      }
      
      #if(!tables.identical)
      {
        #set up sequential numeric code for each of the sorted unique values in each dimnames
        dimnames.1<-dimnames(array.all.sources)[[1]]
        dimnames.1.num<-1:length(dimnames.1)
        #cbind(dimnames.1,dimnames.1.num)
        
        dimnames.2<-dimnames(array.all.sources)[[2]]
        dimnames.2.num<-1:length(dimnames.2)
        #cbind(dimnames.2,dimnames.2.num)
        
        dimnames.3<-dimnames(array.all.sources)[[3]]
        dimnames.3.num<-1:length(dimnames.3)
        #cbind(dimnames.3,dimnames.3.num)
        
        #Map row in table from current study with equivalent row in array containing unique
        #values from every study and check that the identified row in the full array contains
        #the same dimnames values as the row in the current study table. If it does not,
        #then rather than trying to second guess all possible ways this could go wrong
        #simply stop processing and ask user to work with study specific tables to
        #create table statistics he/she requires
        
        index.current.study<-rep(NA,length(table.current.study[,1]))
        index.overall<-rep(NA,(dim.vector.all.sources[1]*dim.vector.all.sources[2]*dim.vector.all.sources[3]))
        
        for(oo in 1:length(table.current.study[,1]))
        {
          d1<-table.current.study[oo,1]
          n1<-dimnames.1.num[dimnames.1==d1]
          
          d2<-table.current.study[oo,2]
          n2<-dimnames.2.num[dimnames.2==d2]
          
          d3<-table.current.study[oo,3]
          n3<-dimnames.3.num[dimnames.3==d3]
          
          index.current.study<-oo
          
          #index.overall applies to empty.table.all.sources.col.1 etc which are all of length dim1*dim2*dim3 because
          #this was created before array.all.sources was replicated to include space for all studies. So calculations
          #of index.overall do not need to take account of ns value
          
          index.overall<-n1+dim.vector.all.sources[1]*(n2-1)+
            dim.vector.all.sources[1]*dim.vector.all.sources[2]*(n3-1)
          
          
          count.current.study.current.row<-table.current.study[oo,4]
          
          #check dimnames all match
          d1.a<-empty.table.all.sources.col.1[index.overall]
          
          d2.a<-empty.table.all.sources.col.2[index.overall]
          
          d3.a<-empty.table.all.sources.col.3[index.overall]
          
          
          if(d1.a!=d1||d2.a!=d2||d3.a!=d3)
          {
            return.message=  "Dimensions of tables not behaving sensibly across studies.Please check the data in each study and calculate counts and percentages, yourself, using the counts from the individual studies"
            message(return.message)
            return(return.message) 
          }
          
          #dimension markers match so allocate count to correct cell
          
          #index.overall.extended.over.all.studies applies to array.all.sources which is of length dim1*dim2*dim3*numsources
          #must therefore add dim.vector.all.sources[1]*dim.vector.all.sources[2]*dim.vector.all.sources[3]*(ns-1)
          #to identify correct row in extended array which has the full dimnames structure replicated for each study
          
          index.overall.extended.over.all.studies<-index.overall+
            dim.vector.all.sources[1]*dim.vector.all.sources[2]*dim.vector.all.sources[3]*(ns-1)
          
          array.all.sources[index.overall.extended.over.all.studies]<-count.current.study.current.row
        }
        
        dimnames(array.all.sources)[4]<-list(as.character(1:numsources))
        #print(array.all.sources)
      }#end of tables not identical loop
      
    }
    #end of ns loop

    
    #Combine across studies if requested
    ####################################
    
    combine.array.all.sources<-as.numeric(array.all.sources[,,,1])
    
    if(numsources>1)
    {
      for(ns in 2:numsources)
      {
        combine.array.all.sources<-combine.array.all.sources+as.numeric(array.all.sources[,,,ns])
      }
    }
    
    combine.array.all.sources<-array(data=combine.array.all.sources,dim=dim(array.all.sources)[1:table.dimensions[1]],
                                     dimnames=dimnames(array.all.sources)[1:table.dimensions[1]])

    
    
    
    #return(array.all.sources)
    
  #end of 3 dims loop
  
  ######################################################################
  
  #clean and process output tables
  
  array.all.sources.temp<-array.all.sources
  
  array.all.sources<-as.numeric(array.all.sources.temp)
  
  array.all.sources<-array(data=array.all.sources,dim=dim(array.all.sources.temp),
                           dimnames=dimnames(array.all.sources.temp))

  output.text.temp<-paste0(",TABLES.COMBINED_all.sources_counts=combine.array.all.sources)")
  
if(split){
  for(ns in numsources:1)
  {
    name.array.study<-paste0("array.study.",ns)
    commas.vect<-rep(",",table.dimensions[1])
    commas.vect<-paste(commas.vect,collapse="")
    calltext<-paste0(name.array.study,"<-array.all.sources[",commas.vect,ns,"]")
    eval(parse(text=calltext))
    
    output.text.temp<-paste0(",TABLE_STUDY.",study.names.valid[ns],"_counts=array.study.",ns, output.text.temp)
  }
}

  ##################################################
  #NOW MOVE TO CALCULATE ROW AND COLUMN PROPORTIONS#
  ##################################################
  
    #start with combined table				  
    if(!split){
    combine.array.all.sources.row.props<-combine.array.all.sources
    combine.array.all.sources.col.props<-combine.array.all.sources
    
    for(st in 1:length(stvar.dimnames.unique))
    {
      numrows<-dim(combine.array.all.sources[,,st])[1]
      numcols<-dim(combine.array.all.sources[,,st])[2]
      
      for(nr in 1:numrows)
      {
        sum.row<-sum(combine.array.all.sources[nr,,st],na.rm=TRUE)
        combine.array.all.sources.row.props[nr,,st]<-signif(combine.array.all.sources[nr,,st]/sum.row,3)
      }
      
      for(nc in 1:numcols)
      {
        sum.col<-sum(combine.array.all.sources[,nc,st],na.rm=TRUE)
        combine.array.all.sources.col.props[,nc,st]<-signif(combine.array.all.sources[,nc,st]/sum.col,3)
      }
      
      
    }
    
    calltext2<-paste0("TABLE.COMBINED_row.props<-combine.array.all.sources.row.props")
    calltext3<-paste0("TABLE.COMBINED_col.props<-combine.array.all.sources.col.props")
    TABLE.COMBINED_row.props <- NULL
    TABLE.COMBINED_col.props <- NULL
    eval(parse(text=calltext2))
    eval(parse(text=calltext3))
    
    output.text.temp<-paste0("output.list <- list(TABLES.COMBINED_all.sources_row.props=TABLE.COMBINED_row.props",
                             ",TABLES.COMBINED_all.sources_col.props=TABLE.COMBINED_col.props",
                             output.text.temp)

    eval(parse(text=output.text.temp))
    }

    #######################	
    #study specific tables#
    #######################
   if(split){
    for(ns in numsources:1)
    {
      #	name.array.study<-paste0("array.study.",ns)
      commas.vect<-rep(",",table.dimensions[1])
      commas.vect<-paste(commas.vect,collapse="")
      calltext<-paste0("study.specific.table<-array.all.sources[",commas.vect,ns,"]")
      study.specific.table<-NULL
      eval(parse(text=calltext))
      
      study.specific.table.row.props<-study.specific.table
      study.specific.table.col.props<-study.specific.table
      
      
      for(st in 1:length(stvar.dimnames.unique))
      {
        numrows<-dim(study.specific.table[,,st])[1]
        numcols<-dim(study.specific.table[,,st])[2]
        
        
        for(nr in 1:numrows)
        {
          sum.row<-sum(study.specific.table[nr,,st],na.rm=TRUE)
          study.specific.table.row.props[nr,,st]<-signif(study.specific.table[nr,,st]/sum.row,3)
        }
        for(nc in 1:numcols)
        {
          sum.col<-sum(study.specific.table[,nc,st],na.rm=TRUE)
          study.specific.table.col.props[,nc,st]<-signif(study.specific.table[,nc,st]/sum.col,3)
        }
      }
      
      calltext4<-paste0("TABLE.STUDY_row.props.",ns,"<-study.specific.table.row.props")
      calltext5<-paste0("TABLE.STUDY_col.props.",ns,"<-study.specific.table.col.props")
      
      #print(calltext4)
      #print(calltext5)
      
      eval(parse(text=calltext4))
      eval(parse(text=calltext5))
      
    }
    
  
    for(ns in numsources:1)
    {
      if(ns>1)
      {	
        output.text.temp<-paste0(",TABLE.STUDY.",study.names.valid[ns],"_row.props=TABLE.STUDY_row.props.",ns,",",
                                 "TABLE.STUDY.",study.names.valid[ns],"_col.props=TABLE.STUDY_col.props.",ns,
                                 output.text.temp)
        
      }
      else	
      {
        output.text.props.counts<-
          paste0("output.list<-list(TABLE.STUDY.",study.names.valid[ns],"_row.props=TABLE.STUDY_row.props.",ns,",",
                 "TABLE.STUDY.",study.names.valid[ns],"_col.props=TABLE.STUDY_col.props.",ns,
                 output.text.temp)
        
      }
    }

    eval(parse(text=output.text.props.counts))
  }
    return.list.first<-list(output.list=output.list,validity.message=validity.message)
  
  #END second dim=3 loop


  #################################################
  # Setup on.exit() to restore options 'warn' value
  #################################################
  
  old_warn_option <- base::getOption("warn")
  on.exit(base::options(warn = old_warn_option), add = TRUE)

  ################################
  #NOW UNDERTAKE CHISQUARED TESTS#
  ################################

  
    #Suppress.chisq.warnings by default
    if(suppress.chisq.warnings)
    {
      options(warn=-1)
    }
    ##################
    #Combined studies#
    ##################
    
    
      results_table <- data.frame()
      numtests<-dim(combine.array.all.sources)[table.dimensions[1]]
      
      chisq.list.temp<-")"
      if(!split){
      for(nt in numtests:1)
      {
        chisqtext<-paste0("chisq.test_TABLES.COMBINED.",nt,"<-stats::chisq.test(combine.array.all.sources[,,nt])")
        eval(parse(text=chisqtext))
        data.name.change <- paste0("chisq.test_TABLES.COMBINED.", nt, "$data.name <- stvar.all.unique.levels[nt]") # tk modified 
        eval(parse(text=data.name.change)) # tk modified
        chisq.list.temp<-paste0(",chisq.test_TABLES.COMBINED_all.sources_counts_table.",nt,"=chisq.test_TABLES.COMBINED.",nt,chisq.list.temp)	
        
      }
      chisq.list.temp <- paste0("chisq.tests <- list(", sub("^,", "", chisq.list.temp))
      eval(parse(text=chisq.list.temp))
      
      results_table <- data.frame(
        lsoa11cd = vapply(chisq.tests, function(x) x$data.name, character(1)),
        X.squared = vapply(chisq.tests, function(x) unname(x$statistic), numeric(1)),
        df = vapply(chisq.tests, function(x) unname(x$parameter), numeric(1)),
        p.value = vapply(chisq.tests, function(x) x$p.value, numeric(1)),
        datasource = "COMBINED",
        row.names = NULL)
      
      if(names_region.valid){
        plot_matrix[[1]] <- shape_list[[1]] |>
          dplyr::left_join(results_table, by = "lsoa11cd") |>
          sf::st_as_sf()
      }
      
      }
      ##################
      #Separate studies#
      ##################
      if(split){
      for(ns in numsources:1)
      {
        input.calltext<-paste0("input.array.source.specific<-array.study.",ns)
        input.array.source.specific <- NULL
        eval(parse(text=input.calltext))
        
        numtests<-dim(input.array.source.specific)[table.dimensions[1]]
        
        
        for(nt in numtests:1)
        {
          
           if(nt>1)
           {
             chisqtext<-paste0("chisq.test_TABLE.STUDY.",ns,"_counts.",nt,"<-stats::chisq.test(input.array.source.specific[,,nt])")
             eval(parse(text=chisqtext))
             data.name.change <- paste0("chisq.test_TABLE.STUDY.", ns,"_counts.",nt, "$data.name <- stvar.all.unique.levels[nt]") # tk modified 
             eval(parse(text=data.name.change)) # tk modified
             chisq.list.temp<-paste0(",chisq.test_TABLE.STUDY.",study.names.valid[ns],"_counts_table.",nt,"=chisq.test_TABLE.STUDY.",ns,"_counts.",nt,chisq.list.temp)	
           }
           else
           {
            
            chisqtext<-paste0("chisq.test_TABLE.STUDY.",ns,"_counts.",nt,"<-stats::chisq.test(input.array.source.specific[,,nt])")
            eval(parse(text=chisqtext))
            data.name.change <- paste0("chisq.test_TABLE.STUDY.", ns,"_counts.",nt, "$data.name <- stvar.all.unique.levels[nt]") # tk modified 
            eval(parse(text=data.name.change)) # tk modified
            chisq.list.text<-paste0("chisq.tests<-list(chisq.test_TABLE.STUDY.",study.names.valid[ns],"_counts_table.",nt,"=chisq.test_TABLE.STUDY.",ns,"_counts.",nt,chisq.list.temp)
           }
        }
        eval(parse(text=chisq.list.text))
        results_server <- data.frame(
          lsoa11cd = vapply(chisq.tests, function(x) x$data.name, character(1)),
          X.squared = vapply(chisq.tests, function(x) unname(x$statistic), numeric(1)),
          df = vapply(chisq.tests, function(x) unname(x$parameter), numeric(1)),
          p.value = vapply(chisq.tests, function(x) x$p.value, numeric(1)),
          datasource = study.names.valid[ns],
          row.names = NULL)
        
       if(names_region.valid){
        plot_matrix[[ns]] <- shape_list[[ns]] |>
          dplyr::left_join(results_server, by = "lsoa11cd") |>
          sf::st_as_sf()
       }
       else
        {
        results_table <- rbind.data.frame(results_table, results_server)
        }
      }
      #END ns loop
      }
      
      #If warnings suppressed now return to default
      if(suppress.chisq.warnings)
      {
        options(warn=0)
      }
    

# Plotting 
  if (!draw.plot && !names_region.valid) {
    
    return(results_table)
    
  } else if (!draw.plot && names_region.valid) {
    
    plot.matrix <- do.call(rbind, plot_matrix) |>
      sf::st_drop_geometry()
    
    return(plot.matrix)
    
  } else if (draw.plot && names_region.valid) {
    
    plot.matrix <- do.call(rbind, plot_matrix)
    
    plotresult <- ggplot2::ggplot(plot.matrix) +
      ggplot2::geom_sf(
        ggplot2::aes(fill = X.squared)
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
      ggplot2::labs(fill = "X.squared") +
      ggplot2::facet_wrap(~datasource) +
      ggplot2::theme_minimal()
    
    return(plotresult)
  }
  }
  

