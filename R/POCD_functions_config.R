#' check rawdata for consistence
#'
#' @param Input Input table
#' @param required.colnames colnames, which are required
#' @return throw out an error if other colnames are needed
#' @export check.colnames.FUN

check.colnames.FUN <- function(Input, # = Config$info.table,
                               required.colnames){# = required.colnames.config.file

  if(any(!(required.colnames %in% colnames(Input)))){
    stop(paste0("Unknown column names in info file detected. The exspected column names are: ", paste(required.colnames, collapse =", "),"."))
  }

}

############################################################################

#' grouping into subsets
#' @description This function creates a list of separate data frames for each group and timepoint. By default no further arguments are necessary, as all input variables have pre-defined values.
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom dplyr rename_at
#' @importFrom dplyr select
#' @importFrom rlang .data
#' @param Input data frame that contains the complete data
#' @param Timepoint.names names of timepoints
#' @param Controls name of control group
#' @param Patients names of treatment groups
#' @param names.col names of all parameter
#' @return list of separate data frames for each group and timepoint.
#' @export grouping.subsets.FUN

grouping.subsets.FUN = function(

                        Input,               #= Config$CogDat$Table,
                        Timepoint.names,     #= Config$parameter$names.Timepoints,
                        Controls,            #= Config$parameter$names.Group$Control,
                        Patients,            #= Config$parameter$names.Group$Patients,
                        names.col)           #= Config$parameter$colnames)
                        #names.CG.cleaning   = Config$parameter$names.CG.cleaning,
                        #names.NCD           = Config$parameter$names.NCD)
                        {

  CTPs = list()

  for (i in Timepoint.names) {
    CTPs[[paste0("Patients.", i)]] = Input %>% subset(Group %in% Patients) %>%
      dplyr::select(.data$SubjectID, .data$Group, grep( i, names.col, value = T))

  }

 for (i in Timepoint.names) {
    CTPs[[paste0("Controls.", i)]] = Input %>% subset(Group %in% Controls) %>%
      dplyr::select(.data$SubjectID, .data$Group, grep( i, names.col, value = T))
 }

  return(CTPs)

}

#############################################################################

#' logarythmic transformation of CTP data
#'
#' @param Input list of data frames containing the CTP information
#' @param CTP.names names of CTP tests
#' @param transformation  info which columns should be transform
#' @return list of dataframes with logarythmic transformed CTP values
#' @export grouping.subsets.FUN

transformation.FUN <- function(Input, #= Config$CogDat$CTPs.control.clean,
                               CTP.names, #= Config$parameter$colnames.CTPs,
                               transformation) #= Config$parameter$CTP.transformation
  {
  CTPs <- Input
  rownames(CTPs) <- Input$SubjectID
  df   <- CTPs[colnames(Input) %in% CTP.names]
  # CTPs   <- CTPs[rowSums(is.na(df)) != ncol(df), ]
  CTPs[ , sapply(rownames(transformation)[transformation[1] == TRUE], function(x) grep(x, colnames(CTPs)))] <-

    apply(CTPs[sapply(rownames(transformation)[transformation[1] == TRUE], function(x) grep(x, colnames(CTPs)))], 2,log)


  return(CTPs)
}

#######################################################################################

#' remove problematic controls group measures
#'
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom dplyr rename_at
#' @importFrom dplyr select
#' @param Input list of data frames with CTPs and cleaning variable
#' @param patterns string pattern, which match with problematic measures
#' @param CG.cleaning name of variable coding for reason not to use
#' @return list of data frames with problematic control group measures removed
#' @export cleaning.cg.FUN

cleaning.cg.FUN = function(Input,# = Config$CogDat$CTPs.raw,
                           patterns = c("do not use"),
                           CG.cleaning)# = Config$parameter$colnames.Comment)
  {
  CTPs = Input

  for (i in names(CTPs)) {
    df = CTPs[[i]]

   if(length(grep(patterns, df[,colnames(df) %in% CG.cleaning])) > 0){

    CTPs[[i]] =  df[-grep(patterns, df[,colnames(df) %in% CG.cleaning]), ]
  } else{
    CTPs[[i]] =  df
  }
}

for (i in grep( "Patients",names(CTPs), value = T)[-1]) {
  CTPs[[i]] = CTPs[[i]][CTPs[[i]]$SubjectID %in% CTPs$Patients.T0$"SubjectID", ]
}

for (i in grep( "Controls",names(CTPs), value = T)[-1]) {
  CTPs[[i]] = CTPs[[i]][CTPs[[i]]$SubjectID %in% CTPs$Controls.T0$"SubjectID", ]
}

return(CTPs)

}

#######################################################################################

#' replace missings due to physical impairment with worst performance
#'
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom dplyr rename_at
#' @importFrom dplyr select
#' @param Input ist of data frames with CTPs and information about health condition
#' @param CTP.names names of CTP tests
#' @param pattern string pattern, which match with problematic measures
#' @param polarity information about the polarity of the CTP tests
#' @return  list of data frames with worst performanced imputated values
#' @export replace.worst.FUN

replace.worst.FUN = function(Input,# = Config$CogDat$CTPs.control.clean,
                             CTP.names,# = Config$parameter$colnames.CTPs,
                             polarity,# = Config$parameter$PoledVars,
                             pattern = "worst performance")
  {

  replace.FUN = function(In){
  CTPs <- In
  df = In[colnames(In) %in% CTP.names]
  # get performance range per CTP
  performance.range = apply(df, 2, range, na.rm = TRUE)

  # define worst performance per CTP
  worst.performance = performance.range[1, ]
  worst.performance[polarity] = performance.range[2, polarity]

  # replace missings due to physical impairment with worst performance
  for (i in colnames(df)) {
    CTPs[,i] <- ifelse(In[grep("Comment", colnames(In))] == pattern & is.na(In[,i]),
                      worst.performance[i], In[,i])
    colnames(CTPs[,i]) <- i

  }
  CTPs <- as.data.frame.matrix(CTPs)

  return(CTPs)
}

CTPs.worst <- lapply(Input, replace.FUN)

return(CTPs.worst)
}

##########################################################################

#' random forest imputation - version 1
#'
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom dplyr rename_at
#' @importFrom dplyr select
#' @param Input list of data frames with CTPs and cleaning variable
#' @param CTP.names names of CTP tests
#' @param seed set seed for reproducible results
#' @return  list of data frames with rendom forest imputed values
#' @export random_forest.FUN

random_forest.FUN = function(Input,#            = Config$CogDat$CTPs.participate,
                             seed             = 145,
                             CTP.names)#        = Config$parameter$colnames.CTPs)
  {

  CTPs <- Input
  #sort data according to SubjectID
  CTPs <- lapply(CTPs, function(i) i[order(i$SubjectID,decreasing = F),])
  CTPs <- lapply(CTPs, function(x) x[!is.na(x$SubjectID),])
  CTPs.red   <- lapply(CTPs, function(x) {data.frame(x[c(colnames(x) %in% CTP.names)], row.names = x$SubjectID)})

  #check if CTPs are numeric or integer and convert integer to factor; necessary for randomforest imputation
  # function to check if Colum is numeric or integer

  is.integer.FUN = function(x) {all(x == round(x), na.rm = T)}

  CTPs.red.missforest  <-   lapply(CTPs.red, function(df){

    check.int <- apply(df, 2,FUN = is.integer.FUN)
    df[,check.int] <- lapply(df[,check.int] , factor)

    if(seed != 0 ){ set.seed(seed)}
    rf <- missForest::missForest(df, variablewise = TRUE)

    # take only imputed data frames without MSE
    rf <- rf$ximp

    #convert all columns to numeric
    rf[,1:length(rf)] <- lapply(rf, as.character)
    rf[,1:length(rf)] <- lapply(rf, as.numeric)

    rf$SubjectID = rownames(rf)

    return(rf)

  })


  for (i in grep("Patients",names(CTPs.red.missforest), value = T)[-1]) {
    CTPs.red.missforest[[i]] = CTPs.red.missforest[[i]][CTPs.red.missforest[[i]]$SubjectID %in% CTPs.red.missforest$Patients.T0$"SubjectID", ]

  }

  for (i in grep("Controls",names(CTPs.red.missforest), value = T)[-1]) {
    CTPs.red.missforest[[i]] = CTPs.red.missforest[[i]][CTPs.red.missforest[[i]]$SubjectID %in% CTPs.red.missforest$Controls.T0$"SubjectID", ]
  }

  CTPs.red.missforest <- lapply(CTPs.red.missforest, function(x) {
    row.names(x) <- x$SubjectID
    x[,-grep("SubjectID", colnames(x))]})

  return(CTPs.red.missforest)

}
