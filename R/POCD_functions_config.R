#' check rawdata for consistence
#'
#' @param Input               Input table
#' @param required.colnames   colnames, which are required
#' @export

check.colnames.FUN <- function(Input             = Config$info.table,
                              required.colnames = required.colnames.config.file) {

  if(any(!(required.colnames %in% colnames(Input)))){
    stop(paste0("Unknown column names in info file detected. The exspected column names are: ", paste(required.colnames, collapse =", "),"."))
  }

}



############################################################################

#' check rawdata for consistence
#'
#' @param Config.colnames    given colnames of the CTP_information_table
#' @param Config.table       CTP_information_table
#' @param rawdata.colnames   given colnames of the rawdata_table
#' @param rawdata.table      rawdata_table
#' @param Timepoint.names     names of the test timepoints
#' @param CTP.number         number of CTP tests
#' @export check.rawdata.FUN



check.rawdata.FUN <- function(Config.colnames  = colnames.Config, Config.table = Config$Config.Table,
                          rawdata.colnames = colnames.CogDat, rawdata.table = Config$CogDat$Input,
                          Timepoint.names  =  Timepoints, CTP.number = CTPs){

  ## check CTP-Inputfile
  # check column names
  if(any(!(Config.colnames %in% names(Config.table)))){
    stop(paste0("Unknown column names in CTP file detected. The exspected column names are: ", paste(Config.colnames, collapse =", "),"."))
  } else {print(paste("Column names of CTP-table are: ", paste(Config.colnames, collapse =", ")))}

   # check colums due to datatype and input
  if(any(!(Config.table[,2] %in% c("-","+")))){
    stop("Wrong input in Column polarity in CTP-Table. Input need to be \"+\" or \"-\" values.")}

  if(!(is.character(Config.table[,1]))){
    stop("Wrong input in Column \"CTP\" in CTP-Table. Input should be of mode character.")}

  # check number of tests
  if(CTP.number != nrow(Config.table)){
    stop(paste0("Your Input: number of CTPs = ",CTP.number," does not match with number of rows in CTP-table"))
    } else{print(paste0("There are ",CTP.number ," CTPs "))}

  #### check SPSS-Inputfile
  # check column names
  if(any(!(rawdata.colnames %in% names(rawdata.table)))){
    stop(paste0("Unknown columnnames in input file detected. ",
                "The exspected column names are: ", paste(rawdata.colnames[!(rawdata.colnames %in% names(rawdata.table))], collapse =", "),"."))
  } else {print(paste("Column names of CTP file are: ", paste(rawdata.colnames, collapse =", ")))}

    # check coloums due to datatype and input
  if(any(!(sapply(rawdata.table[paste0(rep(Timepoint.names, each = CTP.number),"_CTP", 1:CTP.number)], is.numeric)))){
    stop("Wrong input in CTP Columns. The type should be numeric.")
  }
}

############################################################################

#' rename rawdata
#'
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom dplyr rename_at
#' @importFrom dplyr select
#' @param Config.colnames    given colnames of the CTP_information_table
#' @param Config.table       CTP_information_table
#' @param rawdata.colnames   given colnames of the rawdata_table
#' @param rawdata.table      rawdata_table
#' @param Timepoint.names     names of the test timepoints
#' @param CTP.number         number of CTP tests
#' @param Control name of Controlgroup
#' @param Patients names of treatment groups
#' @param health name of variable coding poor health condition
#' @param dead name of variable coding for deceased
#' @param assessed name of variable coding for successfull examination
#' @export rename.rawdata.FUN



rename.rawdata.FUN <- function(Config.colnames  = colnames.Config, Config.table = Config$Config.Table,
                           rawdata.colnames = colnames.CogDat, rawdata.table = Config$CogDat$Input,
                           Timepoint.names =  Timepoints, CTP.number = CTPs,
                           Control = Controlgroup, Patients = Patientsgroup,
                           health = name.poor.health, dead = name.dead, assessed = name.assessed,
                           method = name.method){

  # rename CTP_information_table
  Config.table %<>% rename_at(Config.colnames, ~ c("variable", "polarity"))

  # rename SPSS-table
  Timepoint.names <- paste0("T", 1:length(Timepoint.names)-1)

  col.SubjectID   = "SubjectID"
  col.group       = "Group"
  col.comment1    = paste0(Timepoint.names,"_comment1")
  col.comment2    = paste0(Timepoint.names,"_comment2")
  col.adjustment  = paste0(Timepoint.names,"_Adjustment_recommendation")
  col.CTPs        = paste0(rep(Timepoint.names, each = CTP.number),"_CTP", 1:CTP.number)
  col.NCD         = paste0(Timepoint.names,"_func_subj")

  colnames.new <- c(col.SubjectID,
                       col.group,
                       col.comment1,
                       col.comment2,
                       col.adjustment,
                       col.CTPs,
                       col.NCD)


    rawdata.table <- rawdata.table %>% dplyr::select(rawdata.colnames) %>% rename_at(rawdata.colnames, ~ colnames.new)
  print(paste0("Renamed column names are: ", paste(colnames(rawdata.table) ,collapse = ", ")))

  #sort rawdata according to SubjectID
  rawdata.table <- rawdata.table[order(rawdata.table$SubjectID,decreasing = F),]
  #remove whitespace from SubjectID
  rawdata.table$SubjectID <- trimws(rawdata.table$SubjectID)
  rawdata.table <- data.frame(rawdata.table, stringsAsFactors = F)


  #get parameter for objects
  parameter <- list()
  parameter$names.CTPs <- col.CTPs
  parameter$names.Timepoints <- Timepoint.names
  parameter$names.Group$Controls <- Control
  parameter$names.Group$Patients <- Patients
  parameter$names.comment1 <- col.comment1
  parameter$names.comment2 <- col.comment2
  parameter$names.CG.cleaning <- col.adjustment
  parameter$names.NCD <- col.NCD

  parameter$OldCTPs <- data.frame(sapply(Timepoint.names, function(i) vars_select(parameter$names.CTPs, starts_with(i))),
                                        row.names = NULL,
                                        stringsAsFactors = F
                                        )

  # set polarity of CTPs -> Config$PoledVars
  parameter$PoledVars = grep("-", Config.table$polarity, value = FALSE)

  parameter$N.Timepoints = length(Timepoint.names)
  parameter$N.CTPs = CTP.number
  parameter$names.poor.health = health
  parameter$names.dead = dead
  parameter$names.assessed = assessed
  parameter$names.abort = c("comment1", "comment2")
  parameter$POCDmethod = method
  parameter$N.Subgroups = parameter$N.Timepoints * 2 #(length(parameter$names.Group))

  return(list(Config.table, rawdata.table, parameter))

}

#############################################################################

#' grouping into subsets
#'
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom dplyr rename_at
#' @importFrom dplyr select
#' @param Input              data frame that contains the complete data
#' @param Timepoints         names of timepoints
#' @param Controls name of control group
#' @param Patients names of treatment groups
#' @param CTP.Vars data frame with names of CTP variables
#' @param abort.name names of columns containing reasons for abort
#' @param CG.cleaning.name names of columns containing reasons for abort due to technical or administrative problems
#' @export grouping.subsets.FUN

# create separate data files for each group and time point while choosing only those variables, that are necessary for POCD calculation (only CTPs).
grouping.subsets.FUN = function(

                        Input               = Config$CogDat$Table,
                        Timepoints          = Config$parameter$names.Timepoints,
                        Controls            = Config$parameter$names.Group$Control,
                        Patients            = Config$parameter$names.Group$Patients,
                        names.col           = Config$parameter$colnames)
                        #names.CG.cleaning   = Config$parameter$names.CG.cleaning,
                        #names.NCD           = Config$parameter$names.NCD)
                        {
  # This function creates a list of separate data frames for each group and timepoint. By default no further arguments are necessary, as all input variables have pre-defined values.



  CTPs = list()

  for (i in Timepoints) {
    CTPs[[paste0("Patients.", i)]] = Input %>% subset(Group %in% Patients) %>%
      dplyr::select(SubjectID, Group, grep( i, names.col, value = T))

  }

 for (i in Timepoints) {
    CTPs[[paste0("Controls.", i)]] = Input %>% subset(Group %in% Controls) %>%
      dplyr::select(SubjectID, Group, grep( i, names.col, value = T))
 }

  return(CTPs)

}

#######################################################################################

#' remove problematic controls group measures
#'
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom dplyr rename_at
#' @importFrom dplyr select
#' @param Input        data frame with CTPs and cleaning variable
#' @param Timepoints    indicator of test time point
#' @param patterns string pattern, which match with problematic measures
#' @param CG.cleaning name of variable coding for reason not to use
#' @return CTPs: data frame with problematic control group measures removed
#' @export cleaning.cg.FUN

cleaning.cg.FUN = function(Input = Config$CogDat$CTPs.raw,
                           Timepoints = Config$parameter$names.Timepoints,
                           patterns = c("do not use"),
                           CG.cleaning = Config$parameter$colnames.Comment) {

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
#' @param Input        data frame with CTPs and cleaning variable
#' @param CTP.names         names of CTP tests
#' @param poor.health name of variable coding poor health condition
#' @param patterns string pattern, which match with problematic measures
#' @return  data frame with worst performanced imputated values
#' @export replace.worst.FUN

replace.worst.FUN = function(Input = Config$CogDat$CTPs.control.clean,
                             CTP.names = Config$parameter$colnames.CTPs,
                             polarity = Config$parameter$PoledVars,
                             pattern = "worst performance") {

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
    CTPs[,i] <- ifelse(In[grep("Comment", colnames(In))] == pattern,
                      worst.performance[i], In[,i])
    colnames(CTPs[,i]) <- i

  }
  CTPs <- as.data.frame.matrix(CTPs)

  return(CTPs)
}

CTPs.worst <- lapply(Input, replace.FUN)

# CTPs.worst$doc <- lapply(1:Config$parameter$N.Subgroups, function(x){
#   inp = Input[[x]]
#   df = inp[colnames(inp) %in% CTP.names]
#   data.frame(SubjectID =  inp$SubjectID,
#              Timepoints =  unlist(strsplit(split = "\\.",names(Input[x])))[2],
#              Groups     =  unlist(strsplit(split = "\\.",names(Input[x])))[1],
#              worst_performance   = sapply(1:nrow(inp), function(x) any((inp[x,grep("comment1", colnames(inp))] == poor.health | inp[x,grep("comment2", colnames(inp))] == poor.health ) & is.na(inp[x, colnames(df)]))),
#              worst_performance.count = sapply(1:nrow(inp), function(x) sum((inp[x,grep("comment1", colnames(inp))] == poor.health | inp[x,grep("comment2", colnames(inp))] == poor.health ) & is.na(inp[x, colnames(df)])))
#   )
# })
#
# CTPs.worst$doc <- plyr::join_all(CTPs.worst$doc,type = "full")

return(CTPs.worst)
}

##########################################################################

#' remove cases with no data and deceased patients
#'
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom dplyr rename_at
#' @importFrom dplyr select
#' @param Input        data frame with CTPs and cleaning variable
#' @param CTP.names         names of CTP tests
#' @param modi content of removed data : "full" - all rows with complete NA?s were removed, "onlyDead" - only dead patients were removed
#' @return  data frame with complete NA-rows removed
#' @export remove.na.FUN


remove.na.FUN = function(Input = Config$CogDat$CTPs.worst.perform,
                         CTP.names = Config$parameter$names.CTPs,
                         modi = "full"
                         ) {

  remove.FUN <- function(Inp){
  CTPs <- Inp
if(modi == "full"){
  df <- Inp[colnames(Inp) %in% CTP.names]
  CTPs = CTPs[rowSums(is.na(df)) != ncol(df), ]

  return(CTPs)
}

  if(modi == "onlyDead"){
    CTPs <- CTPs[CTPs[,grep("comment1" , colnames(CTPs), value = T)] != Config$parameter$names.dead,]
    return(CTPs)
  }

  }

  CTPs.remove <- lapply(Input[1:Config$parameter$N.Subgroups], remove.FUN)
  CTPs.remove$doc <- lapply(1:Config$parameter$N.Subgroups, function(x){
    inp <- Input[[x]]
    df = inp[colnames(inp) %in% CTP.names]
    data.frame(
      SubjectID =  inp$SubjectID,
      Timepoints =  unlist(strsplit(split = "\\.",names(Input[x])))[2],
      Groups     =  unlist(strsplit(split = "\\.",names(Input[x])))[1],
      missing    = if(modi == "full"){ rowSums(is.na(df)) == ncol(df)} else  if(modi == "onlyDead"){ 0 }



    )
  })


  CTPs.remove$doc <- plyr::join_all(CTPs.remove$doc,type = "full")

  return(CTPs.remove)
}

# remove.na.Adj.FUN = function(df1, df0 = df1[,1:Config$N$CTPs]) {
#   clean.df = data.frame(na = ifelse(rowSums(is.na(df0)) != ncol(df0), FALSE, TRUE))
#   return(clean.df)
# }

##########################################################################

#' random forest imputation - version 1
#'
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom dplyr rename_at
#' @importFrom dplyr select
#' @param Input        data frame with CTPs and cleaning variable
#' @param Timepoint.names     names of the test timepoints
#' @param CTP.names         names of CTP tests
#' @param seed
#' @return  data frame with rendom forest imputed values
#' @export random_forest.FUN


random_forest.FUN = function(Input            = Config$CogDat$CTPs.participate,
                             Timepoints.names = Config$parameter$names.Timepoints,
                             seed             = 145,
                             CTP.names        = Config$parameter$colnames.CTPs) {

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


  # CTPs.red.missforest$doc <- lapply(1:Config$parameter$N.Subgroups, function(x){
  #   inp = Input[[x]]
  #   df = inp[colnames(inp) %in% CTP.names]
  #   data.frame(SubjectID =   inp$SubjectID,
  #              Timepoints =  unlist(strsplit(split = "\\.",names(Input[x])))[2],
  #              Groups     =  unlist(strsplit(split = "\\.",names(Input[x])))[1],
  #              random_forest       = sapply(1:nrow(df), function(x) any(is.na(df[x,]))),
  #              random_forest.count = sapply(1:nrow(df), function(x) sum(is.na(df[x,])))
  #   )
  #   })
  #
  #
  # CTPs.red.missforest$doc <- plyr::join_all(CTPs.red.missforest$doc,type = "full")

  return(CTPs.red.missforest)

}

##########################################################################

#' seperate documentation and reduce CTP-tables to CTP values only
#'
#' @importFrom magrittr %<>%
#' @importFrom magrittr %>%
#' @importFrom dplyr rename_at
#' @importFrom dplyr select
#' @param Input        data frame with CTPs and cleaning variable
#' @return  data frame with complete NA-rows removed
#' @export doc.FUN


doc.FUN <- function(Input = Config$CogDat){

  names <- grep("CTPs", names(Input), value = T)

 #lapply(names, function(x) {

  documentation <- lapply(names, function(x){ Input[[x]][9] })
  documentation <- documentation[sapply(1:length(documentation), function(x) length(documentation[[x]][[1]])>0)]
  documentation <- unlist(documentation, recursive = F)

  Input$documentation <- plyr::join_all(documentation,type = "full", by = c("SubjectID", "Timepoints", "Groups"))
  Input$documentation <- plyr::join_all(list(Input$documentation, Input$Table[,c("SubjectID", "Group")]), by = "SubjectID") %>% dplyr::rename("Treatment" = Group)


  Input[paste0(names, ".interim")] <- Input[names]
  Input[names] <- lapply(names, function(x) Input[[x]][grep(paste(Config$parameter$names.Timepoints,collapse = "|"),names(Input[[x]]),value = T )])
  Input[names] <- lapply(names, function(x) {
    lapply(Input[[x]], function(i) {

      i <- i[!is.na(i$SubjectID),]
      rownames(i) <- i$SubjectID
     i <- i[,colnames(i) %in% Config$parameter$names.CTPs]
     return(i)

    })

  }  )


  return(Input)
   }




###not implemented yet
#
# ```{r random forest imputation - version 2, eval=FALSE, include=FALSE}
# # this version performs a random forest imputation for the whole data set (all time points combined)
#
# # function for removing cases with no measures at a time point
# remove.na.2.FUN = function(df.full, df.clean){
#   df.clean2 = df.full[rownames(df.full) %in% rownames(df.clean), ]
#   return(df.clean2)
# }
#
#
# # create new data.frame with all time points
# Config$CogDat$all.tp.CTPs$Patients = subset(Config$CogDat$Table, Group != Config$Group$Controlgroup) %>%
#   select((one_of(c(Config$Varnames$OldCTPs, "Group"))))
#
# Config$CogDat$all.tp.CTPs$Controls = subset(Config$CogDat$Table, Group == Config$Group$Controlgroup) %>%
#   select((one_of(c(Config$Varnames$OldCTPs, "Group"))))
#
# Config$CogDat$all.tp.CTPs$combined = rbind(Config$CogDat$all.tp.CTPs$Patients, Config$CogDat$all.tp.CTPs$Controls)
#
#
# # impute missing values with random forest model
# Config$CogDat$missForest.all.tp.CTPs = missForest(Config$CogDat$all.tp.CTPs$combined, variablewise = TRUE)
#
#
# # create subsets for every group and time point
# Config$CogDat$missForest.all.tp.CTPs.grouped = grouping.FUN(Input = Config$CogDat$missForest.all.tp.CTPs$ximp)
#
#
# # remove all patients at that have no data
# buffer = list()
#
# for (i in 1:length(Config$CogDat$missForest.all.tp.CTPs.grouped)) {
#   buffer[[i]] = remove.na.2.FUN(Config$CogDat$missForest.all.tp.CTPs.grouped[[i]], Config$CogDat$clean.CTPs[[i]])
# }
#
# names(buffer) = names(Config$CogDat$missForest.all.tp.CTPs.grouped)
#
#
# # fill up follow up measures with empty rows from baseline (necessary for some calculations)
# Config$CogDat$clean.missForest.all.tp.CTPs.grouped = lapply(buffer, id.FUN)
#
# # patients
# for (i in 2:Config$parameter$N.Timepoints) {
#   Config$CogDat$clean.missForest.all.tp.CTPs.grouped[[i]] = merge(Config$CogDat$clean.missForest.all.tp.CTPs.grouped[[i]],
#                                                                   Config$CogDat$clean.missForest.all.tp.CTPs.grouped$Patients.T0,
#                                                                   by = "SubjectID",
#                                                                   all.y = TRUE)
#
#   # make SubjectID the rowname again
#   rownames(Config$CogDat$clean.missForest.all.tp.CTPs.grouped[[i]]) =
#     Config$CogDat$clean.missForest.all.tp.CTPs.grouped[[i]]$SubjectID
#
#   # remove unneeded columns
#   Config$CogDat$clean.missForest.all.tp.CTPs.grouped[[i]] =
#     Config$CogDat$clean.missForest.all.tp.CTPs.grouped[[i]][, c(2:(Config$N$CTPs + 1))]
# }
#
# # controls
# for (i in (Config$parameter$N.Timepoints + 2):(Config$parameter$N.Timepoints + Config$parameter$N.Timepoints)) {
#   Config$CogDat$clean.missForest.all.tp.CTPs.grouped[[i]] = merge(Config$CogDat$clean.missForest.all.tp.CTPs.grouped[[i]],
#                                                                   Config$CogDat$clean.missForest.all.tp.CTPs.grouped$Controls.T0,
#                                                                   by = "SubjectID",
#                                                                   all.y = TRUE)
#
#   # make SubjectID the rowname again
#   rownames(Config$CogDat$clean.missForest.all.tp.CTPs.grouped[[i]]) =
#     Config$CogDat$clean.missForest.all.tp.CTPs.grouped[[i]]$SubjectID
#
#   # remove unneeded columns
#   Config$CogDat$clean.missForest.all.tp.CTPs.grouped[[i]] =
#     Config$CogDat$clean.missForest.all.tp.CTPs.grouped[[i]][, c(2:(Config$N$CTPs + 1))]
# }
#
# # keep only columns with CTPs
# Config$CogDat$clean.missForest.all.tp.CTPs.grouped[[1]] =
#   Config$CogDat$clean.missForest.all.tp.CTPs.grouped[[1]][, c(1:Config$N$CTPs)]
# Config$CogDat$clean.missForest.all.tp.CTPs.grouped[[Config$parameter$N.Timepoints+1]] =
#   Config$CogDat$clean.missForest.all.tp.CTPs.grouped[[Config$parameter$N.Timepoints+1]][, c(1:Config$N$CTPs)]
#
# rm(buffer, i)
# ```
#
#
# ```{r random forest imputation - version 3, eval=FALSE, include=FALSE}
# # this version performs a random forest imputation separated for patients and controls and for pre- and postoperative (combined) time points
#
# # create separate data.frames with pre- and postoperative data
# Config$CogDat$prepo.CTPs$Patients.preop = subset(Config$CogDat$Table, Group != Config$Group$Controlgroup) %>%
#   select((one_of(Config$Varnames$OldCTPs[, 1])))
#
# Config$CogDat$prepo.CTPs$Patients.postop = subset(Config$CogDat$Table, Group != Config$Group$Controlgroup) %>%
#   select((one_of(c(Config$Varnames$OldCTPs[, -1], "Group"))))
#
# Config$CogDat$prepo.CTPs$Controls.preop = subset(Config$CogDat$Table, Group == Config$Group$Controlgroup) %>%
#   select((one_of(Config$Varnames$OldCTPs[, 1])))
#
# Config$CogDat$prepo.CTPs$Controls.postop = subset(Config$CogDat$Table, Group == Config$Group$Controlgroup) %>%
#   select((one_of(c(Config$Varnames$OldCTPs[, -1], "Group"))))
#
#
# # impute missing values with random forest model
# Config$CogDat$missForest.prepo.CTPs = lapply(Config$CogDat$prepo.CTPs, missForest, variablewise = TRUE)
#
#
# # create subsets for every group and time point
# Config$CogDat$missForest.prepo.CTPs.grouped.postop = rbind(Config$CogDat$missForest.prepo.CTPs$Patients.postop$ximp,
#                                                            Config$CogDat$missForest.prepo.CTPs$Controls.postop$ximp)
#
# # Config$CogDat$missForest.prepo.CTPs.grouped.postop = grouping.FUN(Input = Config$CogDat$missForest.prepo.CTPs.grouped.postop,
# #                                                            Timepoints.Controls = Config$parameter$N.Timepoints - 1,
# #                                                            Timepoints.Patients = Config$parameter$N.Timepoints - 1,
# #                                                            CTP.Vars            = Config$Varnames$OldCTPs[, -1])
#
# Config$CogDat$missForest.prepo.CTPs.grouped = grouping.FUN(Input = Config$CogDat$missForest.prepo.CTPs.grouped.postop)
#
# Config$CogDat$missForest.prepo.CTPs.grouped$Patients.T0 = Config$CogDat$missForest.prepo.CTPs$Patients.preop$ximp
#
# Config$CogDat$missForest.prepo.CTPs.grouped$Controls.T0 = Config$CogDat$missForest.prepo.CTPs$Controls.preop$ximp
#
#
# # remove all patients at that have no data
# buffer = list()
#
# for (i in 1:length(Config$CogDat$missForest.prepo.CTPs.grouped)) {
#   buffer[[i]] = remove.na.2.FUN(Config$CogDat$missForest.prepo.CTPs.grouped[[i]], Config$CogDat$clean.CTPs[[i]])
# }
#
# names(buffer) = names(Config$CogDat$missForest.prepo.CTPs.grouped)
#
#
# # fill up follow up measures with empty rows from baseline (necessary for some calculations)
# Config$CogDat$clean.missForest.prepo.CTPs.grouped = lapply(buffer, id.FUN)
#
# # patients
# for (i in 2:Config$parameter$N.Timepoints) {
#   Config$CogDat$clean.missForest.prepo.CTPs.grouped[[i]] = merge(Config$CogDat$clean.missForest.prepo.CTPs.grouped[[i]],
#                                                                  Config$CogDat$clean.missForest.prepo.CTPs.grouped$Patients.T0,
#                                                                  by = "SubjectID",
#                                                                  all.y = TRUE)
#
#   # make SubjectID the rowname again
#   rownames(Config$CogDat$clean.missForest.prepo.CTPs.grouped[[i]]) =
#     Config$CogDat$clean.missForest.prepo.CTPs.grouped[[i]]$SubjectID
#
#   # remove unneeded columns
#   Config$CogDat$clean.missForest.prepo.CTPs.grouped[[i]] =
#     Config$CogDat$clean.missForest.prepo.CTPs.grouped[[i]][, c(2:(Config$N$CTPs + 1))]
# }
#
# # controls
# for (i in (Config$parameter$N.Timepoints + 2):(Config$parameter$N.Timepoints + Config$parameter$N.Timepoints)) {
#   Config$CogDat$clean.missForest.prepo.CTPs.grouped[[i]] = merge(Config$CogDat$clean.missForest.prepo.CTPs.grouped[[i]],
#                                                                  Config$CogDat$clean.missForest.prepo.CTPs.grouped$Controls.T0,
#                                                                  by = "SubjectID",
#                                                                  all.y = TRUE)
#
#   # make SubjectID the rowname again
#   rownames(Config$CogDat$clean.missForest.prepo.CTPs.grouped[[i]]) =
#     Config$CogDat$clean.missForest.prepo.CTPs.grouped[[i]]$SubjectID
#
#   # remove unneeded columns
#   Config$CogDat$clean.missForest.prepo.CTPs.grouped[[i]] =
#     Config$CogDat$clean.missForest.prepo.CTPs.grouped[[i]][, c(2:(Config$N$CTPs + 1))]
# }
#
# # keep only columns with CTPs
# Config$CogDat$clean.missForest.prepo.CTPs.grouped[[1]] =
#   Config$CogDat$clean.missForest.prepo.CTPs.grouped[[1]][, c(1:Config$N$CTPs)]
# Config$CogDat$clean.missForest.prepo.CTPs.grouped[[Config$parameter$N.Timepoints+1]] =
#   Config$CogDat$clean.missForest.prepo.CTPs.grouped[[Config$parameter$N.Timepoints+1]][, c(1:Config$N$CTPs)]
#
# rm(buffer, i)
# ```
