#' calculation of change scores
#'
#' @param Input list of data frames with poled CTPs
#' @param Timepoint.names names of the test timepoints
#' @param Groups Groups of Treatment used for the sublists e.g. Patients and Controls
#' @return list of data frames with calculated change scores
#' @export calc.changeScore.FUN

calc.changeScore.FUN <- function(Input,#  = Config$CogDat$CTP.Poled.missForest,
                                 Timepoint.names,# = Config$parameter$names.Timepoints,
                                 Groups)# = names(Config$parameter$names.Group))
  {
  CS <- Input
  for (i in Groups){
  df <- CS[grepl( i, names(CS))]

 CS$Delta[[i]] = lapply(
    df[grep(paste(Timepoint.names[-1], collapse = "|"), names(df))],
    function(x) {
      df[[grep(Timepoint.names[1], names(df))]] <- df[[grep(Timepoint.names[1], names(df))]][rownames(df[[grep(Timepoint.names[1], names(df))]]) %in% rownames(x),]
      magrittr::subtract(x, df[[grep(Timepoint.names[1], names(df))]])})


  }
  return(CS$Delta)
}

############################################################################

#' corrected change scores
#'
#' @param Input list of data frames with uncorrected change score of patients group
#' @param Control list of data frames with uncorrected change score of control group
#' @param Timepoint.names names of the test timepoints
#' @param CTP.names names of CTP tests
#' @return list of data frames with corrected change scores
#' @export correct.changeScore.FUN

correct.changeScore.FUN <- function(Input,#  = Config$CogDat$no_imp$CS.Delta$Patients,
                                    Control,# = Config$CogDat$no_imp$CS.Delta$Controls,
                                    Timepoint.names,# = Config$parameter$names.Timepoints,
                                    CTP.names)#        = Config$parameter$colnames.CTPs)
  {
  CS <- list()

  calc.CS.FUN <- function(Timepoint){

    CS[[Timepoint]] <- data.frame(sapply(grep(Timepoint, CTP.names, value = T), function(i) {

      magrittr::subtract( Input[[grep(Timepoint,names(Input))]][[i]],
                          mean(Control[[grep(Timepoint,names(Input))]][[i]], na.rm = T))}),
                        row.names = rownames(Input[[grep(Timepoint,names(Input))]]))

    return(CS)
  }

  CS <- lapply(Timepoint.names, calc.CS.FUN)
  CS <- unlist(CS, recursive = F)

  return(CS)
  }

############################################################################

#' calculate Z- scores
#'
#' @param Input.all list of data frames with corrected change score of patients group
#' @param Control.all list of data frames with corrected change score of control group
#' @param Timepoint.names names of the test timepoints
#' @param CTP.names names of CTP tests
#' @return list of data frames with calculated z scores
#' @export calc.Zscore.FUN

calc.Zscore.FUN <- function (Input.all,#   = Config$CogDat$CTP.Poled.missForest[grep("Patients.T0", names(Config$CogDat$CTP.Poled.missForest))],
                             Control.all,# = Config$CogDat$CTP.Poled.missForest[grep("Controls.T0", names(Config$CogDat$CTP.Poled.missForest))],
                             Timepoint.names,# = "T0",
                             CTP.names)#        = Config$parameter$colnames.CTPs)
  {
  Zscore = list()
  Zcalc.FUN <- function(Timepoint){

  Input   <- Input.all[[grep(Timepoint, names(Input.all))]]
  Control <- Control.all[[grep(Timepoint, names(Control.all))]]


  Zscore[[Timepoint]] <- data.frame(sapply(grep(Timepoint, CTP.names, value = T),
                                            function(i) { magrittr::subtract(Input[[i]], mean(Control[[i]], na.rm = T)) / stats::sd(Control[[i]], na.rm = T)}),
                          row.names = rownames(Input))

  return(Zscore)
  }

  Zscore <- lapply(Timepoint.names, Zcalc.FUN)
  Zscore <- unlist(Zscore, recursive = F)

  return(Zscore)

}

############################################################################
#' calculate combined Z- scores
#'
#' @param Input.all list of data frames containing z score for each CTP of patients group
#' @param Control.all list of data frames containing z score for each CTP of control group
#' @param Timepoint.names names of the test timepoints
#' @return list of data frames with combine scores
#' @export calc.combined.Zscore.FUN

calc.combined.Zscore.FUN <- function(Input.all,#  = Config$CogDat$imp$Zscore.Baseline.Patients,
                                     Control.all,# = Config$CogDat$imp$Zscore.Baseline.Control,
                                     Timepoint.names)# = "T0"
  {
  RCI <- list()

  Zcomb.FUN <- function (Timepoint) {
    Input   <- Input.all[[grep(Timepoint, names(Input.all))]]
    Control <- Control.all[[grep(Timepoint, names(Control.all))]]


  # take only subjects into account that have a complete testing
  RCI.Complete = stats::na.omit(Input)
  RCI.Complete.Control = stats::na.omit(Control)
  # calculate the sum of single RCI scores for every subject with complete testing
  RCI.Sum = rowSums(RCI.Complete)
  RCI.Sum.Control = rowSums(RCI.Complete.Control)
  # calculate SD of RCI sum for the control group
  RCI.SD = stats::sd(RCI.Sum.Control)

  # create template for new data frame
  RCI.Combined = RCI

  # calculate combined Z score using the standard deviation of the RCI sum of the controls. If the control group has less test time points than patients, use the one from the first follow up measure.
  RCI.Combined[[Timepoint]] = data.frame(Zscore.combined = sapply(RCI.Sum, "/", RCI.SD), row.names = rownames(RCI.Complete))
  return(RCI.Combined)
}
  RCI.Combined <- lapply(Timepoint.names, Zcomb.FUN)
  RCI.Combined <- unlist(RCI.Combined, recursive = F)
  return(RCI.Combined)
}

############################################################################
#' calculate POCD criteria 7, 8, 9 according to Rasmussen et al. (2001)
#' @description function for calculation POCD criteria 7, 8, 9 according to Rasmussen et al. (2001)
#' @param RCI.All list of data frames containing z scores for each CTP
#' @param RCI.Combined list of data frames containing combined z scores for each CTP
#' @param Timepoint.names names of the test timepoints
#' @param thresh threshhold; Z values below this threshhold are classified as "hit"; default is -1.96
#' @param hit how many "hits" are necessary for single POCD diagnosis; default is 2
#' @param Diagnose name of Diagnosis by a hit , e.g. PreCI or POCD
#' @return  List with data frames containing final POCD diagnosis
#' @export rci.diagnosis.FUN

# function for calculation POCD criteria 7, 8, 9 according to Rasmussen et al. (2001)
rci.diagnosis.FUN = function(RCI.All,# = Config$CogDat$imp$Zscore.Baseline.Patients,
                             thresh = -1.96,
                             hit = 2,
                             RCI.Combined,# = Config$CogDat$imp$Zscore.combined.Baseline.Patients,
                             Timepoint.names,# = "T0",
                             Diagnose = "PreCI")
  {

  POCD.Diagnosis.9 <- list()

  POCD.diagn.FUN <- function(Timepoint){

    RCI  <- RCI.All[[grep(Timepoint, names(RCI.All))]]
    RCI.comb <- RCI.Combined[[grep(Timepoint, names(RCI.Combined))]]

    # take only subjects into account that have a complete testing
  RCI.Complete = stats::na.omit(RCI)

  # check if Z score is < -1.96
  POCD.Complete =  RCI.Complete < thresh

  # POCD criterium 7: deterioration in 2 or more CTPs < -1.96
  POCD.Diagnosis.7 = data.frame(
    SubjectID = row.names(POCD.Complete),
    ZcoreHit  = rowSums(POCD.Complete) >= hit,
    CTPs_hit  = sapply( 1:nrow(POCD.Complete), function(x) { paste( names(which( POCD.Complete[x, 1 : ncol(POCD.Complete)]  == TRUE)), collapse = " - " )})

)

  # POCD criterium 8: combined Z score < -1.96
  POCD.Diagnosis.8 = data.frame(
    SubjectID = row.names(RCI.comb),
    Zscore.combined  = RCI.comb < thresh)

  # POCD criterium 9: either deterioration in 2 or more CTPs < -1.96 or combined Z score < -1.96
  POCD.Diagnosis.9[[Timepoint]] = plyr::join_all(list(POCD.Diagnosis.8, POCD.Diagnosis.7), by="SubjectID") %>% plyr::mutate(Diag = rowSums(.data[2:3])>=1) %>% plyr::rename(c("Diag" = Diagnose))
  colnames(POCD.Diagnosis.9[[Timepoint]])[-1] <- paste0(Timepoint,".", colnames(POCD.Diagnosis.9[[Timepoint]])[-1])
  colnames(POCD.Diagnosis.9[[Timepoint]][]) <- paste0(Timepoint,".", colnames(POCD.Diagnosis.9[[Timepoint]])[-1])
  return(POCD.Diagnosis.9)
  }

  POCD.Diagnosis.9 <- lapply(Timepoint.names,POCD.diagn.FUN)
  POCD.Diagnosis.9 <- unlist(POCD.Diagnosis.9, recursive = F)
  return(POCD.Diagnosis.9)

}

############################################################################

#' calculate NCD
#'
#' @param Input.all list fo data frames containing CTPs /RCIs of patients group
#' @param Control.all list fo data frames containing CTPs/RCIs of control group
#' @param Timepoint.names names of the test timepoints
#' @param CTP.names names of CTP tests
#' @param CTP.number numer of used CTP test
#' @param impairment list of data frames containing the information about the impairment of the patients in logicals (TRUE, FALSE NA)
#' @param NCD.names name of column containing the variables for impairment
#' @return list of data frames comtaining NCD diagnosis for each patient
#' @export calc.NCD.FUN

calc.NCD.FUN <- function (   Input.all,#   = Config$CogDat$CTP.Poled.missForest[grep("Patients", names(Config$CogDat$CTP.Poled.missForest))],
                             Control.all,# = Config$CogDat$CTP.Poled.missForest[grep("Controls", names(Config$CogDat$CTP.Poled.missForest))],
                             impairment,# = Config$CogDat$CTPs.control.clean.interim[grep("Patients", names(Config$CogDat$CTPs.control.clean.interim))],
                             Timepoint.names,# = "T0",
                             CTP.names,#        = Config$parameter$colnames.CTPs,
                             CTP.number,#       = Config$parameter$N.CTPs,
                             NCD.names)#        = Config$parameter$colnames.NCD
  {
  Score = list()
  NCDcalc.FUN <- function(Timepoint){

    Input   <- Input.all[[grep(Timepoint, names(Input.all))]]
    Input   <- Input[!is.na(rowSums(Input)),]
    Control <- Control.all[[grep(Timepoint, names(Control.all))]]
    Control <- Control[!is.na(rowSums(Control)),]
    imp     <- impairment[[grep(Timepoint, names(impairment))]]
    #imp     <- imp[[grep(NCD.names, colnames(imp))]]
    Input <- merge(x = Input, y = imp[, c("SubjectID", colnames(imp)[colnames(imp) %in% NCD.names])],
                    by.x = "row.names", by.y = "SubjectID",
                    all.x = T  )


    Score[[Timepoint]] <- data.frame(sapply(grep(Timepoint, CTP.names, value = T), function(i) {

      Score[[i]][Input[[i]] < magrittr::subtract(mean(Control[[i]], na.rm = T), stats::sd(Control[[i]], na.rm = T)) &
         Input[[i]] > magrittr::subtract(mean(Control[[i]], na.rm = T), 2*stats::sd(Control[[i]], na.rm = T))] <- 1
      Score[[i]][Input[[i]] <= magrittr::subtract(mean(Control[[i]], na.rm = T), 2*stats::sd(Control[[i]], na.rm = T))] <- 2
      Score[[i]][Input[[i]] >= magrittr::subtract(mean(Control[[i]], na.rm = T), stats::sd(Control[[i]], na.rm = T))] <- 0

        return(Score[[i]])
        }), row.names = Input$Row.names)

    Score[[Timepoint]]$score <- apply(Score[[Timepoint]][,1: CTP.number], 1, max)
    Score[[Timepoint]]$imp   <- Input[[grep("Impairment", colnames(Input))]]

        #Score[[Timepoint]] <- data.frame( score = apply(Score[[Timepoint]][,1: CTP.number], 1, max),
        #                                  impairment = Score[[Timepoint]]$imp)

    Score[[Timepoint]]$NCD[Score[[Timepoint]]$score == 2 & Score[[Timepoint]]$imp == TRUE] <- "major"

    Score[[Timepoint]]$NCD[Score[[Timepoint]]$score == 2 & Score[[Timepoint]]$imp == FALSE |
                           Score[[Timepoint]]$score == 1 & Score[[Timepoint]]$imp == TRUE  |
                           Score[[Timepoint]]$score == 1 & is.na(Score[[Timepoint]]$imp)   |
                           Score[[Timepoint]]$score == 1 & Score[[Timepoint]]$imp == FALSE ] <- "mild"

    Score[[Timepoint]]$NCD[Score[[Timepoint]]$score == 0 & Score[[Timepoint]]$imp == TRUE |
                           Score[[Timepoint]]$score == 0 & is.na(Score[[Timepoint]]$imp)    |
                           Score[[Timepoint]]$score == 0 & Score[[Timepoint]]$imp == FALSE  ] <-   FALSE

    Score[[Timepoint]]$NCD[Score[[Timepoint]]$score == 2 & is.na(Score[[Timepoint]]$imp)] <- NA


colnames(Score[[Timepoint]])[which(colnames(Score[[Timepoint]]) %in% c("score","imp","NCD"))] <- c(paste0(Timepoint, c("_score","_impair","_NCD" )))
Score[[Timepoint]]$SubjectID <- row.names(Score[[Timepoint]])
    return(Score)
  }

  NCDscore <- lapply(Timepoint.names, NCDcalc.FUN)
  NCDscore <- unlist(NCDscore, recursive = F)

    return(NCDscore)

}
