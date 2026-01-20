
#' @title CNAfingerprint prediction function
#'
#' @param cnf CNA features
#' @param target Model Prediction Targets, only support "Oxa" now
#'
#' @return CNAfingerprint prediction result
#' @importFrom dplyr mutate
#' @importFrom dplyr rename_with
#' @importFrom dplyr case_when
#' @importFrom mlr3 as_task_classif
#' @export
#'
#' @examples
#' \donttest{
#' features <- readRDS(system.file("extdata", "examplefeatures.rds",package = "CNAfingerprint", mustWork = TRUE))
#' score <- CNAfingerprint(features,target="Oxa")
#' }
CNAfingerprint <- function(cnf,target="OXA"){
  library(mlr3verse)
  library(dplyr)

  if(target=="OXA"){
    cat("You are running CNAfingerprint to predict Oxaliplatin-based response for mCRCs")

    library(mlr3verse)

    path <- system.file("extdata", "xgboostle.rds", package = "CNAfingerprint", mustWork = TRUE)
    xgboostle <- readRDS(path)

    target_cnf <- c("sample","CN[>8]","cna_burden","chr[8]AMPlevel","chr[20]AMPlevel",
                    "chr[11]AMPlevel","E:LL:9+:BB","chr[5]AMPlevel")
    valdata <- cnf %>%
      as.data.frame() %>%
      dplyr::select(all_of(target_cnf)) %>%
      rename_with(~ case_when(
        .x == "CN[>8]" ~ "feature34",
        .x == "cna_burden" ~ "feature6",
        .x == "chr[8]AMPlevel" ~ "feature253",
        .x == "chr[20]AMPlevel" ~ "feature301",
        .x == "chr[11]AMPlevel" ~ "feature265",
        .x == "E:LL:9+:BB" ~ "feature121",
        .x == "chr[5]AMPlevel" ~ "feature241",
        TRUE ~ .x
      ))
    #
    valdata$lable <- sample(c(0, 1), size = nrow(valdata), replace = TRUE)
    valdata$lable <- factor(valdata$lable,levels = c(0,1))
    task_val <- as_task_classif(valdata[,2:ncol(valdata)], target = "lable",positive = '1')

    set.seed(2025)
    pred <- xgboostle$predict(task_val)
    valdata$CNAfingerprint <- pred$prob[,1]
    valdata <- valdata %>% dplyr::select(all_of(c("sample","CNAfingerprint")))

  }else if(target=="BCG"){
    cat("You are running CNAfingerprint to predict clinical response of BCG perfusion therapy for NMIBC")

    library(mlr3verse)

    path <- system.file("extdata", "bcg_model.rds", package = "CNAfingerprint", mustWork = TRUE)
    xgboostle <- readRDS(path)

    target_cnf <-readRDS(system.file("extdata", "bcg_features.rds", package = "CNAfingerprint", mustWork = TRUE))

    valdata <- cnf %>% as.data.frame()
    colnames(valdata)[2:ncol(valdata)] <-paste0("feature",2:ncol(valdata))

    valdata$lable <- sample(c(0, 1), size = nrow(valdata), replace = TRUE)
    valdata$lable <- factor(valdata$lable,levels = c(0,1))
    task_val <- as_task_classif(valdata[,2:ncol(valdata)], target = "lable",positive = '1')

    set.seed(2025)
    pred <- xgboostle$predict(task_val$select(target_cnf))
    valdata$CNAfingerprint <- pred$prob[,1]
    valdata <- valdata %>% dplyr::select(all_of(c("sample","CNAfingerprint")))

  }else if(target=="BCR"){
    cat("You are running CNAfingerprint to predict biochemical recurrence risk in prostate cancer")

    library(mlr3verse)

    path <- system.file("extdata", "BCRranger_2026_6.rds", package = "CNAfingerprint", mustWork = TRUE)
    xgboostle <- readRDS(path)

    target_cnf <-c("chr[5+8]DELratio",
            "chr[10+18]DELratio",
            "BoChr[1+16]",
            "BoChr[6+8+12]",
            "BP10MB[3-4]",
            "ML:LD:4:AA",
            "M:LD:3:AA",
            "ML:LL:5-8:AA",
            "SS[5-6 & 7-8]")

    feature_BCR <- function(df,target_cnf){
      df$`chr[10+18]DELratio` <- df$`chr[10]DELratio`+df$`chr[18]DELratio`
      df$`chr[8]AMP+[18]DELratio` <- df$`chr[8]AMPratio`+df$`chr[18]DELratio`
      df$`chr[5+8]DELratio` <- df$`chr[8]DELratio`+df$`chr[5]DELratio`
      df$`BoChr[1+16]` <- df$`BoChr[1]`+df$`BoChr[16]`
      df$`BoChr[6+8+12]` <- df$`BoChr[6]`+df$`BoChr[8]`+df$`BoChr[12]`
      df$`SS[5-6 & 7-8]` <- df$`SS[>7 & <=8]` + df$`SS[>5 & <=6]`
      df$`ML:LL:5-8:AA` <-  df$`M:LL:5-8:AA`+ df$`L:LL:5-8:AA`
      df$`ML:LD:4:AA` <-  df$`M:LD:4:AA`+ df$`L:LD:4:AA`
      df$`BP10MB[3-4]` <- df$`BP10MB[3]`+df$`BP10MB[4]`

      df <- df[,c("sample",select)]
      colnames(df)[2:10] <- paste0("feature",1:9)
      df$Recurrence_event <- sample(c(0, 1), size = nrow(df), replace = TRUE)
      df$Recurrence_time <- 50
      return(df)
      }

    valdata <- cnf %>% as.data.frame()
    valdata<- feature_BCR(valdata)

    task_val <- mlr3proba::as_task_surv(valdata[,2:12], 
                                     time = "Recurrence_time",
                                     event = "Recurrence_event")

    set.seed(2025)
    pred <- xgboostle$predict(task_val)
    
    df <- data.frame(sample=valdata[,"sample"],
                     CNAfingerprint_crank=pred$crank,
                     t(pred$distr[1:length(task_val$row_ids)]$survival(c(12,24,36,48,60,120))))
    colnames(df)[3:8] <- paste0("SurProb",c(12,24,36,48,60,120),"M")
    valdata <- df

  }else{
    message("Currently, CNAfingerprint is only used to predict response to oxaliplatin chemotherapy")
    message("Stay for other applications!")
    stop("please setting target='OXA' or target='BCG'")
  }

  return(valdata)
}
