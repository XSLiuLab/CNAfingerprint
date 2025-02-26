
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

  }else{
    message("Currently, CNAfingerprint is only used to predict response to oxaliplatin chemotherapy")
    message("Stay for other applications!")
    stop("please setting target='Oxa'")
  }

  return(valdata)
}
