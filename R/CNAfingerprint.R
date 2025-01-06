
#' @title CNAfingerprint prediction function
#'
#' @param cnf CNA features
#' @param target Model Prediction Targets, only support "Oxa" now
#'
#' @return CNAfingerprint prediction result
#' @importFrom dplyr mutate
#' @importFrom dplyr rename_with
#' @importFrom dplyr case_when
#' @export
#'
#' @examples
#' \donttest{
#' features <- readRDS(system.file("extdata", "examplefeatures.rds",package = "CNAfingerprint", mustWork = TRUE))
#' features$AFP=c(5.6,18)
#' features$`CA19-9` = c(5.6,19)
#' score <- CNAfingerprint(features,target="Oxa")
#' }
CNAfingerprint <- function(cnf,target="Oxa"){

  if(target=="Oxa"){
    cat("You are running CNAfingerprint to predict Oxaliplatin-based response in mCRCs")
    path <- system.file("extdata", "xgb_last.rds", package = "CNAfingerprint", mustWork = TRUE)
    xgb_last <- readRDS(path)

    target_cnf <- c("CN[>8]","E:LL:9+:BB","AFP","L:LL:9+:BB","CNCP[3]","cna_burden",
                    "CA19-9","CN[2]")
    valdata <- cnf[,target_cnf]

    # feature name
    data <- valdata %>%
      mutate(across(everything(), as.numeric)) %>%
      rename_with(~ case_when(
        .x == "CN[>8]" ~ "feature84",
        .x == "E:LL:9+:BB" ~ "feature52",
        .x == "AFP" ~ "feature7",
        .x == "L:LL:9+:BB" ~ "feature42",
        .x == "CNCP[3]" ~ "feature74",
        .x == "cna_burden" ~ "feature15",
        .x == "CA19-9" ~ "feature5",
        .x == "CN[2]" ~ "feature81",
        TRUE ~ .x
      ))
    #
    set.seed(123)
    pred <- stats::predict(xgb_last, as.matrix(data))
    valdata$CNAfingerprint <- pred
  }else{
    message("Currently, CNAfingerprint is only used to predict response to oxaliplatin chemotherapy")
    message("Stay tuned for other applications!")
    stop("please setting target='Oxa'")
  }

  return(valdata)
}
