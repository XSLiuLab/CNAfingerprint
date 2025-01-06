#' @title CNA features calling function
#'
#' @param data data.frame with colnames: "chromosome", "start", "end", "segVal","sample"
#' @param hg the ref genome, support hg38 only.
#'
#' @return CNA features
#'
#'
#' @export
#'
#' @examples \donttest{
#' exampleSeg <- readRDS(system.file("extdata", "exampleSeg.rds",package = "CNAfingerprint", mustWork = TRUE))
#' features <- CNF_call(exampleSeg,hg="hg38")
#' }
#'

CNF_call <- function(data,hg = "hg38"){
  #library(sigminer)
  coln <- c("chromosome", "start", "end", "segVal","sample")
  #
  if(hg != "hg38"){
    stop("Error: The genome must be 'hg38', use function hg19to38 to transform.")
  }
  if(!all(coln %in% colnames(data))){
    stop("Error:The colnames should be: chromosome, start, end, segVal,sample")
  }

  CN <- data[,coln]

  # Wang et al.
  cn <- sigminer::read_copynumber(CN,
                        seg_cols = c("chromosome", "start", "end", "segVal"),
                        samp_col = "sample",
                        genome_build = "hg38", complement = FALSE)
  #
  cn_sum <- cn@summary.per.sample
  rownames(cn_sum) <- cn_sum$sample
  tally_w <- sigminer::sig_tally(cn, method = "W", cores = 10)
  # if(!is.null(note)){
  #   saveRDS(tally_w,file = paste0("./data/results/",note,"_tally_W.rds"))
  # }

  cn_fet <- tally_w[["nmf_matrix"]]
  feature_m <- data.frame(cn_fet)
  colnames(feature_m) <- colnames(cn_fet)
  feature_m$sample <- rownames(cn_fet)
  # if only one sample,rownames(cn_fet)will be null
  colon_cnf <- merge(cn_sum,feature_m,by = "sample", all = T)

  # Yao et al.
  #HRD(cutoff 0.2)
  #library(HRDCNA)
  # score_yhz <- HRDCNA::HRDprediction(data = feature_m)
  # colon_cnf <- merge(colon_cnf,score_yhz,by = "sample", all = T)

  # Tao et al.
  tally_X_noLOH <- sigminer::sig_tally(cn,
                             method = "X",
                             add_loh = FALSE,
                             cores = 10
  )
  # if(!is.null(note)){
  #   saveRDS(tally_X_noLOH, file = paste0("./data/results/",note,"_tally_X_noLOH.rds"))
  # }

  cn_fet <- tally_X_noLOH[["all_matrices"]][["simplified_matrix"]]
  feature_m <- data.frame(cn_fet)
  colnames(feature_m) <- colnames(cn_fet)
  feature_m$sample <- rownames(cn_fet)

  colon_cnf <- merge(colon_cnf,feature_m,by = "sample", all = T)

  # global CNA status(amp and del states of every chr)
  gCNA <- data.frame(matrix(nrow = 0,ncol = 44))
  coln <- c()
  for(i in 1:22) {
    coln <- c(coln, paste0("chr[", i, "]AMP"), paste0("chr[", i, "]DEL"))
  }
  chrname <- paste0("chr",1:22)

  for(i in unique(data$sample)){
    df <- data[data$sample==i,]
    onechr <- c(i)
    for(j in chrname){
      df1<- df[df$chromosome==j,]
      if(any(df1$segVal >2)){
        amp=1
      }else{
        amp=0
      }

      if(any(df1$segVal <2)){
        del=1
      }else{
        del=0
      }
      onechr <- c(onechr,amp,del)
    }
    onechr <- data.frame(t(onechr))
    colnames(onechr) <- c("sample",coln)

    gCNA <- rbind(gCNA,onechr)
  }
  colon_cnf <- merge(colon_cnf,gCNA,by="sample")
  return(colon_cnf)
}
