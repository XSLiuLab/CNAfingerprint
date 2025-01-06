# CNAfingerprint

CNAfingerprint, a tumor genomic DNA copy number alteration (CNA) based tool for predicting oxaliplatin-based chemotherapy clinical response.

## Installation 
Requirements:

R (>= 3.5.0);

dplyr;

sigminer;

xgboost;

carte.


```r
devtools::install_github("XSLiuLab/CNAfingerprint")

library(CNAfingerprint)
```

##Usegae


```r

exampleSeg <- readRDS(system.file("extdata", "exampleSeg.rds",package = "CNAfingerprint", mustWork = TRUE))
features <- CNF_call(exampleSeg,hg="hg38")

# 
features$AFP=c(5.6,18)
features$`CA19-9` = c(5.6,19)


score <- CNAfingerprint(features,target="Oxa")

```

### note

1,Be careful about sample ploidy selection.

2,Hg38 is the only supported version of the genome, you can convert hg19 to hg38
with the following methods:

```r
library(GenomicRanges)
library(liftOver)

cytoarm <- cytobandToArm(ucsc.hg38.cytoband)
ch <- import.chain('hg19ToHg38.over.chain')

exampleSeg <- readRDS(system.file("extdata", "exampleSeg.rds",package = "CNAfingerprint", mustWork = TRUE))

segdata <- exampleSeg[exampleSeg$sample=="sample2",]
seg <- GRanges(seqnames = segdata$chromosome,
  ranges = IRanges(start = segdata$start,end = segdata$end),
  strand = "*",
  TCN = segdata$segVal)
  
hg38.seg <- liftOver(seg, ch)

```
download chain file from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/


## Citation
Junyong Weng, Jinyu Wang, Ziyu Tao, et.al.Copy number alteration fingerprint predicts the clinical response of oxaliplatin-based chemotherapy in metastatic colorectal cancer.


## Contributors
CNAfingerprint was developed by Ziyu Tao and Jinyu Wang. Please contact Jinyu Wang: wangjy10@shanghaitech.edu.cn for any questions or suggestions. Thank you for your use and feedback.

---

Cancer Biology Group @ShanghaiTech

Research group led by Xue-Song Liu in ShanghaiTech University
