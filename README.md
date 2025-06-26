# CNAfingerprint

Copy number alterations (CNAs) are a predominant source of genetic alterations in human cancer and play an important role in cancer progression. 

In the previous work ([Wang et al.](https://doi.org/10.1371/journal.pgen.1009557) and [Tao et al.](https://doi.org/10.1093/bib/bbad053) ), 
we developed a mechanism-agnostic method to categorize CNA based on various fragment properties, 
which reflect the consequences of mutagenic processes and can be extracted from  different types of data (WGS, SNP array), also the low-cost shallow WGS, 
this allowed us to get information about the tumor at a lower price!

So we developed CNAfingerprint, a tumor genomic DNA copy number alteration (CNA) based tool for assisting in clinical diagnosis and medication guidance of tumors.


As of today, CNAfingerprint can be implemented as follows:

1, target="OXA", clinical response of Oxaliplatin-based chemotherapy in metastatic colorectal cancer (mCRC).

2, target="BCG", clinical response of bacillus Calmette-Guérin (BCG) perfusion therapy in non-muscle invasive bladder cancer (NMIBC).

3. ...

we also developed a pan cancer homologous recombination deficiency (HRD) predictor, find it [here](https://github.com/XSLiuLab/HRDCNA)!


## Installation 

Requirements:

R (>= 3.5.0);

dplyr;

sigminer;

xgboost;

mlr3

carte.


```r
devtools::install_github("XSLiuLab/CNAfingerprint")

library(CNAfingerprint)
```

##Usegae


```r

exampleSeg <- readRDS(system.file("extdata", "exampleSeg.rds",package = "CNAfingerprint", mustWork = TRUE))
features <- CNF_call(exampleSeg,hg="hg38")

# clinical response of Oxaliplatin-based chemotherapy in mCRC

score <- CNAfingerprint(features,target="OXA")

# clinical response of Oxaliplatin-based chemotherapy in NMIBC

score <- CNAfingerprint(features,target="BCG")


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

Yao, H. et al. Copy number alteration features in pan-cancer homologous recombination deficiency prediction and biology. Commun Biol 6, 527 (2023). https://doi.org/10.1038/s42003-023-04901-3


Jy Wang. et.al. Copy number alteration fingerprint predicts the clinical response of oxaliplatin-based chemotherapy in metastatic colorectal cancer. [Under review]


Tt Cai. et.al. Copy number alteration feature predicts the clinical response of BCG perfusion therapy in non-muscle invasive bladder cancer. [Under review]

## Contributors

CNAfingerprint was developed by Zy Tao, Jy Shen and Jy Wang. Please contact Jy Wang: wangjy10@shanghaitech.edu.cn for any questions or suggestions. Thank you for your use and feedback.

---

Cancer Biology Group @ShanghaiTech

Research group led by Xue-Song Liu in ShanghaiTech University
