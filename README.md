# BayICE
BayICE: A hierarchical Bayesian deconvolution model with stochastic search variable selection.

To install this package, please copy and paste following codes into your R session:

1. install.packages("devtools")
2. library(devtools)
3. install_github("AshTai/BayICE")

## Example
```R
library(CloneDeMix)
data("ESCC_chr1")
# ESCC_chr1$tumor is a G by N count matrix of tumor samples, where G is the number of loci and N is the sample size.
# ESCC_chr1$normal is a DNA count vector of a normal sample. 

res <- CloneDeMix(tumor=ESCC_chr1$tumor, normal=ESCC_chr1$normal,threshold = 10^-5, iterC = 10^3,
  CNVstate = c(0:10), method = "aic")
head(res$CNV); head(res$MCP)
```

## Contact information
Anshun Tai ([daansh13@gmail.com](mailto:daansh13@gmail.com))
https://anshuntai.weebly.com
