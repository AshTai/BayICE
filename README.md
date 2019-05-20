# BayICE
BayICE: A hierarchical Bayesian deconvolution model with stochastic search variable selection.

To install this package, please copy and paste following codes into your R session:

1. install.packages("devtools")
2. library(devtools)
3. install_github("AshTai/BayICE")

## Example
```R
library(BayICE)
data("NSCLC")
# Get the gene-specific expression for generating simulated data from real cancer data
n.base <- rowMeans(mdata.u[,grep("N",colnames(mdata.u))])
n.base <- n.base[n.base>1 & n.base < 1000]

sim.data <- MN_sim(n.base=n.base)

# For sequencing data, we recommend the log-transformation of raw data as the BayICE input.
res <- BayICE_iter(ref.set=log(sim.data$ref+1), mix.set=log(sim.data$mix+1),ref.id=sim.data$type,iter=2000)
```

## Contact information
Anshun Tai ([daansh13@gmail.com](mailto:daansh13@gmail.com))
https://anshuntai.weebly.com
