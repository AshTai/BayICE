% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayICE.R
\name{BayICE}
\alias{BayICE}
\title{BayICE: A hierarchical Bayesian deconvolution model with stochastic search variable selection.}
\usage{
BayICE(ref.set, mix.set, ref.id, iter = 10^4, parC = 10^-5, burnin = 0.6,
  thinning = 3)
}
\arguments{
\item{mix.set}{A G by M expression matrix of bulk samples, where G is the gene number and M is the total number of bulk samples.}

\item{ref.id}{A vector of cell types.}

\item{iter}{The number of iterations. Default is 10000.}

\item{parC}{A parameter of the spike-and-slab distribution. Default is 10^-5.}

\item{res.set}{A G by N expression matrix of reference, where G is the gene number and N is the total number of reference.}
}
\value{
A list of estimated cell proportions (w),
the probability of informative genes (p),
the fraction of informative genes per cell type (pi), and
the cell-specific profiles (mean).
}
\description{
The main function of BayICE.
}
\examples{
data("NSCLC")
n.base <- rowMeans(mdata.u[,grep("N",colnames(mdata.u))])
n.base <- n.base[n.base>1 & n.base < 1000]
sim.data <- MN_sim(n.base=n.base)
res <- BayICE(ref.set=log(sim.data$ref+1), mix.set=log(sim.data$mix+1),ref.id=sim.data$type,iter=2000)
}
\author{
An-Shun Tai \email{daansh13@gmail.com}
}
