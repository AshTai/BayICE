#' BayICE: A hierarchical Bayesian deconvolution model with stochastic search variable selection.
#'
#' Multinomial-based simulator
#'
#' @name MN_sim
#' @author An-Shun Tai \email{daansh13@gmail.com}
#' @param seed seed. Default=123.
#' @param effect_size The normalized effect size. Default=0.5.
#' @param n.base The gene-specific baseline expression.
#' @return The expression matrix of mixed samples and reference. A vector of cell types.
#' @export
#' @examples
#' data("NSCLC")
#' n.base <- rowMeans(mdata.u[,grep("N",colnames(mdata.u))])
#' n.base <- n.base[n.base>1 & n.base < 1000]
#' sim.data <- MN_sim(n.base=n.base)

MN_sim <- function(seed=123,effect_size=0.5,n.base){
  require(MCMCpack)
  set.seed(seed=seed)
  sim.typeN <- 4; sim.refN.perT <- 20
  sim.id.type <- rep(1:sim.typeN,each=sim.refN.perT)
  sim.gN <- 5000; sim.refN <- sim.typeN*sim.refN.perT; sim.mixN <- 90
  sim.fgN <- 300

  sim_DM <- matrix(0,sim.typeN,sim.refN); sim_DM[cbind(sim.id.type,1:sim.refN)] <- 1
  sim_w <- t(rdirichlet(sim.mixN,alpha=rep(1,sim.typeN))*rep((1-seq(0.1,0.9,0.1)),each=10))
  sim_w <- rbind(sim_w,rep((seq(0.1,0.9,0.1)),each=10))


  sampling.alpha <- sample(1:length(n.base),sim.gN)
  sim_alpha <- n.base[sampling.alpha]

  ###################
  # S2
  effect_size <- effect_size
  sim_beta <- matrix(0,sim.gN,sim.typeN+1)
  sim_type.gene.index <- rep(0,sim.gN)
  sim_type.gene.index[1:sim.fgN] <- sample(1:(sim.typeN+1),sim.fgN,replace=T)
  for(i in 1:(sim.typeN+1)){
    sim_beta[sim_type.gene.index==i,i] <- sample(c(effect_size,-effect_size),size = sum(sim_type.gene.index==i),replace=T)*
      sim_alpha[sim_type.gene.index==i]
  }
  sim_beta[,5] <- sim_beta[,5] + sim_alpha*rnorm(sim.gN,0,0.03)


  noise.factor.mix <- matrix(runif(sim.gN*sim.mixN,0.8,1.2),sim.gN,sim.mixN)
  noise.factor.ref <- matrix(runif(sim.gN*sim.refN,0.8,1.2),sim.gN,sim.refN)

  sim_mu.mix <- (sim_alpha%*%matrix(1,1,sim.mixN))*noise.factor.mix+
    sim_beta%*%sim_w
  sim_mu.ref <- (sim_alpha%*%matrix(1,1,sim.refN))*noise.factor.ref+
    sim_beta[,-(sim.typeN+1)]%*%sim_DM
  sim_obs.mix <- apply(matrix(1:sim.mixN),1,
                       function(x)rmultinom(1,size=5*10^6,prob = sim_mu.mix[,x]))
  sim_obs.ref <- apply(matrix(1:sim.refN),1,
                       function(x)rmultinom(1,size=5*10^6,prob = sim_mu.ref[,x]))
  structure(list(mix=sim_obs.mix,ref=sim_obs.ref,type=sim.id.type))
}
