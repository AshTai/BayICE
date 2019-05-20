#' BayICE: A hierarchical Bayesian deconvolution model with stochastic search variable selection.
#'
#' The main function of BayICE.
#'
#' @name BayICE
#' @author An-Shun Tai \email{daansh13@gmail.com}
#' @param res.set A G by N expression matrix of reference, where G is the gene number and N is the total number of reference.
#' @param mix.set A G by M expression matrix of bulk samples, where G is the gene number and M is the total number of bulk samples.
#' @param ref.id A vector of cell types.
#' @param iter The number of iterations. Default is 10000.
#' @param parC A parameter of the spike-and-slab distribution. Default is 10^-5.
#' @return A list is containing cell proportions (w), gene selection (t), pi, and the cell-specific profiles.
#' @export
#' @examples
#' data("NSCLC")
#' n.base <- rowMeans(mdata.u[,grep("N",colnames(mdata.u))])
#' n.base <- n.base[n.base>1 & n.base < 1000]
#' sim.data <- MN_sim(n.base=n.base)
#' res <- BayICE_iter(ref.set=log(sim.data$ref+1), mix.set=log(sim.data$mix+1),ref.id=sim.data$type,iter=2000)

BayICE_iter <- function(ref.set,mix.set,ref.id,iter=10^4,parC=10^-5){
  require(MCMCpack)
  require(invgamma)
  require(quadprog)
  require(MASS)
  require(tibble)
  require(mvtnorm)
#  require(bigmemory)
#  options(bigmemory.typecast.warning=FALSE)

  if(nrow(ref.set)!=nrow(mix.set)){
    stop("the gene lists of reference set and mixed set are inconsistent")
  }else{n.gene <- nrow(ref.set)}
  id.numb <- 1
  id.type <- rep(NA,length(ref.id))
  for(id in unique(ref.id)){
    id.type[ref.id==id] <- id.numb
    id.numb <- id.numb+1
  }

  n.type <- length(unique(id.type));
  n.mix <- ncol(mix.set);
  n.ref <- ncol(ref.set);
  N.total <- n.mix+n.ref;
  type.design <- matrix(0,n.type,n.ref); type.design[cbind(id.type,1:n.ref)] <- 1
  uniq.id.type <- as.numeric(unique(id.type))

  ff.mean.ini0 <- function(x){
    if(sum(id.type==x)==1){
      return(ref.set[,id.type==x])
    }else{
      return(rowMeans(ref.set[,id.type==x]))
    }
  }
  mean.ini0 <- apply(as.matrix(uniq.id.type),1,ff.mean.ini0)
  N.sigma2 <- rowSums((ref.set-mean.ini0%*%type.design)^2)/(n.ref-n.type)
  ref <- ref.set*sqrt(N.total/N.sigma2)
  mix <- mix.set*sqrt(N.total/N.sigma2)

  #note: set data storage
  prior.t <- prior.mean <- array(NA,dim = c(n.gene,(n.type+1),iter))
  prior.pi <- matrix(NA,iter,(n.type+1))
  prior.w <- array(dim=c(n.type+1,n.mix,iter))

  #note: set initial value
  hyper.par <- list(c=parC,a0=5,b0=50,c0=1,d0=1)

  ff.mean.ini <- function(x){
    if(sum(id.type==x)==1){
      return(ref[,id.type==x])
    }else{
      return(rowMeans(ref[,id.type==x]))
    }
  }
  mean.ini <- apply(as.matrix(uniq.id.type),1,ff.mean.ini)

  #note: identify the initial signature set
  p.v<- matrix(NA,nrow(ref),n.type)
  for(i in 1:n.type){
    cell.v <- rep("A",n.ref); cell.v[id.type==uniq.id.type[i]] <- "B"
    if(sum(id.type==uniq.id.type[i])==1){
      p.v[,i] <- apply(ref,1,function(x)t.test(x[cell.v=="A"],mu=x[cell.v=="B"])$p.value)
    }else{
      p.v[,i] <- apply(ref,1,function(x)t.test(x[cell.v=="A"],x[cell.v=="B"])$p.value)
    }
  }

  t.ini <- matrix(0,n.gene,n.type+1)
  for(i in 1:n.type){
    t.ini[order(p.v[,i])[1:200],i] <- 1
  }

  #note: set the initial proportions
  nL.w <- function(x,y,z){
    w.m <- x
    value.nL <- sum((y-z%*%matrix(w.m))^2)
    value.nL
  }
  optim.w <- function(y,z,eta){
    ff <- constrOptim(rep(1/(n.type+1),n.type), f=nL.w,grad=NULL,y=y,z=z,ui = rbind(-rep(1,n.type), diag(n.type)), ci = c(-1,rep(0,n.type)))
    opt.w <- ff$par
    c(opt.w,1-sum(opt.w))
  }
  res.lm.const <- apply(mix[rowSums(t.ini)!=0,],2,optim.w,z=mean.ini[rowSums(t.ini)!=0,])
  res.lm.const[n.type+1,] <- 0.2
  res.lm.const <- apply(res.lm.const,2,function(x)x/sum(x))

  beta.miss.ini <- rowMeans( t(t(mix-mean.ini%*%res.lm.const[-c(n.type+1),])/res.lm.const[c(n.type+1),])) - rowMeans(as.matrix(ref[,id.type==1]))
  beta.ini <- (mean.ini)-rowMeans(as.matrix(ref[,id.type==1]))

  mean.miss <- t(t((mix-mean.ini%*%res.lm.const[-c(n.type+1),]))/res.lm.const[c(n.type+1),])
  p.v.miss <- apply(matrix(1:n.gene),1,function(x)t.test(ref[x,],mean.miss[x,])$p.value)
  t.ini[order(p.v.miss)[1:200],n.type+1] <- 1

  # set initial variance
  # mark #
  var.ref.ini <- rowMeans((cbind(ref-mean.ini%*%type.design))^2)
  var.mix.ini <- apply(mix,1,var)


  prior <- list(
    alpha=rowMeans(as.matrix(ref[,id.type==1])),
    beta=beta.ini,
    beta.miss=beta.miss.ini,
    psi2=(cbind(beta.ini,beta.miss.ini)+10^-3)^2,
    t=t.ini,
    w=res.lm.const,
    pi=rep(0.3,n.type+1),
    sigma2.ref=var.ref.ini/N.total,
    sigma2.mix=var.mix.ini/N.total
  )

  gene_weight <- rowSums(t.ini)!=0
  iter.c <- 1; stop.sign <- F; count <- 1
  run.time <- Sys.time()
  while(stop.sign==F){
    if(count==floor(iter.c/iter*100)){
      tt <- difftime(Sys.time(),run.time,units ="mins")
      print(paste(count,"%","--","Expected running time =",round((100-count)*tt,3),"mins"))  ;
      count <- count +1
      run.time <- Sys.time()
    }

    # update alpha
    #alpha_mean <- ( rowSums(ref-cbind(prior$beta)%*%type.design)+rowSums(mix-cbind(prior$beta,prior$beta.miss)%*%prior$w) )/(n.ref+n.mix)
    alpha_mean <- rowMeans(ref[,id.type==1])
    prior$alpha <- rnorm(n.gene,alpha_mean,sqrt(prior$sigma2.ref*N.total/(n.ref+n.mix)))

    # update beta
    beta.sampling <- function(d){
      w1.beta <- (prior$sigma2.ref*N.total/sum(id.type==d))^-1
      w2.beta <- (prior$sigma2.mix*N.total/sum(prior$w[d,]^2))^-1
      w3.beta <- ( ( (prior$t[,d]+(1-prior$t[,d])*hyper.par$c)*(prior$psi2[,d]) )^-1)*((prior$sigma2.mix*N.total)^-1)
      pos.mean.beta <-
        w1.beta/(w1.beta+w2.beta+w3.beta)*
        rowMeans( as.matrix( (ref-prior$alpha%*%matrix(1,1,n.ref))[,id.type==d] ) )+
        w2.beta/(w1.beta+w2.beta+w3.beta)*
        (mix-prior$alpha%*%matrix(1,1,n.mix)-cbind(prior$beta,prior$beta.miss)[,-d]%*%prior$w[-d,])%*%matrix(prior$w[d,])/sum(prior$w[d,]^2)
      pos.var.beta <- 1/(w1.beta+w2.beta+w3.beta)
      rnorm(n.gene,pos.mean.beta,sqrt(pos.var.beta))
    }
    prior$beta <- apply(as.matrix(uniq.id.type),1,beta.sampling)
    prior$beta[,1] <- 0

    # update beta.miss
    beta.miss.sampling <- function(d){
      w2.beta <- (prior$sigma2.mix*N.total/sum(prior$w[n.type+1,]^2))^-1
      w3.beta <- ( ( (prior$t[,(n.type+1)]+(1-prior$t[,(n.type+1)])*hyper.par$c)*
                       (prior$psi2[,(n.type+1)]) )^-1)*((prior$sigma2.mix*N.total)^-1)
      pos.mean.beta <-
        w2.beta/(w2.beta+w3.beta)*
        (mix-prior$alpha%*%matrix(1,1,n.mix)-cbind(prior$beta,prior$beta.miss)[,-(n.type+1)]%*%prior$w[-(n.type+1),])%*%matrix(prior$w[(n.type+1),])/sum(prior$w[(n.type+1),]^2)
      pos.var.beta <- 1/(w2.beta+w3.beta)
      rnorm(n.gene,pos.mean.beta,sqrt(pos.var.beta))
    }
    prior$beta.miss <- beta.miss.sampling(1)

    effect_size <- cbind(prior$beta,prior$beta.miss)
    prior.mean[,,iter.c] <- effect_size+prior$alpha

    # update sigma2.ref
    rss.ref <- rowSums( (ref-prior$alpha%*%matrix(1,1,n.ref)-cbind(prior$beta)%*%type.design)^2)/N.total
    prior$sigma2.ref <- rinvgamma(n.gene,(n.ref)/2-1,(rss.ref)/2)

    # update sigma2
    rss.mix <- rowSums( (mix-prior$alpha%*%matrix(1,1,n.mix)-cbind(prior$beta,prior$beta.miss)%*%prior$w)^2 )/N.total
    rss.beta <- rowSums( effect_size^2/(prior$t+hyper.par$c*(1-prior$t))/prior$psi2 )/N.total
    prior$sigma2.mix <- rinvgamma(n.gene,(n.mix+n.type+1)/2-1,(rss.mix+rss.beta)/2)


    # sampling t
    density.beta1 <- dnorm(effect_size,0,sqrt(prior$psi2*((prior$sigma2.mix*N.total)%*%matrix(1,1,n.type+1)) ),log = T) +
      log(matrix(1,n.gene,1)%*%prior$pi)
    density.beta2 <-  dnorm(effect_size,0,sqrt(prior$psi2*((prior$sigma2.mix*N.total)%*%matrix(1,1,n.type+1))*hyper.par$c),log = T) +
      log(matrix(1,n.gene,1)%*%(1-prior$pi))
    prob.beta <-  1/(1+exp(density.beta2-density.beta1))
    prior$t <-  matrix(0,n.gene,n.type+1)
    t.sample <- matrix(rbinom(n=n.gene*(n.type+1),1,prob=prob.beta),n.gene,(n.type+1))
    prior$t[t.sample==1] <- 1
    prior$t[,1] <- 0
    prior.t[,,iter.c] <- prior$t

    # update psi2
    par.invG.psi2 <- (effect_size)^2/(prior$t+hyper.par$c*(1-prior$t))/prior$sigma2.mix%*%matrix(1,1,n.type+1)/N.total
    prior$psi2 <- matrix(rinvgamma(n.gene*(n.type+1),hyper.par$a0+(1)/2,hyper.par$b0+par.invG.psi2),
                         n.gene,n.type+1)

    # update pi
    prior$pi <- rbeta((n.type+1),hyper.par$c0+colSums(prior$t==1),
                      hyper.par$d0+n.gene-colSums(prior$t==1))
    prior.pi[iter.c,] <- prior$pi
    prior.pi[,1] <- 0.5

    # update component prop
    gene_index <- rowSums(prior$t)!=0
    gene_weight <-  (gene_weight*(iter.c-1)+gene_index)/iter.c
    w_s <- apply(matrix(1:n.mix),1,function(x)rdirichlet(1,alpha =prior$w[,x]*100+1))

    mu.mix.1 <- prior$alpha%*%matrix(1,1,n.mix)+
                       cbind(prior$beta,prior$beta.miss)%*%w_s
    mu.mix.2 <- prior$alpha%*%matrix(1,1,n.mix)+
                       cbind(prior$beta,prior$beta.miss)%*%prior$w

    density1 <- (dnorm(mix,mean=mu.mix.1,sd=sqrt(prior$sigma2.mix*N.total)%*%matrix(1,1,n.mix),log = T))
    density2 <- (dnorm(mix,mean=mu.mix.2,sd=sqrt(prior$sigma2.mix*N.total)%*%matrix(1,1,n.mix),log = T))
    density1[density1==-Inf] <- 0; density2[density2==-Inf] <- 0
    r1 <- colSums(density1*gene_weight)
    r2 <- colSums(density2*gene_weight)
    accept_res <- exp(r1-r2)>runif(n.mix)

    prior$w[,accept_res] <- w_s[,accept_res]
    prior.w[,,iter.c] <- prior$w

    gc()

    iter.c <- iter.c+1
    if(iter.c > iter){stop.sign <- T}
  }

  #output
    structure(list(w=prior.w,t=prior.t,pi=prior.pi,mean=prior.mean/sqrt(N.total/N.sigma2)))
}


###
BayICE_s <- function(ref.set,mix.set,ref.id,iter=10^4,parC=10^-5,burnin=0.6,thinning=3){
  require(MCMCpack)
  require(invgamma)
  require(quadprog)
  require(MASS)
  require(tibble)
  require(mvtnorm)
  #  require(bigmemory)
  #  options(bigmemory.typecast.warning=FALSE)

  # a sequence of indexes
  seqMCMC0 <- c(round(iter*burnin):iter)
  seqMCMC <- seqMCMC0[seq(1,length(seqMCMC0),by=thinning)]

  if(nrow(ref.set)!=nrow(mix.set)){
    stop("the gene lists of reference set and mixed set are inconsistent")
  }else{n.gene <- nrow(ref.set)}
  id.numb <- 1
  id.type <- rep(NA,length(ref.id))
  for(id in unique(ref.id)){
    id.type[ref.id==id] <- id.numb
    id.numb <- id.numb+1
  }

  n.type <- length(unique(id.type));
  n.mix <- ncol(mix.set);
  n.ref <- ncol(ref.set);
  N.total <- n.mix+n.ref;
  type.design <- matrix(0,n.type,n.ref); type.design[cbind(id.type,1:n.ref)] <- 1
  uniq.id.type <- as.numeric(unique(id.type))

  ff.mean.ini0 <- function(x){
    if(sum(id.type==x)==1){
      return(ref.set[,id.type==x])
    }else{
      return(rowMeans(ref.set[,id.type==x]))
    }
  }
  mean.ini0 <- apply(as.matrix(uniq.id.type),1,ff.mean.ini0)
  N.sigma2 <- rowSums((ref.set-mean.ini0%*%type.design)^2)/(n.ref-n.type)
  ref <- ref.set*sqrt(N.total/N.sigma2)
  mix <- mix.set*sqrt(N.total/N.sigma2)

  #note: set data storage
  prior.t <- prior.mean <- 0
  prior.pi <- 0
  prior.w <- 0

  #note: set initial value
  hyper.par <- list(c=parC,a0=5,b0=50,c0=1,d0=1)

  ff.mean.ini <- function(x){
    if(sum(id.type==x)==1){
      return(ref[,id.type==x])
    }else{
      return(rowMeans(ref[,id.type==x]))
    }
  }
  mean.ini <- apply(as.matrix(uniq.id.type),1,ff.mean.ini)

  #note: identify the initial signature set
  p.v<- matrix(NA,nrow(ref),n.type)
  for(i in 1:n.type){
    cell.v <- rep("A",n.ref); cell.v[id.type==uniq.id.type[i]] <- "B"
    if(sum(id.type==uniq.id.type[i])==1){
      p.v[,i] <- apply(ref,1,function(x)t.test(x[cell.v=="A"],mu=x[cell.v=="B"])$p.value)
    }else{
      p.v[,i] <- apply(ref,1,function(x)t.test(x[cell.v=="A"],x[cell.v=="B"])$p.value)
    }
  }

  t.ini <- matrix(0,n.gene,n.type+1)
  for(i in 1:n.type){
    t.ini[order(p.v[,i])[1:200],i] <- 1
  }

  #note: set the initial proportions
  nL.w <- function(x,y,z){
    w.m <- x
    value.nL <- sum((y-z%*%matrix(w.m))^2)
    value.nL
  }
  optim.w <- function(y,z,eta){
    ff <- constrOptim(rep(1/(n.type+1),n.type), f=nL.w,grad=NULL,y=y,z=z,ui = rbind(-rep(1,n.type), diag(n.type)), ci = c(-1,rep(0,n.type)))
    opt.w <- ff$par
    c(opt.w,1-sum(opt.w))
  }
  res.lm.const <- apply(mix[rowSums(t.ini)!=0,],2,optim.w,z=mean.ini[rowSums(t.ini)!=0,])
  res.lm.const[n.type+1,] <- 0.2
  res.lm.const <- apply(res.lm.const,2,function(x)x/sum(x))

  beta.miss.ini <- rowMeans( t(t(mix-mean.ini%*%res.lm.const[-c(n.type+1),])/res.lm.const[c(n.type+1),])) - rowMeans(as.matrix(ref[,id.type==1]))
  beta.ini <- (mean.ini)-rowMeans(as.matrix(ref[,id.type==1]))

  mean.miss <- t(t((mix-mean.ini%*%res.lm.const[-c(n.type+1),]))/res.lm.const[c(n.type+1),])
  p.v.miss <- apply(matrix(1:n.gene),1,function(x)t.test(ref[x,],mean.miss[x,])$p.value)
  t.ini[order(p.v.miss)[1:200],n.type+1] <- 1

  # set initial variance
  # mark #
  var.ref.ini <- rowMeans((cbind(ref-mean.ini%*%type.design))^2)
  var.mix.ini <- apply(mix,1,var)


  prior <- list(
    alpha=rowMeans(as.matrix(ref[,id.type==1])),
    beta=beta.ini,
    beta.miss=beta.miss.ini,
    psi2=(cbind(beta.ini,beta.miss.ini)+10^-3)^2,
    t=t.ini,
    w=res.lm.const,
    pi=rep(0.3,n.type+1),
    sigma2.ref=var.ref.ini/N.total,
    sigma2.mix=var.mix.ini/N.total
  )

  gene_weight <- rowSums(t.ini)!=0
  iter.c <- 1; stop.sign <- F; count <- 1
  run.time <- Sys.time()
  while(stop.sign==F){
    if(count==floor(iter.c/iter*100)){
      tt <- difftime(Sys.time(),run.time,units ="mins")
      print(paste(count,"%","--","Expected running time =",round((100-count)*tt,3),"mins"))  ;
      count <- count +1
      run.time <- Sys.time()
    }

    # update alpha
    #alpha_mean <- ( rowSums(ref-cbind(prior$beta)%*%type.design)+rowSums(mix-cbind(prior$beta,prior$beta.miss)%*%prior$w) )/(n.ref+n.mix)
    alpha_mean <- rowMeans(ref[,id.type==1])
    prior$alpha <- rnorm(n.gene,alpha_mean,sqrt(prior$sigma2.ref*N.total/(n.ref+n.mix)))

    # update beta
    beta.sampling <- function(d){
      w1.beta <- (prior$sigma2.ref*N.total/sum(id.type==d))^-1
      w2.beta <- (prior$sigma2.mix*N.total/sum(prior$w[d,]^2))^-1
      w3.beta <- ( ( (prior$t[,d]+(1-prior$t[,d])*hyper.par$c)*(prior$psi2[,d]) )^-1)*((prior$sigma2.mix*N.total)^-1)
      pos.mean.beta <-
        w1.beta/(w1.beta+w2.beta+w3.beta)*
        rowMeans( as.matrix( (ref-prior$alpha%*%matrix(1,1,n.ref))[,id.type==d] ) )+
        w2.beta/(w1.beta+w2.beta+w3.beta)*
        (mix-prior$alpha%*%matrix(1,1,n.mix)-cbind(prior$beta,prior$beta.miss)[,-d]%*%prior$w[-d,])%*%matrix(prior$w[d,])/sum(prior$w[d,]^2)
      pos.var.beta <- 1/(w1.beta+w2.beta+w3.beta)
      rnorm(n.gene,pos.mean.beta,sqrt(pos.var.beta))
    }
    prior$beta <- apply(as.matrix(uniq.id.type),1,beta.sampling)
    prior$beta[,1] <- 0

    # update beta.miss
    beta.miss.sampling <- function(d){
      w2.beta <- (prior$sigma2.mix*N.total/sum(prior$w[n.type+1,]^2))^-1
      w3.beta <- ( ( (prior$t[,(n.type+1)]+(1-prior$t[,(n.type+1)])*hyper.par$c)*
                       (prior$psi2[,(n.type+1)]) )^-1)*((prior$sigma2.mix*N.total)^-1)
      pos.mean.beta <-
        w2.beta/(w2.beta+w3.beta)*
        (mix-prior$alpha%*%matrix(1,1,n.mix)-cbind(prior$beta,prior$beta.miss)[,-(n.type+1)]%*%prior$w[-(n.type+1),])%*%matrix(prior$w[(n.type+1),])/sum(prior$w[(n.type+1),]^2)
      pos.var.beta <- 1/(w2.beta+w3.beta)
      rnorm(n.gene,pos.mean.beta,sqrt(pos.var.beta))
    }
    prior$beta.miss <- beta.miss.sampling(1)

    effect_size <- cbind(prior$beta,prior$beta.miss)
    prior.mean0 <- effect_size+prior$alpha

    # update sigma2.ref
    rss.ref <- rowSums( (ref-prior$alpha%*%matrix(1,1,n.ref)-cbind(prior$beta)%*%type.design)^2)/N.total
    prior$sigma2.ref <- rinvgamma(n.gene,(n.ref)/2-1,(rss.ref)/2)

    # update sigma2
    rss.mix <- rowSums( (mix-prior$alpha%*%matrix(1,1,n.mix)-cbind(prior$beta,prior$beta.miss)%*%prior$w)^2 )/N.total
    rss.beta <- rowSums( effect_size^2/(prior$t+hyper.par$c*(1-prior$t))/prior$psi2 )/N.total
    prior$sigma2.mix <- rinvgamma(n.gene,(n.mix+n.type+1)/2-1,(rss.mix+rss.beta)/2)


    # sampling t
    density.beta1 <- dnorm(effect_size,0,sqrt(prior$psi2*((prior$sigma2.mix*N.total)%*%matrix(1,1,n.type+1)) ),log = T) +
      log(matrix(1,n.gene,1)%*%prior$pi)
    density.beta2 <-  dnorm(effect_size,0,sqrt(prior$psi2*((prior$sigma2.mix*N.total)%*%matrix(1,1,n.type+1))*hyper.par$c),log = T) +
      log(matrix(1,n.gene,1)%*%(1-prior$pi))
    prob.beta <-  1/(1+exp(density.beta2-density.beta1))
    prior$t <-  matrix(0,n.gene,n.type+1)
    t.sample <- matrix(rbinom(n=n.gene*(n.type+1),1,prob=prob.beta),n.gene,(n.type+1))
    prior$t[t.sample==1] <- 1
    prior$t[,1] <- 0
    prior.t0 <- prior$t

    # update psi2
    par.invG.psi2 <- (effect_size)^2/(prior$t+hyper.par$c*(1-prior$t))/prior$sigma2.mix%*%matrix(1,1,n.type+1)/N.total
    prior$psi2 <- matrix(rinvgamma(n.gene*(n.type+1),hyper.par$a0+(1)/2,hyper.par$b0+par.invG.psi2),
                         n.gene,n.type+1)

    # update pi
    prior$pi <- rbeta((n.type+1),hyper.par$c0+colSums(prior$t==1),
                      hyper.par$d0+n.gene-colSums(prior$t==1))
    prior.pi0 <- prior$pi
    prior.pi0[1] <- 0.5

    # update component prop
    gene_index <- rowSums(prior$t)!=0
    gene_weight <-  (gene_weight*(iter.c-1)+gene_index)/iter.c
    w_s <- apply(matrix(1:n.mix),1,function(x)rdirichlet(1,alpha =prior$w[,x]*100+1))

    mu.mix.1 <- prior$alpha%*%matrix(1,1,n.mix)+
      cbind(prior$beta,prior$beta.miss)%*%w_s
    mu.mix.2 <- prior$alpha%*%matrix(1,1,n.mix)+
      cbind(prior$beta,prior$beta.miss)%*%prior$w

    density1 <- (dnorm(mix,mean=mu.mix.1,sd=sqrt(prior$sigma2.mix*N.total)%*%matrix(1,1,n.mix),log = T))
    density2 <- (dnorm(mix,mean=mu.mix.2,sd=sqrt(prior$sigma2.mix*N.total)%*%matrix(1,1,n.mix),log = T))
    density1[density1==-Inf] <- 0; density2[density2==-Inf] <- 0
    r1 <- colSums(density1*gene_weight)
    r2 <- colSums(density2*gene_weight)
    accept_res <- exp(r1-r2)>runif(n.mix)

    prior$w[,accept_res] <- w_s[,accept_res]
    prior.w0 <- prior$w

    if(iter.c %in% seqMCMC){
      prior.mean <- prior.mean + prior.mean0
      prior.t <- prior.t + prior.t0
      prior.pi <- prior.pi + prior.pi0
      prior.w <- prior.w + prior.w0
    }

    gc()

    iter.c <- iter.c+1
    if(iter.c > iter){stop.sign <- T}
  }

  prior.mean <- prior.mean/length(seqMCMC)
  prior.t <- prior.t/length(seqMCMC)
  prior.pi <- prior.pi/length(seqMCMC)
  prior.w <- prior.w/length(seqMCMC)
  #output
  structure(list(w=prior.w,t=prior.t,pi=prior.pi,mean=prior.mean/sqrt(N.total/N.sigma2)))
}

