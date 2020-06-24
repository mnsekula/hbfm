#####################################
## User Functions
#####################################


#' Stochastic EM algorithm for initial parameter estimation
#' 
#' Function used to implement the stochastic EM algorithm defined in the
#' "A sparse Bayesian factor model for the construction of gene co-expression 
#' networks from single-cell RNA sequencing count data" manuscript.
#' 
#' This algorithm should be used before \code{hbfm.fit} as it generates initial parameter values
#' for use in \code{hbfm.fit}.
#'
#' @param Y data.frame or matrix of gene expression counts where the rows correspond to genes and
#'          columns correspond to cells; Y must contain integers and have row names
#' @param Fac number of factors to consider in the model; only a single number is accepted
#' @param M.stoc total number of stochastic EM iterations
#' @param M.int initial number of MCMC draws before maximization
#' @param M.eval number of iterations to be used for parameter estimation; the final \code{M.eval} iterations from
#'               the \code{M.stoc} number of stochastic EM iterations will be considered
#' @param M.ll.seq intervals for calculating marginal log-likelihood before final \code{M.eval} iterations; used 
#'                 to reduce computational time
#' @param H number of lambda draws for marginal log-likelihood calculation
#' @param seed seed for random number generation 
#' @param verbose if TRUE, \code{stoc.em} status is displayed at every 100th iteration
#'
#' @return hbfm.par-class object containing:
#' 
#' \itemize{ 
#'   \item{Y}{: data.frame or matrix of gene expression counts}
#'   \item{Fac}{: number of factors considered in the model}
#'   \item{beta}{: initial beta parameter vector for input into \code{hbfm.fit}}
#'   \item{theta}{: initial theta parameter vector for input into \code{hbfm.fit}}
#'   \item{alpha}{: initial alpha parameter matrix for input into \code{hbfm.fit}}
#'   \item{lambda}{: initial lambda parameter matrix for input into \code{hbfm.fit}}
#'   \item{phi}{: initial phi parameter vector for input into \code{hbfm.fit}}
#'   \item{h1}{: location hyperparameter of lognormal prior for phi to input into \code{hbfm.fit}}
#'   \item{h2}{: scale hyperparameter of lognormal prior for phi to input into \code{hbfm.fit}}
#'   \item{mll}{: calculated marginal log-likelihoods at given iterations; used for convergence diagnostics}
#' }
#' 
#' @examples
#' \dontrun{
#' ## Load dataset
#' data(gene.mat)
#' 
#' ## Run stochastic EM first
#' ## Consider F=5 factors
#' fit1 <- stoc.em(Y=gene.mat, Fac = 5)
#' 
#' }
#' 
#' @export
#'


stoc.em <- function(Y, Fac, M.stoc = 2000, M.int = 100, M.eval = 200, M.ll.seq = 10, H=50, seed = 123, verbose = FALSE){# Start of function
  ###### Hyperparmeters from manuscript
  ###### Add as arguments (later)
  GIG.adj = 0.9
  phi.scale = 0.25 #lognormal scale parameter used when drawing candidate values of phi (only used in beginning M.int iterations)
  a <- 0.001 #Shape hyperparameter for beta
  b <- 1000 #Scale hyperparameter for beta
  a_star <- 1#Shape1 hyperparameter for theta
  b_star <- 1#Shape2 hyperparameter for theta
  h1_var <- 100 #Variance hyperparameter for h1
  h2_ab <- 1 #Hyperparameters for h2
  ######################
  
  Y <- data.matrix(Y)
  G <- nrow(Y)
  N <- ncol(Y)
  max_iter <- M.int + 1:M.stoc #iterations for maximization steps
  M <- M.int + length(max_iter) #total number of iterations
  
  Yg_sum <- rowSums(Y)
  # Determine which iterations to calculate marginal log likelihood
  log_lik_save <- c(seq(1,M-M.eval,by=M.ll.seq),(M-M.eval):M) 
  
  ## Intialize starting values
  set.seed(seed)
  beta <- rgamma(G,2,scale=1/2)
  theta <- rep(0.4,Fac)
  alpha <- matrix(NA,nrow = G,ncol = Fac)
  h1 <- runif(1, min = -0.5, max = -0.2)
  h2 <- runif(1, min = 0.84, max = 0.9)
  phi <- rlnorm(Fac,h1,sqrt(h2))
  for(f in 1:Fac){
    alpha[,f] <- rbinom(G,1,theta[f])*sample(c(-1,1),G,replace=TRUE)
  } 
  l_Y <- log(Y+1) - log(beta) #look into -log(beta)
  mod1 <- lm(l_Y ~ alpha)
  lambda <- exp(t(coef(mod1)[2:(Fac+1),]))
  colnames(lambda) <- NULL
  
  
  ## Generate parameter matrices and arrays
  beta_samp <- matrix(NA,nrow = M, ncol = G)
  phi_samp <- matrix(NA, nrow = M, ncol = Fac)
  theta_samp <- matrix(NA,nrow = M, ncol = Fac)
  lambda_samp <- array(NA,dim=c(N,Fac,M))
  alpha_samp <- array(NA,dim=c(G,Fac,M))
  log_lik_int_samp <- array(NA, dim = c(G,N,M)) #DOES THIS NEED TO BE AN ARRAY?
  #log_lik_cell_samp <- matrix(NA, nrow = length((M-999):M), ncol = N) #Save last 1000 iterations for WAIC calculations (REMOVE?)
  log_lik_tot_samp <- rep(NA,M)
  
  
  #####################################
  ## Stochastic EM code
  #####################################
  for(m in 1:M){
    ## beta sampling
    beta_shape <- Yg_sum + a
    beta_scale <- (colSums(lPOWa(a= alpha, tl = t(lambda), p = phi)) + 1/b)^-1
    if(m %in% max_iter){
      # Posterior mode of betas
      beta <- (beta_shape-1)*beta_scale
    }else{
      # MCMC sampling
      beta <- rgamma(G,beta_shape,scale=beta_scale)
    }
    
    ## theta sampling
    sum_abs_alpha <- colSums(abs(alpha))
    theta_a_shape <- a_star + sum_abs_alpha
    theta_b_shape <- b_star + G - sum_abs_alpha
    theta <- rbeta(Fac, theta_a_shape, theta_b_shape)
    
    ## alpha sampling
    lambda_Y_prod <- Y%*%log(lambda)
    if(m %in% max_iter){
      # Maximize alpha's
      for(f in sample(Fac)){
        psi_f <- lPOWa(a=alpha[,-f], tl=t(lambda[,-f]), p = phi[-f])
        ## A = 0, B = 1, C = -1
        logA <- log(1-theta[f]) - beta*colSums(psi_f)
        logB <- log(theta[f]) - log(2) + lambda_Y_prod[,f] - phi[f]/2*rowSums(Y) - beta*colSums((exp(-phi[f]/2)*lambda[,f])*psi_f)
        logC <- log(theta[f]) - log(2) - lambda_Y_prod[,f] - phi[f]/2*rowSums(Y) - beta*colSums((exp(-phi[f]/2)/lambda[,f])*psi_f)
        log_ABC <- rbind(logA, logB, logC)
        max_ABC <- apply(log_ABC, 2, function(x){ 
          y <- which(x==max(x))
          return(ifelse(y==1,0,ifelse(y==2,1,-1)))
        })
        alpha[,f] <- max_ABC
      }# End alpha sampling
    }else{
      # MCMC sampling
      for(f in sample(Fac)){
        psi_f <- lPOWa(a=alpha[,-f], tl=t(lambda[,-f]), p = phi[-f])
        log_thetas <- log(theta[f]) - log(2) - log(1-theta[f])
        ## A = 0, B = 1, C = -1
        logB_m_logA <- log_thetas + lambda_Y_prod[,f] - phi[f]/2*rowSums(Y) - beta*colSums((exp(-phi[f]/2)*lambda[,f]-1)*psi_f)
        logC_m_logA <- log_thetas - lambda_Y_prod[,f] - phi[f]/2*rowSums(Y) - beta*colSums((exp(-phi[f]/2)/lambda[,f]-1)*psi_f)
        logB_m_logC <- logB_m_logA - logC_m_logA 
        expB_p_expC <- exp(logB_m_logA) + exp(logC_m_logA)
        expB_m_C <- exp(logB_m_logC)
        alpha_P_abs1 <- expB_p_expC/(expB_p_expC+1)
        alpha_P1 <- expB_m_C/(expB_m_C+1)
        ## Adjust for NAN values
        alpha_P_abs1[which(is.nan(alpha_P_abs1))] <- 1
        alpha_P1[which(is.nan(alpha_P1))] <- 1
        alpha[,f] <- rbinom(G,1,alpha_P_abs1)
        abs1_ind <- which(alpha[,f] != 0)
        if(length(abs1_ind) > 0){
          alpha[abs1_ind,f] <- 2*rbinom(length(abs1_ind),1,alpha_P1[abs1_ind])-1
        }
      }# End alpha sampling
    }
    
    ## lambda sampling
    for(f in sample(Fac)){
      num_ind <- sum(alpha[,f] != 0) #number of non-zero alphas
      if(num_ind > 0){
        lambda_lam <- colSums(alpha[,f]*Y) + 0.01 #0.01 is adjustment
        a1_beta <- (alpha[,f] == 1)*beta*exp(-phi[f]/2)
        an1_beta <- (alpha[,f] == -1)*beta*exp(-phi[f]/2)
      }
      for(i in 1:N){
        if(num_ind == 0){
          ## no non-zero alphas
          ## draw new lambda values from LN
          cand_lamb <- rlnorm(1,0,sqrt(phi[f]))
          MH_logprob_lamb <- 0
        }else{
          other_fs <- lVECa(a=t(alpha[,-f]), tl=lambda[i,-f], p = phi[-f]) #cpp function
          lambda_psi <- 2*sum(a1_beta*other_fs) + 0.01 #0.01 is adjustment
          lambda_chi <- 2*sum(an1_beta*other_fs)
          cand_lamb <- GIGrvg::rgig(1,lambda=lambda_lam[i], chi=GIG.adj*lambda_chi, psi=GIG.adj*lambda_psi)
          MH_logprob_lamb <- -1/2*(1-GIG.adj)*(lambda_psi*(cand_lamb - lambda[i,f])+lambda_chi*(1/cand_lamb - 1/lambda[i,f])) -
            1/(2*phi[f])*(log(cand_lamb)^2-log(lambda[i,f])^2)
        }
        
        if(is.na(MH_logprob_lamb)){
          MH_logprob_lamb <- -Inf
        }
        if(runif(1) < exp(MH_logprob_lamb)){
          lambda[i,f] <- cand_lamb
        }
      }# End cell loop
    }# End lambda sampling
    
    # Hyperparameters for phi sampling
    h1_prec <- Fac/h2 + 1/h1_var
    h1 <- rnorm(1,mean(log(phi))*(Fac/h2)/h1_prec,(h1_prec)^(-1/2)) 
    h2 <- MCMCpack::rinvgamma(1,(Fac/2+h2_ab-1)+1, (h2_ab + sum((log(phi)-h1)^2)/2)) 
    
    # phi sampling
    num_ind <- colSums(alpha != 0)
    if(m %in% max_iter){
      for(f in sample(Fac)){
        # optimize function
        psi_f <- lPOWa(a=alpha[,-f], tl=t(lambda[,-f]), p = phi[-f])
        phi[f] <- optimize(phi_log_lik, c(0.01, 10), maximum = TRUE, tol=0.01, psi_f = psi_f, fac=f, alpha=alpha, beta=beta, 
                           lambda =lambda, h1=h1, h2=h2, Y=Y, N=N, G=G)$maximum
      }
    }else{
      # MH sampling
      for(f in sample(Fac)){
        if(num_ind[f] > 0){ #Only use cand phi if factor is active
          cand <- rlnorm(1,log(phi[f]), phi.scale)
          psi_f <- lPOWa(a=alpha[,-f], tl=t(lambda[,-f]), p = phi[-f])
          loglik_cand <- phi_log_lik(cand = cand, psi_f = psi_f, fac=f, alpha=alpha, beta=beta, lambda =lambda, h1=h1, h2=h2, Y=Y, N=N, G=G)
          loglik_phi <- phi_log_lik(cand = phi[f], psi_f = psi_f, fac=f, alpha=alpha, beta=beta, lambda =lambda, h1=h1, h2=h2, Y=Y, N=N, G=G)
          MH_logprob <- loglik_cand + dlnorm(phi[f],log(cand),phi.scale,log=TRUE) - loglik_phi - dlnorm(cand,log(phi[f]),phi.scale,log=TRUE)
          if(runif(1) < exp(MH_logprob)){
            phi[f] <- cand
          }
        }else{
          phi_null <- rlnorm(1,h1,sqrt(h2)) #Pull new value from prior
          phi[f] <- phi_null
        }
      }
    }
    
    beta_samp[m,] <- beta
    theta_samp[m,] <- theta
    alpha_samp[,,m] <- alpha
    lambda_samp[,,m] <- lambda
    phi_samp[m,] <- phi
    
    if(m %in% log_lik_save){
      log_lik_int_samp[,,m] <- margLik(nh=H,a=alpha,b=beta,p=phi,Y=Y) #DOES THIS NEED TO BE AN ARRAY?
      log_lik_tot_samp[m] <- sum(log_lik_int_samp[,,m])
    }
    
    ## Get optimized parameter estimates
    if(m == M){
      beta <- apply(beta_samp[tail(max_iter,M.eval),],2,mean)
      theta <- apply(theta_samp[tail(max_iter,M.eval),],2,mean)
      alpha <- apply(alpha_samp[,,tail(max_iter,M.eval)],1:2,num_mode)
      lambda <- apply(lambda_samp[,,tail(max_iter,M.eval)],1:2,mean)
      phi <- apply(phi_samp[tail(max_iter,M.eval),],2,mean)
    }
    
    ## Remove later
    if(verbose & m%%100 == 0){
      print(paste("Sample", m, "completed"))
    }
    
  }# End Stocastic EM
  
  res <- list(Y=Y, Fac=Fac, beta=beta, theta=theta, alpha=alpha, lambda=lambda, phi=phi, h1=h1, h2=h2, mll=log_lik_tot_samp)
  class(res) <- "hbfm.par"
  ## also return log_lik_tot_samp
  return(res)
  
  
} #End of function





#' Hierarchical Bayesian factor model MCMC sampler
#' 
#' Function used to implement the MCMC sampler for the hierarchical Bayesian factor model (HBFM) 
#' defined in the "A sparse Bayesian factor model for the construction of gene co-expression 
#' networks from single-cell RNA sequencing count data" manuscript.
#' 
#' @param stoc.em.param object of hbfm.par-class created by the \code{stoc.em} function
#' @param M total number of MCMC iterations
#' @param M.save number of iterations to be used for correlation estimation; the final \code{M.save} iterations from
#'               the \code{M} number of MCMC iterations will be saved for further analysis
#' @param M.ll.seq intervals for calculating marginal log-likelihood before final \code{M.save} iterations; used 
#'                 to reduce computational time
#' @param H number of lambda draws for marginal log-likelihood calculation
#' @param phi.scale lognormal scale parameter used when drawing candidate values of phi in
#'                  Metropolis-Hastings step; can be used to adjust acceptance rate of phi samples
#' @param par.samp if TRUE, samples of each parameter in \code{hbfm.fit} for each MCMC iteration are saved; if FALSE
#'                 only the calculated parameters of correlation and marginal log-likelihood are saved
#' @param seed seed for random number generation 
#' @param verbose if TRUE, \code{hbfm.fit} status is displayed at every 100th iteration
#'
#' @return hbfm.fit-class object containing:
#' 
#'  \itemize{ 
#'   \item{Y}{: data.frame or matrix of gene expression counts}
#'   \item{Fac}{: number of factors considered in the model}
#'   \item{samples}{: samples of estimated parameters and/or calculated parameters from MCMC iterations}
#'   \itemize{
#'   \item{mll}{: calculated marginal log-likelihoods}
#'   \item{corr}{: samples of gene-gene correlation matrix}
#'   \item{beta}{: samples of beta parameters; included only when \code{par.samp = TRUE}}
#'   \item{theta}{: samples of theta parameters; included only when \code{par.samp = TRUE}}
#'   \item{alpha}{: samples of alpha parameters; included only when \code{par.samp = TRUE}}
#'   \item{lambda}{: samples of lambda parameters; included only when \code{par.samp = TRUE}}
#'   \item{phi}{: samples of phi parameters; included only when \code{par.samp = TRUE}}
#'   \item{h1}{: samples of location hyperparameter of lognormal prior for phi; included only when \code{par.samp = TRUE}}
#'   \item{h2}{: samples of scale hyperparameter of lognormal prior for phi; included only when \code{par.samp = TRUE}}
#'   }
#' }
#' 
#' @examples
#' \dontrun{
#' ## Load dataset
#' data(gene.mat)
#' 
#' ## Run stochastic EM first
#' ## Consider F=5 factors
#' fit1 <- stoc.em(Y=gene.mat, Fac = 5)
#' 
#' ## Run MCMC sampler with initial parameter values from stoc.em
#' fit.res1 <- hbfm.fit(fit1)
#' }
#' 
#' 
#' @export
#'

hbfm.fit <- function(stoc.em.param, M = 4000, M.save = 1000, M.ll.seq = 10, H=50, phi.scale = 0.25, par.samp = FALSE, seed = 123, verbose = FALSE){# Start of function
  ###### Hyperparmeters from manuscript
  ###### Add as arguments (later)
  GIG.adj = 0.9
  a <- 0.001 #Shape hyperparameter for beta
  b <- 1000 #Scale hyperparameter for beta
  a_star <- 1#Shape1 hyperparameter for theta
  b_star <- 1#Shape2 hyperparameter for theta
  h1_var <- 100 #Variance hyperparameter for h1
  h2_ab <- 1 #Hyperparameters for h2
  ######################
  
  Y <- data.matrix(stoc.em.param$Y)
  G <- nrow(Y)
  N <- ncol(Y)
  Fac <- stoc.em.param$Fac
  
  Yg_sum <- rowSums(Y)
  gene_names <- rownames(Y)
  # Determine which iterations to calculate marginal log likelihood
  log_lik_save <- c(seq(1,M-M.save,by=M.ll.seq),(M-M.save):M) 
  
  ## Obtain starting values from stoc.em.param
  beta <- stoc.em.param$beta
  theta <- stoc.em.param$theta
  alpha <- stoc.em.param$alpha
  h1 <- stoc.em.param$h1
  h2 <- stoc.em.param$h2
  phi <- stoc.em.param$phi
  lambda <- stoc.em.param$lambda
  
  ## Generate parameter matrices and arrays
  beta_samp <- matrix(NA,nrow = M, ncol = G)
  phi_samp <- matrix(NA, nrow = M, ncol = Fac)
  theta_samp <- matrix(NA,nrow = M, ncol = Fac)
  lambda_samp <- array(NA,dim=c(N,Fac,M))
  alpha_samp <- array(NA,dim=c(G,Fac,M))
  corr_samp <- array(NA, dim = c(G,G,M.save), dimnames = list(rownames(Y),rownames(Y)))
  log_lik_int_samp <- array(NA, dim = c(G,N,M))
  #log_lik_cell_samp <- matrix(NA, nrow = length((M-M.save+1):M), ncol = N) #Save for WAIC calculations (REMOVE?)
  log_lik_tot_samp <- rep(NA,M)
  
  set.seed(seed)
  #####################################
  ## MCMC sampling code
  #####################################
  for(m in 1:M){
    ## beta sampling
    beta_shape <- Yg_sum + a
    beta_scale <- (colSums(lPOWa(a= alpha, tl = t(lambda), p = phi)) + 1/b)^-1
    beta <- rgamma(G,beta_shape,scale=beta_scale)
    
    ## theta sampling
    sum_abs_alpha <- colSums(abs(alpha))
    theta_a_shape <- a_star + sum_abs_alpha
    theta_b_shape <- b_star + G - sum_abs_alpha
    theta <- rbeta(Fac, theta_a_shape, theta_b_shape)
    
    ## alpha sampling
    lambda_Y_prod <- Y%*%log(lambda)
    for(f in sample(Fac)){
      psi_f <- lPOWa(a=alpha[,-f], tl=t(lambda[,-f]), p = phi[-f])
      log_thetas <- log(theta[f]) - log(2) - log(1-theta[f])
      ## A = 0, B = 1, C = -1
      logB_m_logA <- log_thetas + lambda_Y_prod[,f] - phi[f]/2*rowSums(Y) - beta*colSums((exp(-phi[f]/2)*lambda[,f]-1)*psi_f)
      logC_m_logA <- log_thetas - lambda_Y_prod[,f] - phi[f]/2*rowSums(Y) - beta*colSums((exp(-phi[f]/2)/lambda[,f]-1)*psi_f)
      logB_m_logC <- logB_m_logA - logC_m_logA 
      expB_p_expC <- exp(logB_m_logA) + exp(logC_m_logA)
      expB_m_C <- exp(logB_m_logC)
      alpha_P_abs1 <- expB_p_expC/(expB_p_expC+1)
      alpha_P1 <- expB_m_C/(expB_m_C+1)
      ## Adjust for NAN values
      alpha_P_abs1[which(is.nan(alpha_P_abs1))] <- 1
      alpha_P1[which(is.nan(alpha_P1))] <- 1
      alpha[,f] <- rbinom(G,1,alpha_P_abs1)
      abs1_ind <- which(alpha[,f] != 0)
      if(length(abs1_ind) > 0){
        alpha[abs1_ind,f] <- 2*rbinom(length(abs1_ind),1,alpha_P1[abs1_ind])-1
      }
    }# End alpha sampling
    
    ## lambda sampling
    for(f in sample(Fac)){
      num_ind <- sum(alpha[,f] != 0) #number of non-zero alphas
      if(num_ind > 0){
        lambda_lam <- colSums(alpha[,f]*Y) + 0.01 #0.01 is adjustment
        a1_beta <- (alpha[,f] == 1)*beta*exp(-phi[f]/2)
        an1_beta <- (alpha[,f] == -1)*beta*exp(-phi[f]/2)
      }
      for(i in 1:N){
        if(num_ind == 0){
          ## no non-zero alphas
          ## draw new lambda values from LN
          cand_lamb <- rlnorm(1,0,sqrt(phi[f]))
          MH_logprob_lamb <- 0
        }else{
          other_fs <- lVECa(a=t(alpha[,-f]), tl=lambda[i,-f], p = phi[-f]) #cpp function
          lambda_psi <- 2*sum(a1_beta*other_fs) + 0.01 #0.01 is adjustment
          lambda_chi <- 2*sum(an1_beta*other_fs)
          cand_lamb <- GIGrvg::rgig(1,lambda=lambda_lam[i], chi=GIG.adj*lambda_chi, psi=GIG.adj*lambda_psi)
          MH_logprob_lamb <- -1/2*(1-GIG.adj)*(lambda_psi*(cand_lamb - lambda[i,f])+lambda_chi*(1/cand_lamb - 1/lambda[i,f])) -
            1/(2*phi[f])*(log(cand_lamb)^2-log(lambda[i,f])^2)
        }
        
        if(is.na(MH_logprob_lamb)){
          MH_logprob_lamb <- -Inf
        }
        if(runif(1) < exp(MH_logprob_lamb)){
          lambda[i,f] <- cand_lamb
        }
      }# End cell loop
    }# End lambda sampling
    
    # Hyperparameters for phi sampling
    h1_prec <- Fac/h2 + 1/h1_var
    h1 <- rnorm(1,mean(log(phi))*(Fac/h2)/h1_prec,(h1_prec)^(-1/2)) #non-adjusted
    h2 <- MCMCpack::rinvgamma(1,(Fac/2+h2_ab-1)+1, (h2_ab + sum((log(phi)-h1)^2)/2)) #non-adjusted
    
    # phi sampling
    num_ind <- colSums(alpha != 0)
    # MH sampling
    for(f in sample(Fac)){
      if(num_ind[f] > 0){ #Only use cand phi if factor is active
        cand <- rlnorm(1,log(phi[f]), phi.scale)
        psi_f <- lPOWa(a=alpha[,-f], tl=t(lambda[,-f]), p = phi[-f])
        loglik_cand <- phi_log_lik(cand = cand, psi_f = psi_f, fac=f, alpha=alpha, beta=beta, lambda =lambda, h1=h1, h2=h2, Y=Y, N=N, G=G)
        loglik_phi <- phi_log_lik(cand = phi[f], psi_f = psi_f, fac=f, alpha=alpha, beta=beta, lambda =lambda, h1=h1, h2=h2, Y=Y, N=N, G=G)
        MH_logprob <- loglik_cand + dlnorm(phi[f],log(cand),phi.scale,log=TRUE) - loglik_phi - dlnorm(cand,log(phi[f]),phi.scale,log=TRUE)
        if(runif(1) < exp(MH_logprob)){
          phi[f] <- cand
        }
      }else{
        phi_null <- rlnorm(1,h1,sqrt(h2)) #Pull new value from prior
        phi[f] <- phi_null
      }
    }
    
    beta_samp[m,] <- beta
    theta_samp[m,] <- theta
    alpha_samp[,,m] <- alpha
    lambda_samp[,,m] <- lambda
    phi_samp[m,] <- phi
    
    ## Calculate correlation in M.save iterations
    if(m %in% (M-M.save+1):M){
      corr_lmug_lmug <-apply(alpha,1,cov_lmu,phi=phi,t_a=t(alpha))
      var_lmu_ig <- apply(alpha,1,var_lmu,phi=phi)
      corr_samp[,,(m-(M-M.save))] <- corr_lmug_lmug/sqrt(var_lmu_ig * matrix(rep(var_lmu_ig,G),nrow=G,byrow=TRUE))
    }
    
    if(m %in% log_lik_save){
      log_lik_int_samp[,,m] <- margLik(nh=H,a=alpha,b=beta,p=phi,Y=Y)
      log_lik_tot_samp[m] <- sum(log_lik_int_samp[,,m])
      
      #INCLUDE LATER IF WAIC IS OPTION FOR MODEL SELECTION
      # if(m %in% (M-M.save+1):M){
      #   log_lik_cell_samp[m-(M-1000),] <- colSums(log_lik_int_samp[,,m])
      # }
    }
    
    ## Remove later
    if(verbose & m%%100 == 0){
      print(paste("Sample", m, "completed"))
    }
    
  }# End Stocastic EM
  
  calc_samples <- list(mll=log_lik_tot_samp,corr=corr_samp) #samples for calculated parameters (marginal log likelihood, correlation)
  
  if(par.samp){
    par_samples <- list(beta_samp, theta_samp, alpha_samp, lambda_samp, phi_samp,mll=log_lik_tot_samp,corr=corr_samp) #samples for model parameters
    all_samples <- c(par_samples,calc_samples)
  }else{
    all_samples <- calc_samples #only include samples for calculated parameters
  }
  
  res <- list(Y=Y,Fac=Fac,samples=all_samples)
  class(res) <- "hbfm.fit"
  return(res)
  
  
} #End of function





#' Estimate the gene-gene correlation matrix from MCMC samples
#' 
#' Function used to combine and analyze MCMC results from one or more sets of
#' correlation samples generated by the \code{hbfm.fit} function.
#' 
#' @param hbfm.list list where each element contains an hbfm.fit-class object; each element of the list contains an object from a different MCMC chain
#'
#' @return hbfm.corr-class object containing:
#' 
#'  \itemize{ 
#'   \item{corr}{: symmetric estimated correlation matrix generated from MCMC samples}
#'   \item{CI.low}{: symmetric matrix representing the lower bound of the 95\% credible intervals (CI) generated from the MCMC samples}
#'   \item{CI.upp}{: symmetric matrix representing the upper bound of the 95\% credible intervals (CI) generated from the MCMC samples}
#'   \item{CI.eval}{: symmetric logical matrix generated from 95\% credible intervals (CI) from MCMC samples; TRUE means 95\% CI does not include 0 (i.e., significant correlation)}
#'   \item{p.val}{: symmetric matrix consisting of approximate "p-values" generated from MCMC samples}
#' }
#' 
#' @details 
#' The correlation estimate is determined by calculating the average of the posterior samples for each gene-gene pair.
#' 
#' To determine whether the correlation is significant, a 95\% credible interval (CI) is determined from the posterior samples. If the CI
#' includes 0, the correlation is deemed to be non-significant and the "CI.eval" element is FALSE. If the CI does not include 0, the correlation is
#' considered significant and the "CI.eval" element is TRUE.
#' 
#' An alternative measurement of significance is the approximate "p-value" calculation, which is output in the "p.val" matrix. This "p-value" is determined 
#' by finding the smallest "a" value such that the 100(1-a)\% CI contains 0. The corresponding "a" value represents the proportion of posterior distribution 
#' that is outside the smallest credible interval that contains 0.
#'
#' In most cases a significant 95\% CI from "CI.eval" will correspond with an approximate "p-value" < 0.05 from "p.val".
#' 
#' @examples
#' \dontrun{
#' ## Load dataset
#' data(gene.mat)
#' 
#' ## Run stochastic EM first
#' ## Consider F=5 factors
#' fit1 <- stoc.em(Y=gene.mat, Fac = 5)
#' 
#' ## Run MCMC sampler with initial parameter values from stoc.em
#' fit.res1 <- hbfm.fit(fit1)
#' 
#' ## Obtain estimated gene-gene correlations from MCMC samples
#' fit.corr <- corr.est(list(fit.res1))
#' print(fit.corr)
#' 
#' }
#' 
#' @export

corr.est <- function(hbfm.list){
  M.samps <- numeric()
  M.genes <- numeric()
  for(l in 1:length(hbfm.list)){
    M.samps <- c(M.samps, dim(hbfm.list[[l]]$samples$corr)[3])
    M.genes <- c(M.genes, dim(hbfm.list[[l]]$samples$corr)[2])
  }
  corr_mat_samps <- array(NA, dim = c(M.genes[1],M.genes[1],sum(M.samps)))
  M.start <- 1
  for(l in 1:length(hbfm.list)){
    corr_mat_samps[,,M.start:(M.start - 1 + M.samps[l])] <- hbfm.list[[l]]$samples$corr
    M.start <- M.start + M.samps[l]
  }
  gene_names <- rownames(hbfm.list[[1]]$samples$corr)
  corr_mat_est <- apply(corr_mat_samps, 1:2, mean, na.rm=TRUE)
  corr_mat_sd <- apply(corr_mat_samps, 1:2, sd, na.rm=TRUE)
  CI_low_all <- apply(corr_mat_samps,1:2,quantile, probs = c(0.025), na.rm=TRUE)
  CI_up_all <- apply(corr_mat_samps,1:2,quantile, probs = c(0.975), na.rm=TRUE)
  CI_eval_all <- (CI_low_all * CI_up_all) > 0
  p_val <- apply(corr_mat_samps,1:2, p.approx)
  diag(p_val) <- diag(CI_eval_all) <-  NA
  rownames(corr_mat_est) <- colnames(corr_mat_est) <- gene_names
  rownames(CI_eval_all) <- colnames(CI_eval_all) <- gene_names
  rownames(CI_low_all) <- colnames(CI_low_all) <- gene_names
  rownames(CI_up_all) <- colnames(CI_up_all) <- gene_names
  rownames(p_val) <- colnames(p_val) <- gene_names
  
  res <- list(corr = corr_mat_est, CI.low = CI_low_all, CI.upp = CI_up_all, CI.eval = CI_eval_all, p.val = p_val)
  class(res) <- "hbfm.corr"
  return(res)
}





#' Calculate DIC of hbfm
#'
#' Function to calculate the Deviance Information Criterion (DIC) from \code{hbfm.fit} output.
#'
#' @param hbfm.list list where each element contains an hbfm.fit-class object; each element of the list contains an object from a different MCMC chain
#'
#' @return DIC of hbfm with \code{Fac} number of factors
#' 
#' @details 
#' For each of the final \code{M.save} iterations of the MCMC, \code{hbfm.fit} calculates the log marginal likelihood, log(ML).
#'  
#' The DIC is calculated with the formula: \eqn{DIC = mean(D) + pD}; where \eqn{D = -2*log(ML)} and \eqn{pD = var(D)/2}.
#' 
#' @references 
#' Gelman, A., Carlin, J. B., Stern, H. S., & Rubin, D. B. (2004). \emph{Bayesian Data Analysis: 2nd Edition}. Chapman and Hall/CRC.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' ## Load dataset
#' data(gene.mat)
#' 
#' ## Run stochastic EM first
#' ## Consider F=5 factors
#' fit1 <- stoc.em(Y=gene.mat, Fac = 5)
#' 
#' ## Run MCMC sampler with initial parameter values from stoc.em
#' fit.res1 <- hbfm.fit(fit1)
#' 
#' ## Obtain estimated gene-gene correlations from MCMC samples
#' hbfm.DIC(list(fit.res1))
#' 
#' }


hbfm.DIC <- function(hbfm.list){
  mll.samps.tot <- numeric()
  for(l in 1:length(hbfm.list)){
    M <- length(hbfm.list[[l]]$samples$mll)
    M.save <- dim(hbfm.list[[l]]$samples$corr)[3]
    mll.samps.tot <- c(mll.samps.tot, hbfm.list[[l]]$samples$mll[(M-M.save+1):M])
  }
  D <- -2*mll.samps.tot
  D_bar <- mean(D)
  p_D <- var(D)/2
  DIC <- D_bar + p_D
  return(DIC)
}





#' Rank hbfm chains
#'
#' Function to rank the MCMC chains by average marginal log-likelihood.
#'
#' @param hbfm.list list where each element contains an hbfm.fit-class object; each element of the list contains an object from a different MCMC chain
#'
#' @return data.frame that orders the input hbfm.fit-class objects based on average marginal log-likelihood; the numbers in the chain.num column
#' represent the elements (chains) from hbfm.list
#' 
#' 
#' @details 
#' The \code{rank.chains} function is used to rank the input chains based on average marginal log-likelihood in order to identify low performing chains.
#'  
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' ## Load dataset
#' data(gene.mat)
#' 
#' ## Run stochastic EM first
#' ## Consider F=5 factors
#' fit1 <- stoc.em(Y=gene.mat, Fac = 5)
#' 
#' ## Run MCMC sampler with initial parameter values from stoc.em
#' fit.res1 <- hbfm.fit(fit1)
#' 
#' ## Run a second chain
#' fit2 <- stoc.em(Y=gene.mat, Fac = 5, seed = 234)
#' fit.res2 <- hbfm.fit(fit2, seed = 234)
#' 
#' ## Obtain estimated gene-gene correlations from MCMC samples
#' rank.chains(list(fit.res1, fit.res2))
#' 
#' }
#' 


rank.chains <- function(hbfm.list){
  mll.avg <- numeric()
  for(l in 1:length(hbfm.list)){
    M <- length(hbfm.list[[l]]$samples$mll)
    M.save <- dim(hbfm.list[[l]]$samples$corr)[3]
    mll.avg <- c(mll.avg, mean(hbfm.list[[l]]$samples$mll[(M-M.save+1):M]))
  }
  c.ranks <- order(mll.avg, decreasing = TRUE)
  rankings <- data.frame(chain.num = c.ranks, avg.mll = mll.avg[c.ranks])
  return(rankings)
}





#####################################
## S3 methods
#####################################

#' Prints hbfm.corr object
#' @export
#' @noRd
`print.hbfm.corr` <-
  function(x, ...){
    x.frame <- data.frame(gene.name1 = character(),
                          gene.name2 = character(),
                          corr = double(),
                          CI.low = double(),
                          CI.upp = double(),
                          CI.sig = logical(),
                          p.val = character(),
                          stringsAsFactors = FALSE)
    k<-1
    for(i in 1:(nrow(x$corr)-1)){
      for(j in (i+1):ncol(x$corr)){
        x.frame[k,1] <- rownames(x$corr)[i]
        x.frame[k,2] <- colnames(x$corr)[j]
        x.frame[k,3] <- round(x$corr[i,j],4)
        x.frame[k,4] <- round(x$CI.low[i,j],4)
        x.frame[k,5] <- round(x$CI.upp[i,j],4)
        x.frame[k,6] <- ifelse(x$CI.eval[i,j],"  *  ","")
        x.frame[k,7] <- formatC(x$p.val[i,j],format = "e", digits = 2)
        k<-k+1
      }
    }
    print(x.frame)
  }

#' Summary for hbfm.corr object
#' @export
#' @noRd
`summary.hbfm.corr` <-
  function(object, ...){
    print(object)
  }





#####################################
## Internal calculation functions
#####################################

#' Covariance structure for log(mu_gi)
#'
#' @param a_g gene-specific alpha vector; length Fac
#' @param phi phi vector; length fac
#' @param t_a transpose of alpha matrix; Fac x G
#'
#' @return
#' covariance between gene g and the other genes
#' @export
#' @noRd
#' 
cov_lmu <- function(a_g, phi, t_a){
  a_x_a <- a_g*t_a
  colSums(phi*a_x_a)
}


#' Variance structure for log(mu_gi)
#'
#' @param a_g gene-specific alpha vector; length Fac
#' @param phi phi vector; length fac
#'
#' @return
#' variance for gene g
#' 
#' @export
#' @noRd


var_lmu <- function(a_g, phi){
  sum(phi*a_g^2)
}


#' Approximate p-values from correlation posterior 
#' 
#' Find the smallest "a" such that the 100(1-a)% credible interval
#' contains 0.
#'
#' @param x samples from posterior
#'
#' @return
#' approximate p-value
#' @export
#' @noRd


p.approx <- function(x){
  if(median(x) < 0){
    n <- length(which(x > 0))
  } else if(median(x) > 0){
    n <- length(which(x < 0))
  } else{
    # median is 0
    n <- sum(abs(x) > 0)/2
  }
  
  # Check for 0 on boundaries
  if(n == 0 & min(abs(x)) == 0){
      n <- sum(abs(x) == 0)
  }
  
  if(n == 0){
    p.val1 <- .01/(length(x)+2)
  }else{
    p.val1 <- n/(length(x))
  }
  p.val <- 2*p.val1
  return(p.val)
}





#####################################
## Internal MCMC functions
#####################################

#' Calculate proportional posterior log-likelihood for phi candidate value
#'
#' Used for phi optimization step in stoc.em
#'
#' @param cand value to replace in phi vector
#' @param psi_f results from lPOWa function
#' @param fac factor number to evaluate
#' @param alpha current alpha estimates
#' @param beta current beta estimates
#' @param lambda current lambda estimates
#' @param h1 current h1 (hyperparameter) estimate
#' @param h2 current h2 (hyperparameter) estimate
#' @param Y matrix of gene expression counts
#' @param N number of cells
#' @param G number of genes
#'
#' @return
#' Proportional posterior log-likelihood for phi candidate value
#' @export
#' @noRd

phi_log_lik <- function(cand, psi_f, fac, alpha, beta, lambda, h1, h2, Y, N, G){
  loglik_cand <- -cand/2*sum(abs(alpha[,fac])*Y)-
    sum(beta*exp(-cand/2*abs(alpha[,fac]))*colSums(matrix(lambda[,fac]^rep(alpha[,fac],each = N), nrow = N, ncol = G)*psi_f))+
    (-(N/2+1))*log(cand) +
    -1/cand*sum(log(lambda[,fac])^2)/2 -
    ((log(cand)-h1)^2)/(2*h2)
  return(loglik_cand)
}

#' Determine the mode of a vector
#'
#' @param x vector of elements
#'
#' @return
#' mode
#' @export
#' @noRd
#' 


num_mode <- function(x){
  #In case of tie, give preference to most recent value
  ux <- unique(rev(x))
  ux[which.max(tabulate(match(x, ux)))]
}






#####################################
## Datasets
#####################################

#' Example gene expression matrix
#'
#' An example dataset consisting of gene expression counts for G = 10 genes (rows)
#' across N = 100 cells (columns) 
#'
#' @format A matrix with 10 rows and 100 columns
#' 
#' @examples
#' data(gene.mat)
#' rownames(gene.mat)
#' colnames(gene.mat)
#' 
"gene.mat"

