#' @name BNP.eq.predict
#' 
#' @title Prediction step for Bayesian non-parametric model for test equating
#' 
#' @description The 
#' 
#' @param model A 'BNP.eq' object.
#' @param from Numeric. A vector of indices indicating from which patterns equating should be performed.
#'        The covariates involved are integrated out.
#' @param into Numeric. A vector of indices indicating into which patterns equating should be performed.
#'        The covariates involved are integrated out.
#' @param alpha Numeric. Significance for credible bands.
#'
#' @details Lorem ipsum dolor sit amet, consectetur adipiscing elit. 
#' Vivamus finibus vitae eros quis dictum. Donec lacus risus, facilisis quis 
#' tincidunt et, tincidunt sed mi. Nullam ullamcorper eros est, sed 
#' fringilla metus volutpat eu. Etiam ornare nulla id lorem posuere, eu vehicula 
#' urna vestibulum. Quisque luctus, diam ac mattis faucibus, leo felis tincidunt urna, eu 
#' tempus massa neque nec nibh. Aliquam erat volutpat. Fusce tempor mattis enim quis pretium. 
#' Aliquam volutpat luctus felis, nec fringilla enim tincidunt sed. 
#' Nam nec leo quis erat lobortis vulputate ac at neque
#' 
#' @return  
#'  A 'BNP.eq.predict' object, which is list containing the following items:
#' @return pdf A list of PDF's.
#' @return cdf A list of CDF's.
#' @return equ Numeric. Equated values.
#' @return grid Numeric. Grid used to evaluate pdf's and cdf's.
#' 
#' @references asdasd
#' 
#' @author Daniel Leon \email{dnacuna@uc.cl}, Felipe Barrientos \email{afb26@stat.duke.edu}.
#' 
#' @keywords BNP equating, Bayesian non-parametrics, equating
BNP.eq.predict <- function(model, from=NULL, into=NULL, alpha=0.05){
  if(class(model) != "BNP.eq")
    stop("Fitted object must be BNP.eq.")
  res<-print('The BNP.eq.predict() function is currently not available in SNSequate')
  return(res)
}



#' @name BNP.eq
#' 
#' @title Bayesian non-parametric model for test equating
#' 
#' @description The Bayesian nonparametric (BNP) approach (Ghoshand Ramamoorthi, 2003; Hjort et al., 2010) 
#' starts by focusing on spaces of distribution functions, so that uncertainty is expressed on F itself. 
#' The prior distribution p(F) is defined on the space F of all distribution functions defined on X . If X 
#' is an infinite set then F is infinite-dimensional, and the corresponding prior model 
#' p(F) on F is termed nonparametric. The prior probability model is also referred to
#' as a random probability measure (RPM), and it essentially corresponds to a distribution 
#' on the space of all distributions on the set X . Thus Bayesian nonparametric models 
#' are probability models defined on a function space (Muller and Quintana, 2004).
#' 
#' Gonzalez et al. (2015) proposed a Bayesian non-parametric approach for equating. The main
#' idea consists of introducing covariate dependent BNP models for a collection of 
#' covariate-dependent equating transformations
#' 
#'  \eqn{ \left\{ \boldsymbol{\varphi}_{\boldsymbol{z}_f, \boldsymbol{z}_t} (\cdot): 
#'            \boldsymbol{z}_f, \boldsymbol{z}_t \in \mathcal{L}
#'        \right\} 
#'   }
#' 
#' 
#' 
#' 
#' @param scores_x Vector.  Scores of form X.
#' @param scores_y Vector.  Scores of form Y.
#' @param range_scores Vector of length 2.  Represent the minimum and maximum scores in the test.
#' @param design Character.  Only supports 'EG' design now.
#' @param covariates Data.frame.  A data frame with factors, containing covariates 
#'        for test X and Y, stacked in that order.
#' @param prior List.  Prior information for BNP model. 
#'        For more information see \code{\link{DPpackage}}.  
#' @param mcmc List.  MCMC information for BNP model. 
#'        For more information see \code{\link{DPpackage}}.
#' @param normalize Logical.  Whether normalize or not the 
#'        response variable. This is due to Berstein's polynomials. Default is TRUE.
#'
#' @details Lorem ipsum dolor sit amet, consectetur adipiscing elit. 
#' Vivamus finibus vitae eros quis dictum. Donec lacus risus, facilisis quis 
#' tincidunt et, tincidunt sed mi. Nullam ullamcorper eros est, sed 
#' fringilla metus volutpat eu. Etiam ornare nulla id lorem posuere, eu vehicula 
#' urna vestibulum. Quisque luctus, diam ac mattis faucibus, leo felis tincidunt urna, eu 
#' tempus massa neque nec nibh. Aliquam erat volutpat. Fusce tempor mattis enim quis pretium. 
#' Aliquam volutpat luctus felis, nec fringilla enim tincidunt sed. 
#' Nam nec leo quis erat lobortis vulputate ac at neque
#' 
#' @return  
#'   A 'BNP.eq' object, which is list containing the following items:
#' @return Y Response variable.
#' @return X Design Matrix.
#' @return fit DPpackage object. Fitted model with raw samples.
#' @return max_score Maximum score of test.
#' @return patterns A matrix describing the different patterns formed
#'         from the factors in the covariables.
#' @return patterns_freq The normalized frequency of each pattern.
#' 
#' @references asdasd
#' 
#' @author Daniel Leon \email{dnacuna@uc.cl}, Felipe Barrientos \email{afb26@stat.duke.edu}.
#' 
#' @keywords BNP equating, Bayesian non-parametrics, equating
BNP.eq <- function(scores_x, scores_y, range_scores=NULL, design="EG", covariates=NULL, 
                   prior=NULL, mcmc=NULL, normalize=TRUE){

  ## Response
  Y <- c(scores_x, scores_y)

  if(is.null(range_scores)){
    max_score=max(Y)
  } else{
    max_score=range_scores[2]
  }
  ## Bernstein poly assume response on [0, 1]
  if(normalize)
    Y <- Y / max_score

  ## Covariates.
  form_cov <- factor( c( rep("X", length(scores_x)),
                         rep("Y", length(scores_y))) )
  X <- model.matrix(~ ., data=cbind(Form=form_cov, covariates))

  ## Count patterns in the design matrix. Basically it
  ## extract each unique pattern from the discrete covariates
  Grid <- count(X)
  patterns <- Grid[, 1:(ncol(Grid)-1)]
  patterns_freq <- as.numeric(Grid$freq)
  patterns_freq <- patterns_freq/sum(patterns_freq)

#  # Fitting single-atoms LDBPP model
#  # mcmc parameters
#  if(is.null(mcmc))
#    mcmc <-  list(nburn = 5000, nskip = 20, ndisplay = 10, nsave = 1000)

  ## Default params, not needed for our purposes
  # Predictions will be made later
#  grid <- 0.5
#  npred <- 1
#  xpred <-  patterns[1, ]  ## Just to get results. Prediction isn't made by DPpackge

  # List of prior information for DPpackage
#  if(is.null(prior))
#    prior <- list(maxn = 25, a1=1, a2=1,
#                  lambda = 25, nu = 6,
#                  psiinv = diag(1000, ncol(patterns)),
#                  m0 = rep(0, ncol(patterns)),
#                  S0 = diag(1000, ncol(patterns)))

  # State
#  state <- NULL

  # fitting the model
#  cat("Calling 'DPpackage' to sample model: \n\n")
#  fit <- tLDBPPdensity(formula=Y~X[, -1],xpred=xpred,  # We take out the intercept because DPpackage adds it
#                       prior=prior,
#                       mcmc=mcmc,
#                       state=state,status=TRUE,
#                       grid=grid,
#                       compute.band=FALSE,type.band="PD")

#  ret <- list(scores_x=scores_x, scores_y=scores_y,
#              Y=Y, X=X,fit=fit, max_score=max_score,
#              patterns=patterns, patterns_freq=patterns_freq)
#  class(ret) <- "BNP.eq"

  ret<-print('The BNP.eq() is currently not available in SNSequate')
    ret
    }
