#' Perform restricted MLE  (RMLE) to estimate parameters under the constraint defined by the boundary of null hypothesis
#'
#' Perform restricted MLE  (RMLE) to estimate parameters under the constraint defined by the boundary of null hypothesis, \eqn{\mu_T - \mu_R = \eta\sigma_R} where \eqn{\eta} is the margin multiplier.

#' @param nT sample size for test data
#' @param nR sample size for reference data
#' @param smplMuT sample mean for test data
#' @param smplMuR sample mean for reference data
#' @param smplSigmaT sample standard deviation for test data
#' @param smplSigmaR sample standard devivation for reference data
#' @param vecT a vector of observations for test product
#' @param vecR a vector of observations for reference product
#' @param eta the margin multipler
#' @return a list containing the RMLE for the means and standard deviations for both test and reference data
# @export

RMLE_equivTest <- function(nT,nR,smplMuT,smplMuR,smplSigmaT,smplSigmaR,vecT,vecR,eta)
{
  #set initial estimator
  muT = smplMuT
  muR = smplMuR

  sigmaT = smplSigmaT
  sigmaR = smplSigmaR

  oldPrmtr = c(muT,muR,sigmaT,sigmaR)

  #prepare iteration
  dif = 1
  iterMax = 200
  iter = 1
  #iteration
  while(dif>1.0e-4 & iter < iterMax)
  {
    #print(c(iter,oldPrmtr))
    muT = muR + eta*sigmaR
    sigmaT = sqrt( mean( (vecT - muT)^2 ) )

    muR = (nT*(smplMuT-eta*sigmaR)/sigmaT^2 + nR*smplMuR/sigmaR^2)/(nT/sigmaT^2 + nR/sigmaR^2)

    poly = polynom::polynomial(c(-mean( (vecR-muR)^2), 0, 1,-nT/nR*eta*(smplMuT-muR)/sigmaT^2, nT/nR*eta^2/sigmaT^2 ))
    pz = solve(poly)
    #print(pz)
    i = 4
    rootFound = FALSE
    while(i >=1)
    {
      if(Im(pz[i])==0 & Re(pz[i])>0)
      {
       sigmaR = Re(pz[i])
       rootFound = TRUE
       break
      }
      i = i - 1
    }
    if(!rootFound)
      print("Cannot find real root for sigmaR")

    newPrmtr = c(muT,muR,sigmaT,sigmaR)
    if( mean( abs(oldPrmtr - newPrmtr)) < 1.0e-4)
      break
    oldPrmtr = newPrmtr
    iter = iter + 1
  }

  list(muT = muT,
       muR = muR,
       sigmaT = sigmaT,
       sigmaR = sigmaR
       )
}



#' Equivalence test by Modified Wald test with standard error estimated by  RMLE (MWCMLE)
#'
#' Equivalence test by Modified Wald test with standard error estimated by  RMLE (MWCMLE).
#'
#'See \insertCite{weng2018improved;textual}{equivalenceTest}.
#'
# Test statistics:
# \deqn{a+b \hat{\mu}_T - \hat{\mu}_R+f\hat{\sigma}_R}
#
# \deqn{W_L=\frac{\hat{\mu}_T - \hat{\mu}_R+f\hat{\sigma}_R}{\sqrt{\frac{\tilde{\sigma}^2_{T,L}}{n_T} +\left(\frac{1}{n_R}+\frac{f^2 V_{nR}}{n_R-1}\right)\tilde{\sigma}_{RL})^2 }}
# }
#
#' @inheritParams equivTestFixedMargin
#' @return a list containing the test result
#' @export
#'
#' @references
#' {
#' \insertRef{weng2018improved}{equivalenceTest}
#' }
equivTestMWCMLE <- function(vecT,vecR,alpha=0.05,marginX=1.7,method="MWCMLE")
{
  eta = marginX
  nT = length(vecT)
  nR = length(vecR)

  smplMuT = mean(vecT)
  smplMuR = mean(vecR)

  smplSigmaT = sd(vecT)
  smplSigmaR = sd(vecR)

  k = sqrt((nR-1)/2)*exp(lgamma((nR-1)/2) - lgamma(nR/2))

  est = RMLE_equivTest(nT,nR,smplMuT,smplMuR,smplSigmaT,smplSigmaR,vecT,vecR,-eta)
  tstatL = (smplMuT - smplMuR + eta*smplSigmaR)/sqrt(est$sigmaT^2/nT + (1.0/nR + eta^2*(1-1/k^2))*est$sigmaR^2)
  estL = est

  est = RMLE_equivTest(nT,nR,smplMuT,smplMuR,smplSigmaT,smplSigmaR,vecT,vecR,eta)
  tstatU = (smplMuT - smplMuR - eta*smplSigmaR)/sqrt(est$sigmaT^2/nT + (1.0/nR + eta^2*(1-1/k^2))*est$sigmaR^2)
  estU = est

  qntl = qnorm(alpha,lower.tail=F)

  rslt = ifelse(tstatL>qntl & tstatU< -qntl,1,0)
  ret = list(method=method,
             alpha = alpha,
             estL=estL,
             estU=estU,
             tstatL = tstatL,
             tstatU = tstatU,
             qntl = qntl,
             rslt = rslt)
  invisible(ret)
}



