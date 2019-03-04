#author: Chao Wang, (wan9c9@gmail.com)

#' Compute the Satterthwaite approximation of degree of freedom for t distribution
#'
#' Compute the Satterthwaite approximation of degree of freedom for t distribution.
#'
#' @param s1 sample standard deviation for group 1
#' @param n1 sample size for group 1
#' @param n1s adjusted sample size for group 1
#' @param s2 sample standard deviation  for group 2
#' @param n2 sample size for group 2
#' @param n2s adjusted sample size for group 2
#' @return degree of freedom
dfSatterthwaite <- function(s1,n1,n1s,
                            s2,n2,n2s)
{
  invisible((s1^2/n1s + s2^2/n2s)^2/((s1^2/n1s)^2/(n1-1) + (s2^2/n2s)^2/(n2-1)))
}

#' Create summary information of a dataset
#'
#' Create a list of summary statistics of a dataset for equivalence test.
#'
#' @param smpl a vector representing the dataset
#' @return a list of objects summarizing the dataset
#' @export
#' @examples
#' vecT = rnorm(n=20)
#' s = createEquivTestSmpl(vecT)
#'
createEquivTestSmpl <- function(smpl)
{
  if(is.vector(smpl))
  ret = list(obs = smpl,
       nL = length(smpl),
       muLot = smpl,
       mu = mean(smpl),
       sigma = sd(smpl),
       CV = sd(smpl)/mean(smpl),
       min = min(smpl),
       max = max(smpl))
  invisible(ret)
}


#' Conduct the equivalence test with fixed margin
#'
#' Conduct the equivalence test with fixed margin.
#'
#' @param vecT the sample data for test product, can be a vector of observed values or a list returned by \code{createEquivTestSmpl}
#' @param vecR the sample data for reference product, can be a vector of observed values or a list returned by \code{createEquivTestSmpl}
#' @param alpha the nominal size, default = 0.05
#' @param marginX the margin multiplier, default = 1.5
#' @param sampleSizeX the sample size adjustment coefficient, default = 1.5
#' @param qa a string representing the name of the quality attribute, default = ""
#' @param sigmaTOverride a numeric value to override the estimate for standard deviation of the test product
#' @param labelT the name of the test product, default = "Proposed"
#' @param labelR the name of the reference product, default = "Reference"
#' @param show.message a logic value indicating whether messages are to be shown, default = FALSE
#' @param method a string indicating the method used in the equivalence test.
#' @return a list of objects summarizing the data and test results, in particular, \code{rslt} = 1 if  \eqn{H_0} is rejected, and \code{rslt} = 0 if \eqn{H_0} is not rejected.
#' @references
#' {
#' \insertRef{tsong2017development}{equivalenceTest}
#' }
#' @importFrom Rdpack reprompt
#' @export
#' @examples
#' vecT = rnorm(20,-1.5,1)
#' vecR = rnorm(20,0,1)
#' et = equivTestFixedMargin(vecT, vecR)
#'
equivTestFixedMargin <- function(vecT,vecR,alpha=0.05,marginX=1.5,sampleSizeX=1.5,qa="",sigmaTOverride=NULL,labelT="Proposed",labelR="Reference",show.message=FALSE,method="Fixed Margin"){
  smplT = vecT
  smplR = vecR
  if(is.numeric(smplT))
  {
    smplT = createEquivTestSmpl(smplT)
  }
  if(is.numeric(smplR))
  {
    smplR = createEquivTestSmpl(smplR)
  }

  if(!is.null(sigmaTOverride))
  {
    sigmaTOriginal = smplT$sigma
    smplT$sigma = sigmaTOverride
  }

  margin = marginX*smplR$sigma
  dm = smplT$mu - smplR$mu
  #adjust lot sample size if unbalanced
  nRLs = ifelse(smplR$nL>sampleSizeX*smplT$nL,sampleSizeX*smplT$nL,smplR$nL)
  nTLs = ifelse(smplT$nL>sampleSizeX*smplR$nL,sampleSizeX*smplR$nL,smplT$nL)
  dsigma = sqrt(smplR$sigma^2/nRLs + smplT$sigma^2/nTLs)
  tdf =  dfSatterthwaite(smplT$sigma,smplT$nL,nTLs,
                         smplR$sigma,smplR$nL,nRLs)
  tq = qt(p=1-alpha, df = tdf, lower.tail= TRUE)

  #else
  #{
  #  sp = ((nRL-1)*smplR$sigma^2 + (nTL-1)*smplT$sigma^2)/(nRL+nTL-2)
  #
  #  dsigma = sqrt(smplR$sigma^2/nRL + smplT$sigma^2/nTL)
  #  tq = qt(p=1-alpha,df = nRL+nTL-2,lower.tail = FALSE)
  #}

  ci = c(dm-tq*dsigma,dm+tq*dsigma)

  rslt=list(method = method,
            smplT=smplT,
            smplR=smplR,
            labelT=labelT,
            labelR=labelR,
            meanDif=dm,
            ci = ci,
            qa = qa,
            margin=c(-margin,margin),
            rslt = ifelse( all(abs(ci)/margin <= 1), 1,0) # 1 reject H0: not equivalent
            )
  twoRowValue <- function(x){paste0("\\multirow{2}{*}{",x,"}")}
  intervalPrint <- function(a,b){paste0("(",a,",",b,")")}
  rslt$etRslt = data.frame(
                      Product=c(labelT,labelR),
                      "# of Lots"=c(smplT$nL,smplR$nL),
                      "Range"=c(paste0("(",smplT$min,"-",smplT$max,")"),
                                paste0("(",smplR$min,"-",smplR$max,")")
                                ),
                      "Mean"=c(smplT$mu,smplR$mu),
                      "Std. Dev."=c(smplT$sigma,smplR$sigma),
                      #"CV(%)" = c(100*smplT$CV,100*smplR$CV),
                      "Mean Diff. (90% CI)" = rep(paste0(round(dm,2)," (",round(ci[1],2),", ",round(ci[2],2),")"),2),
                      "Margin" = rep(intervalPrint(round(-margin,2),round(margin,2)),2),
  "Pass ET?"=rep(ifelse(all(abs(ci)/margin <= 1),"Yes","No"),2),
  check.names = FALSE
  )

  rslt$latexTable =paste(c(
    "\\begin{table}[H]
     \\centering
    \\begin{tabular}{ccccccccc}\n",
    paste(c("Drug","\\# of Lots", "Range", "Mean", "Std. Dev.", "CV(\\%)","Mean Diff. (90\\% CI)","Equiv. Margin","Pass Equiv. Test?"),collapse=" & "),
    "\\\\\n",
    paste0(c(labelT,smplT$nL,paste0("(",smplT$min,"-",smplT$max,")"),round(smplT$mu,2),round(smplT$sigma,2),round(100*smplT$CV,2),
             twoRowValue(paste0(round(dm,2)," (",round(ci[1],2),round(ci[2],2),")")),twoRowValue(intervalPrint(round(-margin,2),round(margin,2))),twoRowValue(ifelse(rslt$rslt,"Yes","No"))),collapse=" & "),"\\\\\n",
    paste(c(labelR,smplR$nL,paste0("(",smplR$min,"-",smplR$max,")"),round(smplR$mu,2),round(smplR$sigma,2),round(100*smplR$CV,2), " ", " ", " "),collapse=" & "),
    "\n",
    "\\end{tabular}
    \\end{table}")
  )

  if(show.message)
  { message("Summary of data")
    if(qa!="")
      message("QA: ",qa)
    message("     Test: #lot: ",smplT$nL, ", mean: ",round(smplT$mu,2), ifelse(is.null(sigmaTOverride),paste0(", sigma: ", round(smplT$sigma,2)),
                                                                             paste0(", sigma (override): ", round(smplT$sigma,2),", sigma (sample est): ", round(sigmaTOriginal,2))))

    message("Reference: #lot: ",smplR$nL, ", mean: ",round(smplR$mu,2), ", sigma: ", round(smplR$sigma,2))


    message("diff in mu: ", round(dm,2), ", est'd sigma: ",round(dsigma,2), ", df: ", round(tdf,2), ", tq:", round(tq,2))
    message(round(100*(1-alpha*2),1),"% CI: [",round(ci[1],2),", ", round(ci[2],2), "]")
    message("Margin: ", round(margin,2))
    message("The null hypothesis that the difference between the means of two samples are not within a margin is ",ifelse(rslt$rslt==1,""," not "), "rejected.")
  }
  #cat(rslt$latexTable)
  invisible(rslt)
  }

