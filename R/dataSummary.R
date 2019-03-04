#' Summarize data for equivalence test
#'
#' Summarize data for equivalence test, can be two datasets or three datasets.
#'
#' @param vecT vector of sample for T
#' @param labelT label for T
#' @param vecR vector of sample for R
#' @param labelR label for R
#' @param vecR1 vector of sample for R1
#' @param labelR1 label for R1
#' @return a data.frame consisting the sample size, min, max, mean, SD, and percentage coefficient of variation for the samples
#' @export
#'
#' @examples
#' vecT = rnorm(10,-1.5,1)
#' vecR = rnorm(10)
#' vecR1 = rnorm(15,1,2)
#' ss = summarizeSample(vecT,"T",vecR,"R",vecR1,"R1")
#print(xtable(s,caption = "Summary of data for ADCC.",digits = 2),include.rownames = FALSE)


summarizeSample <- function(vecT,labelT,
                            vecR,labelR,
                            vecR1=NULL,labelR1="")
{
  briefSummary <- function(x,label)
  {
  list("Product"=label,
       "n" = length(x),
       "Min" = min(x),
       "Max" = max(x),
       "Mean" = mean(x),
       "SD" = sd(x),
       "CV(%)" = 100*sd(x)/abs(mean(x)))
  }
  sT = briefSummary(vecT,labelT)
  sR = briefSummary(vecR,labelR)
  if(is.null(vecR1))
  {
    s =rbind.data.frame(sT,sR)

  }else{
    sR1 = briefSummary(vecR1,labelR1)
    s =rbind.data.frame(sT,sR,sR1)
  }
  colnames(s)[7] = "CV (%)"
  invisible(s)
}
