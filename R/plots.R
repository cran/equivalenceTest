#' Provide a side-by-side scatter plot of two or three datasets for equivalence test.
#'
#' Provide a side-by-side scatter plot of two samples for equivalence test.
#' @param vecT a vector of the sample for test product
#' @param vecR a vector of the sample for reference product
#' @inheritParams equivTestFixedMargin
#' @param vecR1 a vector of the sample for reference product R1
#' @param labelR1 label for reference product R1
#' @export
#' @return NULL
#' @examples
#' vecT = rnorm(20,-1.5,1)
#' vecR = rnorm(20,0,1)
#' vecR1 = rnorm(20,0,1)
#' scatterPlotEquivTestData(vecT,vecR,labelT="T",labelR="R",qa="potency")
#' scatterPlotEquivTestData(vecT,vecR,vecR1,labelT="T",labelR="R",labelR1="R1",qa="potency")

scatterPlotEquivTestData <- function(vecT,vecR,vecR1=NULL,qa="",labelT="Test",labelR="Reference",labelR1="Reference1")
{
  n1 = length(vecT)
  n2 = length(vecR)

  idx1 = seq(0.5,1.5,length.out=n1)
  idx2 = seq(2.5,3.5,length.out = n2)
  if(!is.null(vecR1))
  {
    n3 = length(vecR1)
    idx3 = seq(4.5,5.5,length.out = n3)
    yrange = range(c(vecT,vecR,vecR1))
  }else{
    yrange = range(c(vecT,vecR))
  }

  plot(idx1,vecT,col="red",xlim=c(0,ifelse(is.null(vecR1),4,6)),ylim=yrange,xaxt='n',xlab="",ylab=qa,pch=16,cex=2,cex.axis=1.5,cex.lab=1.5)
  points(idx2,vecR,col="blue",pch=18,cex=2)
  if(!is.null(vecR1))
  {  points(idx3,vecR1,col="black",pch=18,cex=2)
     axis(side = 1,at = c(1,3,5),labels = c(labelT,labelR,labelR1),cex.axis=1.5,cex.lab=1.5)}else{
     axis(side = 1,at = c(1,3),labels = c(labelT,labelR),cex.axis=1.5,cex.lab=1.5)
  }
  invisible()
}



#' Histogram with a fitted normal density function
#'
#' Provide a histogram with a fitted normal density.
#' @param x the data
#' @param main the title of the plot
#' @inheritParams stats::plot
#' @return NULL
#' @export
#' @examples
#' x = rnorm(20)
#' histWNormDensity(x)
#'
histWNormDensity <- function(x,main="")
{
  h <- hist(x, breaks=length(x)*0.4, main=main,density=10, xlim = range(x), col="lightgray", xlab="")
  xfit<-seq(min(x),max(x),length=40)
  yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
  yfit <- yfit*diff(h$mids[1:2])*length(x)
  lines(xfit, yfit, col="black", lwd=2)
  invisible()
}


#' Plot the equivalence test result
#'
#' Plot the equivalence test result including the margin, confidence intervals of the mean difference, and estimated mean difference.
#'
#' @param meanDif difference between mean of test and reference product
#' @param ci confidence interval for mean difference, a vector of two values
#' @param alpha nominal level of the hypothesis test
#' @param margin a vector consisting of lower margin and upper margin
#' @param qaNameLong the quality attribute name
#' @param testDrugName test drug name
#' @param refDrugName reference drug name
#' @param showDrugName logic value indicating if the drug names are to be shown.
#' @param showQA logic value indicating if the quality attribute (QA) is to be shown.
#' @param showCINumbers whether CI values are shown in the figure.
#' @return NULL
#' @export
#' @examples
#' equivTestPlot(0.623,c(-2,2),0.05,c(-9.79,9.79),
#'   "q a","test","reference")
#' equivTestPlot(0.623,c(-2,2),0.05,c(-9.79,9.79),
#'   "Relative Potency","test","reference",showDrugName = TRUE,showQA=TRUE,showCINumbers = TRUE)
#' equivTestPlot(0.5,c(-1.05,2.05),0.05,c(-9.79,9.79),
#'   "Relative Potency","test","reference",showQA=TRUE,showCINumbers = TRUE)

equivTestPlot <- function(meanDif,ci,alpha,margin, qaNameLong,testDrugName="",refDrugName="",showDrugName=FALSE,showQA=FALSE,showCINumbers=FALSE)
{

  plot(x = seq(1.2*margin[1],1.2*margin[2],length=10),y=rep(0.5,10),type="n",xlim=2*margin, ylim=c(0,1.5), col="lightgrey", xaxt='n', yaxt='n', xlab="", ylab="",bty="n")
  #title(paste0("Equivalence Test on ",qaNameLong))
  #rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "lightgrey")
  lines(rep(margin[1], 100), seq(0.3, 0.7, length=100),lty=1,col="black", lwd=5)
  lines(rep(margin[2], 100), seq(0.3, 0.7, length=100),lty=1,col="black", lwd=5)
  lines(seq(ci[1], ci[2], length=100), rep(0.5, 100), lwd=5, col="blue")
  points(meanDif, 0.5, pch=19, col="red",  cex=1.5, lwd=1.2)
  if(showDrugName)
    text((margin[1]+margin[2])/2, 1.1, paste0(testDrugName, " vs. ", refDrugName), font=2,cex=1.5)
  if(showQA)
    text((margin[1]+margin[2])/2, 0.9, paste0(qaNameLong), font=2,cex=1.5)
  text((ci[1]+ci[2])/2, 0.7, paste0(round(100*(1-2*alpha),0),"% CI", ifelse(showCINumbers,paste0("= \n [",round(ci[1],2),", ", round(ci[2],2),"]"),"")), font=2, col="blue",cex=1.5)
  text(margin[1], 0.15, round(margin[1],2),font=2, col="black")
  text(margin[2], 0.15, round(margin[2],2),font=2, col="black")
  invisible()
}


#' Provide a combined plot for equivvalence test
#'
#' Provide a combined plot for equivalence test, including both scatter plot of the sample data and a bar plot indicating the test result, where the null hypothesis is rejected if the red line representing the mean value of the test product lies within a grey rectangle centered at a blue line representing the mean value of the reference product.
#' @param et the list returned by \code{equivTestFixedMargin}
#'
#' @return NULL
#' @export
#' @examples
#' vecR = rnorm(20,0,1)
#' vecT = rnorm(20,-1.5,1)
#' et = equivTestFixedMargin(vecT,vecR)
#' equivTestFixedMarginCombPlot(et)
equivTestFixedMarginCombPlot <- function(et)
{
  vecT = et$smplT$obs
  vecR = et$smplR$obs
  nR = length(vecR)
  nT = length(vecT)
  idxT = seq(0.5,1.5,length.out=nT)
  idxR = seq(2.5,3.5,length.out = nR)


  plot(idxT,vecT,col="red",xlim=c(0,4),ylim=range(c(vecR,vecT)),xaxt='n',xlab="",ylab=et$qa,pch=16,cex=2,cex.axis=1.5,cex.lab=1.5)
  points(idxR,vecR,col="blue",pch=18,cex=2)
  centerIdx = c(1.75,2.25)

  rect(centerIdx[1],et$smplT$mu-et$ci[1] + et$margin[1],
       centerIdx[2],et$smplT$mu-et$ci[2]+et$margin[2],col="grey")
  segments(centerIdx[1],et$smplT$mu,centerIdx[2],et$smplT$mu,col="red",lwd=2)
  segments(centerIdx[1],et$smplR$mu,centerIdx[2],et$smplR$mu,col="blue",lwd=2)

  #segments(centerIdx[1],et$smplT$mu-et$ci[1] + et$margin[1],centerIdx[2],et$smplT$mu-et$ci[1] + et$margin[1],col="grey",lty=2,lwd=2) # lower limit
  #segments(centerIdx[1],et$smplT$mu-et$ci[2]+et$margin[2],centerIdx[2],et$smplT$mu-et$ci[2]+et$margin[2],col="grey",lty=2,lwd=2) # upper limit for mu_T hat
  axis(side = 1,at = c(1,3),labels = c(et$labelProposed,et$labelReference),cex.axis=1.5,cex.lab=1.5)
}
