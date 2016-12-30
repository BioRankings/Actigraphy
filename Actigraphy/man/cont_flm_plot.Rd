\name{cont_flm_plot}
\alias{cont_flm_plot}
\title{Plot Functional Linear Model Analysis Results of a Continuous Type}

\description{This function produces two plots: An effect of a continuous covariate on activity values by time and the F-test for the effect of the continuous covariate.}

\usage{cont_flm_plot(smoothdata, matchresults, flmresults, xlim, ylim, 
	ftest, nperm, lb, xat, legendx, legendy, L, xlab="Time", ylab="Activity")
}

\arguments{
  \item{smoothdata}{List output from the \code{fda.smoothdata} function.}
  \item{matchresults}{List output from the \code{matchid} function.}
  \item{flmresults}{List output from the \code{flm_cate} function.}
  \item{xlim}{X axis limits for activity plot.}
  \item{ylim}{Y axis limits for activity plot.}
  \item{ftest}{A logic value indicating whether to implement F test or not. F test will be implement if \code{ftest is TRUE}.}
  \item{nperm}{The number of permutations for the F-test.}
  \item{lb}{X-axis labels.}
  \item{xat}{X axis label positions.}
  \item{legendx}{X axis position of the left edge of the legend box.}
  \item{legendy}{Y axis position of the upper edge of the legend box.}
  \item{L}{The length of the time points.}
  \item{xlab}{The label for the x-axis.}
  \item{ylab}{The label for the y-axis.}
}

\value{One plot of the estimated group means and a possible second plot of the F-test results.}

\author{William D. Shannon, Tao Li, Hong Xian, Jia Wang, Elena Deych, Carlos Gonzalez}

\examples{
	data(act_29pt)
	data(clinic_29pt_ahi)
	
	colnames(act_29pt) <- sub("X", "", colnames(act_29pt)) 
	data <- as.matrix(act_29pt[,-1])
	
	matchid  <- fda.matchid(data, clinic_29pt_ahi, "contin")
	FDcont <- fda.smoothdata(matchid)
	
	L <- nrow(data)
	lb <- c("Midnight", "6AM", "Noon", "6PM", "Midnight") 
	xat <- c(0, L/4, L/2, 3*L/4, L)
	
	geftFDcont <- flm_cate(FDcont)
	predy <- as.vector(geftFDcont$freg$yhatfdobj$y)
	
	xlim <- c(0, L) 
	ylim <-  c(min(predy), max(predy) + 100)
	
	legendx <- 0
	legendy <- max(predy) - 100
	
	cont.flm.results <- cont_flm_plot(FDcont, matchid, geftFDcont, xlim, 
		ylim, TRUE, 10, lb, xat, legendx, legendy, L)
}