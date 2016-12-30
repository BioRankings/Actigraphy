\name{cat_flm_plot}
\alias{cat_flm_plot}
\title{Plot Functional Linear Model Analysis Results of a Categorical Type}

\description{
	This function produce either one or two plots: An effect of a categorical (factor) covariate on activity 
	values by time and potentially the F-test for the effect of the categorical covariate.
}

\usage{cat_flm_plot(smoothdata, matchresults, flmresults, ftest, 
	nperm, lb, xat, varname, col, ylim, L, xlab="Time", ylab="Activity")
}

\arguments{
  \item{smoothdata}{List output from the \code{fda.smoothdata} function}
  \item{matchresults}{List output from the \code{matchid} function.}
  \item{flmresults}{List output from the \code{flm_cate} function.}
  \item{ftest}{A logic value indicating whether to implement F test or not. F test will be implement if \code{ftest is TRUE}.}
  \item{nperm}{The number of permutations for the F-test.}
  \item{lb}{X axis labels.}
  \item{xat}{X axis label positions.}
  \item{varname}{Name of categorical covariate.}
  \item{col}{Colors for levels of the predictor.}
  \item{ylim}{Y axis limits for activity plot.}
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
	ahi <- clinic_29pt_ahi
	
	ahi$ahicat <- as.factor(ifelse(ahi$AHI >= 0 & ahi$AHI <= 5, 1, 
					ifelse(ahi$AHI > 5 & ahi$AHI <= 15, 2,
					ifelse(ahi$AHI > 15 & ahi$AHI <= 30, 3,
					ifelse(ahi$AHI > 30, 4, 0)))))
	
	matchidb <- fda.matchid(data, ahi[,-2] , "factor", 
		c("normal", "mild", "moderate", "severe"))
	FDcatahi <- fda.smoothdata(matchidb)
	
	L <- nrow(data)
	lb <- c("Midnight", "6AM", "Noon", "6PM", "Midnight") 
	xat <- c(0, L/4, L/2, 3*L/4, L)
	
	geftFDcatahi <- flm_cate(FDcatahi)
	predy <- as.vector(geftFDcatahi$freg$yhatfdobj$y)
	
	xlim <- c(0, L) 
	ylim <- c(min(predy), max(predy) + 100)
	
	cat.flm.results <- cat_flm_plot(FDcatahi, matchidb, geftFDcatahi, 
		TRUE, 5, lb, xat, "AHI", 1:4, ylim, L)
}



