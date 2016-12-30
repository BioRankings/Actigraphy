\name{flm_cate}
\alias{flm_cate}
\title{Functional Linear Model Analysis}

\description{A function that does functional linear model analysis.}

\usage{flm_cate(FD, basistype="fourier", nbasis=9, norder=4)}

\arguments{
  \item{FD}{The list from the function \code{fda.smoothdata}.}
  \item{basistype}{A string specifying either \code{"fourier"} and \code{"bspline"}.}
  \item{nbasis}{The number of basis functions to be used for functional linear model analysis.  Default value is 9.}
  \item{norder}{an integer specifying the order of b-splines, which is one higher than their degree. The default of 4 gives cubic splines.}
}

\value{
	A list consisting of three components as follows:
		\item{freg}{A \code{fRegress} fit object containing the intercept and coefficient functions (check function \code{fRegress} for details)}
		\item{fregstd}{A list containing the standard error functions of the intercept and coefficient functions.} 
}

\author{William D. Shannon, Tao Li, Hong Xian, Jia Wang, Elena Deych, Carlos Gonzalez}

\examples{
	data(act_29pt)
	data(clinic_29pt_ahi)
	
	colnames(act_29pt) <- sub("X", "", colnames(act_29pt)) 
	data <- as.matrix(act_29pt[,-1])
	
	matchid  <- fda.matchid(data, clinic_29pt_ahi, "contin")
	FDcont <- fda.smoothdata(matchid)
	
	geftFDcont <- flm_cate(FDcont)
}


