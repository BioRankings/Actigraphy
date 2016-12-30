\name{fda.matchid}
\alias{fda.matchid}
\title{Match IDs from Clinical and Actigraph Data}

\description{
	A function used to match actigraphy data and clinical covariates by subject IDs and 
	return a list of the data combined by IDs. Only the subjects with both actigraphy and covariate 
	data will be returned by this function.
}

\usage{fda.matchid(mat, acov, type, grouplab)}

\arguments{
  \item{mat}{A data frame with the rows being the time and the columns being the activity, with the column names being the subjects.}
  \item{acov}{A two column data frame that contains only subject IDs and a covariate of interest, respectively.} 
  \item{type}{A string specifying either \code{"contin"} for continuous and \code{"factor"} for categorical covariates.}
  \item{grouplab}{A vector of names of the categories if \code{type is TRUE}.} 
}

\details{
	Note: Only the subjects with both actigraphy and covariate data will be returned by this function.
}

\value{
	A list consisting of two components as follows:
	\item{mat}{A matrix where rows represent the time, columns are the samples, and the 
		column names are the subjects.}
	\item{cov}{A two column matrix that contains the actigraphy data and clinical covariates.} 
}

\author{William D. Shannon, Tao Li, Hong Xian, Jia Wang, Elena Deych, Carlos Gonzalez}

\examples{
	data(act_29pt)
	data(clinic_29pt_ahi)
	
	colnames(act_29pt) <- sub("X", "", colnames(act_29pt)) 
	data <- as.matrix(act_29pt[,-1])
	
	### Example 1: Continuous Covariate
	matchida  <- fda.matchid(data, clinic_29pt_ahi, "contin")
	
	
	### Example 2: Categorical Covariate
	ahi <- clinic_29pt_ahi
	
	ahi$ahicat <- as.factor(ifelse(ahi$AHI >= 0 & ahi$AHI <= 5, 1, 
	       ifelse(ahi$AHI > 5 & ahi$AHI <= 15, 2,
	       ifelse(ahi$AHI > 15 & ahi$AHI <= 30, 3,
	       ifelse(ahi$AHI > 30, 4, 0)))))
	
	matchidb  <- fda.matchid(data, ahi[,-2], "factor", 
		c("normal", "mild", "moderate", "severe"))
}