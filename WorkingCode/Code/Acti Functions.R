
library(fda)		# lots of uses



fda.smoothdata <- function(data, basistype="fourier", nbasis=9, norder=4){
	if(missing(data)) 
		stop("Missing data")
	
	mat <- data$mat
	cov <- data$cov
	L <- nrow(mat)
	
	if(tolower(basistype) == "fourier"){
		fbase <- fda::create.fourier.basis(rangeval=c(0, L), nbasis)
	}else if(tolower(basistype) == "bspline"){
		fbase <- fda::create.bspline.basis(rangeval=c(0, L), nbasis, norder)
	}else{
		stop("basistype must be 'fourier' or 'bspline'.")
	}
	
	fpar <- fda::fdPar(fbase) 
	fd <- fda::smooth.basis(1:L, mat, fpar)
	FD <- list(fd=fd, cov=cov)
	return(FD)
}

fda.matchid <- function(mat, acov, type, grouplab){
	if(missing(mat) || missing(acov) || missing(type)) 
		stop("Missing Arguments")
	if(class(mat) != "data.frame" && class(mat) != "matrix") 
		stop("mat must be a data.frame or a matrix")
	if(class(acov) != "data.frame" && class(acov) != "matrix") 
		stop("acov must be a data.frame or a matrix")
	
	coln <- colnames(mat)
	for(i in 1:length(coln)){
		if(substr(coln[i], 1, 1) == "X"){
			coln <- sub("X", "", coln)
			colnames(mat) <- coln
		}
	}
	numcolnames <- coln
	
	acov <- acov[!is.na(acov[[2]]),]
	id <- acov[[1]]
	bthnum <- is.numeric(id) + is.numeric(coln)
	
	if(bthnum == 1){
		commonid <- intersect(as.numeric(id), as.numeric(numcolnames))	
		sltmat <- as.matrix(mat[, charmatch(as.numeric(commonid), as.numeric(numcolnames))])
	}else{
		commonid <- intersect(as.character(id), as.character(numcolnames))
		sltmat <- as.matrix(mat[, charmatch(commonid, numcolnames)])
	}
	acov <- acov[charmatch(commonid, id),]
	
	if(length(commonid) < 1)
		stop("IDs do not match")
	
	a <- acov
	if(tolower(type) == "factor"){
		a[[2]] <- as.factor(a[[2]])
		m <- data.frame(id=a[[1]], stats::model.matrix(~a[[2]]))	
		
		if(missing(grouplab))
			stop("grouplab is required for 'factor' type.")
		colnames(m) <- c("id", grouplab)
	}else if(tolower(type) == "contin"){
		a[[2]] <- as.numeric(a[[2]])
		m <- data.frame(id=a[[1]], stats::model.matrix(~a[[2]]))
		contcovname <- names(a[2])
		colnames(m) <- c("id", "intercept", paste(contcovname))
	}else{
		stop("type must be 'factor' or 'contin'.")
	}

	return(list(mat=sltmat, cov=m))
}

flm_cate <- function(FD, basistype="fourier", nbasis=9, norder=4){
    if (missing(FD)) 
        stop("Missing FD")

	cov <- FD$cov[,-1]
	grp <- ncol(cov)
	fd <- FD$fd
	L <- length(fd$argvals)
	npt <- ncol(fd$y)
	
	if(tolower(basistype) == "fourier"){
		fbase <- fda::create.fourier.basis(rangeval=c(0, L), nbasis)
	}else if(tolower(basistype) == "bspline"){
		fbase <- fda::create.bspline.basis(rangeval=c(0, L), nbasis, norder)
	}else{
		stop("basistype must be 'fourier' or 'bspline'.")
	}
	
	fpar <- fda::fdPar(fbase)
	
	xfdlist <- vector("list", grp)
	xfdlist[[1]] <- cov[,1] + 0
	
	for(i in 2:grp) 
		xfdlist[[i]] <- cov[,i] + 0
	betalist <- xfdlist
	
	for(i in 1:grp) 
		betalist[[i]] <- fpar 
	
	freg2 <- fda::fRegress(fd$fd, xfdlist, betalist)
	
	preact2 <- stats::predict(freg2$yhatfdobj, 1:L)
	resid2 <- fd$y - preact2[, 1:npt]
	sigma2 <- cov(t(resid2))
	fregstd2 <- fda::fRegress.stderr(freg2, fd$y2cMap, sigma2)
	
    return(list(freg=freg2, fregstd=fregstd2))
}

cont_flm_plot <- function(smoothdata, matchresults, flmresults, xlim, ylim, ftest, 
		nperm, lb, xat, legendx, legendy, L, xlab="Time", ylab="Activity"){

	if(missing(smoothdata) || missing(matchresults) || missing(flmresults) || missing(xlim) || 
			missing(ylim) || missing(ftest) || missing(nperm) || missing(lb) || missing(xat) || 
			missing(legendx) || missing(legendy) || missing(L)) 
		stop("Missing Arguments")
	
	colort <- factor(matchresults$cov[,3])
	ucont <- length(unique(colort))
	covname <- names(matchresults$cov[3])
	maintitle <- paste(ylab, "~", covname, sep="")
	LofLegend <-ifelse(length(levels(colort)) <= 100, 10, 1)
	
	minmedmax <- function(contvar){
		contvar <- contvar
		mincont <- signif(min(contvar), 3)
		medcont <- signif(stats::median(contvar), 3)
		maxcont <- signif(max(contvar), 3)
		contlab <- list(mincont, medcont, maxcont)
		return(contlab)
	}
	
	contlabels <- minmedmax(matchresults$cov[,3])

	par(mfrow=c(2,1), mar=c(4,4,3,1))
	plot(0, 0, xlim=xlim, ylim=ylim, xaxt="n", xlab=xlab, ylab=ylab, type='n', main=maintitle)
	for(i in 1:length(colort)) 
		lines(predict(flmresults$freg$yhatfdobj, c(1:L))[,i], col=grDevices::topo.colors(ucont)[colort[i]])
        
	axis(1, at=xat, labels=lb)
	
	legend("topleft", legend=c(contlabels[[1]], "",  "", "", contlabels[[3]]), pch=15,  title="Range of Values",
			col=grDevices::topo.colors(ucont)[c(1, round(ucont/4), round(ucont/2), round(ucont*3/4), ucont)])

	geftFtestresults <- flm_ftest(smoothdata, "fourier", smoothdata$fd$fd$basis$nbasis, 4, ftest, nperm, xat, lb, 1.5)
        
	return(geftFtestresults)
}

cat_flm_plot <- function(smoothdata, matchresults, flmresults, ftest, nperm, lb, xat, 
		varname, col, ylim, L, xlab="Time", ylab="Activity"){
	
	if(missing(smoothdata) || missing(matchresults) || missing(flmresults) || missing(ftest) || 
			missing(nperm) || missing(lb) || missing(varname) || missing(col) || missing(ylim) || missing(L)) 
		stop("Missing Arguments")
	
	geft <- flmresults
	maintitle<- paste(ylab, "~", varname, sep="")
	leg <- names(matchresults$cov)[-1]
	l <- length(geft$freg$betaestlist)
	beta <- geft$freg$betaestlist
	std <- geft$fregstd$betastderrlist
	
	par(mfrow=c(2,1), mar=c(4,4,3,1))
	plot(0, 0, xlim=c(0, L), ylim=ylim, xlab=xlab, ylab=ylab, type="n", main=maintitle, xaxt="n", cex.main=0.8)
	lines(beta[[1]]$fd, col=col[1], lwd=2)
	for(i in 2:l) 
		lines(beta[[1]]$fd + beta[[i]]$fd, col=col[i], lwd=2)
	
	axis(1, at=xat, labels=lb)
	legend("topleft", leg, col=col, lty=rep(1, length(beta)), cex=0.8, lwd=2)
	
	geftFtestresults <- flm_ftest(smoothdata, smoothdata$fd$fd$basis[[2]], smoothdata$fd$fd$basis$nbasis, 4, ftest, nperm, xat, lb, 1.5)
	
	return(geftFtestresults)
}

flm_ftest <- function(FD, basistype="fourier", nbasis=9, norder=4, ftest=TRUE, nperm=1000, xat, lb, mul=1){ 
	if(missing(FD)) 
		stop("Missing FD")
	
	cov <- FD$cov[,-1]
	grp <- ncol(cov)
	fd <- FD$fd
	L <- length(fd$argvals)
	npt <- ncol(fd$y)
	
	if(tolower(basistype) == "fourier"){
		fbase <- fda::create.fourier.basis(rangeval=c(0, L), nbasis)
	}else if(tolower(basistype) == "bspline"){
		fbase <- fda::create.bspline.basis(rangeval=c(0, L), nbasis, norder)
	}else{
		stop("basistype must be 'fourier' or 'bspline'.")
	}
	
	fpar <- fda::fdPar(fbase)
	
	xfdlist <- vector("list",grp)
	xfdlist[[1]] <- cov[,1]+0
	for(i in 2:grp) 
		xfdlist[[i]] <- cov[,i]+0
	betalist <- xfdlist
	
	for(i in 1:grp) 
		betalist[[i]] <- fpar 
	
	freg2 <- fda::fRegress(fd$fd, xfdlist, betalist)
	
	preact2 <- stats::predict(freg2$yhatfdobj, c(1:L))
	resid2 <- fd$y - preact2[,1:npt]
	sigma2 <- cov(t(resid2))
	fregstd2 <- fda::fRegress.stderr(freg2, fd$y2cMap, sigma2)
	
	Fratio <- NULL
	if(ftest){
		Fratio <- Ftest(fd$fd, xfdlist, betalist, NULL, nperm, 1:L, xaxt="n", mul=mul) 
		axis(1, at=xat, labels=lb)
	}
	
	return(list(freg=freg2, fregstd=fregstd2, Fratio=Fratio))
}

Ftest <- function(yfdPar, xfdlist, betalist, wt=NULL, nperm=1000, argvals=NULL, q=0.05, plotres=TRUE, mul=1, ...){
	if(missing(yfdPar) || missing(xfdlist) || missing(betalist)) 
		stop("Missing Arguments")
	
	set.seed <- 123456789
	Fnull <- rep(0, nperm)
	Fnullvals <- c()
	q <- 1 - q
	begin <- proc.time()
	fRegressList <- fda::fRegress(yfdPar, xfdlist, betalist)
	elapsed.time <- max(proc.time() - begin, na.rm=TRUE)     
	print(paste("Permutation F test running (", nperm, " permutations)", sep=""))
	print(paste("Estimated Computing time = ", round(nperm * elapsed.time), " seconds", sep=""))
	
	yhat <- fRegressList$yhatfdobj
	if(is.list(yhat) && ("fd" %in% names(yhat))) 
		yhat <- yhat$fd
	tFstat <- fda::Fstat.fd(yfdPar, yhat, argvals)
	Fvals <- tFstat$F
	Fobs <- max(Fvals)
	argvals <- tFstat$argvals
	
	if(is.vector(yfdPar)){
		n <- length(yfdPar)
	}else{
		n <- ncol(yfdPar$coefs)
	}
	
	for(i in 1:nperm){
		tyfdPar <- yfdPar[sample(n)]
		fRegressList <- fda::fRegress(tyfdPar, xfdlist, betalist)
		yhat <- fRegressList$yhatfdobj
		if(is.list(yhat) && ("fd" %in% names(yhat))) 
			yhat <- yhat$fd
		tFstat <- fda::Fstat.fd(yfdPar, yhat, argvals)
		Fnullvals <- cbind(Fnullvals, tFstat$F)
		Fnull[i] <- max(Fnullvals[, i])
	}
	
	pval <- mean(Fobs < Fnull)
	qval <- stats::quantile(Fnull, q)
	pvals.pts <- apply(Fvals < Fnullvals, 1, mean)
	qvals.pts <- apply(Fnullvals, 1, quantile, q)
	if(plotres){
		if(is.fd(yfdPar)){
			ylims <- c(min(c(Fvals, qval, qvals.pts)), max(c(Fobs, qval))*mul)
			
			if(is.null(names(yhat$fdnames))){
				xlab <- "argvals"
			}else{
				xlab <- names(yhat$fdnames)[1]
			}
			
			plot(argvals, Fvals, type="l", ylim=ylims, col=2, lwd=2, xlab=" ", ylab="F-statistic", main="Permutation F-Test", ...)
			lines(argvals, qvals.pts, lty=3, col=4, lwd=2)
			abline(h=qval, lty=2, col=4, lwd=2)
			legendstr <- c("Observed Statistic", paste("pointwise", 1 - q, "critical value"), paste("maximum", 1 - q, "critical value"))
			legend("topleft", legend=legendstr, col=c(2,4,4), lty=c(1, 3, 2), lwd=c(2,2,2), cex=0.8)
		}else{
			xlims <- c(min(c(Fnull, Fobs)), max(c(Fnull, Fobs)))
			hstat <- hist(Fnull, xlim=xlims, lwd=2, xlab="F-value", main="Permutation F-Test")
			abline(v=Fobs, col=1, lwd=2)
			abline(v=qval, col=1, lty=2, lwd=2)
			legendstr <- c("Observed Statistic", paste("Permutation", 1 - q, "critical value"))
			legend("topleft", legend=legendstr, lty=c(1, 2), lwd=c(2,2))
		}
	}
	
	return(list(pval=pval, qval=qval, Fobs=Fobs, Fnull=Fnull, Fvals=Fvals, 
			Fnullvals=Fnullvals, pvals.pts=pvals.pts, qvals.pts=qvals.pts, 
			fRegressList=fRegressList, argvals=argvals))
}

