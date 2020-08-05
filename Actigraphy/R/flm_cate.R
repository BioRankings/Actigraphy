flm_cate <-
function(FD, basistype="fourier", nbasis=9, norder=4){
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
