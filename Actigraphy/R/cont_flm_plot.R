cont_flm_plot <-
function(smoothdata, matchresults, flmresults, xlim, ylim, ftest, 
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
