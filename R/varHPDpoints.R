varHPDpoints <-
function(model,factors="units",factors2=NULL,hpdtype="w",jitter=0,...){
	res1=res2=c(rep(0,length(model$VCV[,1])))
	names1=names2=c()
	allnames=names(posterior.mode(model$VCV))
	solnames=names(posterior.mode(model$Sol))
	genes=solnames[-grep(":",solnames)]
	genes=genes[grep("gene",genes)]
	genes=sub("gene","",genes)
	ngenes=length(genes)
	for (f in factors){		
		pattern=paste('^\\w+\\.',f,'$',sep="")
		f1=grep(pattern,allnames)
		res1=res1+model$VCV[,f1]
	}
	
	if (is.null(factors2)==FALSE) {
		for (f in factors2){		
			pattern=paste('^\\w+\\.',f,'$',sep="")
			f1=grep(pattern,allnames)
			res2=res2+model$VCV[,f1]
		}
		stat=res2/(res1+res2)
		means=apply(stat,2,mean)
		hpds=HPDinterval(stat)
	} else {
		res1=res2=c(rep(0,length(model$VCV[,1])))
		for (f in factors){		
			pattern=paste('^\\w+\\.',f,'$',sep="")
			f1=grep(pattern,allnames)
			res1=res1+model$VCV[,f1]
		}
		stat=log(sqrt(exp(res1)),2)
		means=apply(stat,2,mean)
		hpds=HPDinterval(stat)
	}
	points(c(1:ngenes)+jitter,means,...)
	if (hpdtype=="l") {
		lines(hpds[,1],lty=2,...)
		lines(hpds[,2],lty=2,...)
	}	else {
		arrows(c(1:ngenes)+jitter,means,c(1:ngenes)+jitter,hpds[,1],angle=90,code=2,length=0.03,...)
		arrows(c(1:ngenes)+jitter,means,c(1:ngenes)+jitter,hpds[,2],angle=90,code=2,length=0.03,...)
	} 		
}
