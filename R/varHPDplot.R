varHPDplot <-
function(model,factors="units",factors2=NULL, ylimits=NULL,hpdtype="w",jitter=0,plot=T,testlimit=0.05,grid=T,zero=T,...){
	res1=res2=c(rep(0,length(model$VCV[,1])))
	names1=names2=c()
	allnames=names(posterior.mode(model$VCV))
	solnames=names(posterior.mode(model$Sol))
	genes=solnames[-grep(":",solnames)]
	genes=genes[grep("gene",genes)]
	genes=sub("gene","",genes)
	ngenes=length(genes)
	first=1
	for (f in factors){		
		pattern=paste('^\\w+\\.',f,'$',sep="")
		f1=grep(pattern,allnames)
		res1=res1+model$VCV[,f1]
		if (first==1) {
			names1=allnames[f1]
			first=0
		} else {
			names1=paste(names1," + ",allnames[f1],sep="")
		}
	}
	
	if (is.null(factors2)==FALSE) {
		first=1
		for (f in factors2){		
			pattern=paste('^\\w+\\.',f,'$',sep="")
			f1=grep(pattern,allnames)
			res2=res2+model$VCV[,f1]
			if (first==1) {
				names2=allnames[f1]
				first=0
			} else {
				names2=paste(names2," + ",allnames[f1],sep="")
			}
		}
		stat=res2/(res1+res2)
		means=apply(stat,2,mean)
		hpds=HPDinterval(stat)
		ps=mcmc.pval(stat,testlim=testlimit,sided=1)
		names.final=paste(rep("contribution of ",length(names2)),names2,sep="")		
		ylabel="proportion of variance"
	} else {
		first=1
		res1=res2=c(rep(0,length(model$VCV[,1])))
		for (f in factors){		
			pattern=paste('^\\w+\\.',f,'$',sep="")
			f1=grep(pattern,allnames)
			res1=res1+model$VCV[,f1]
			if (first==1) {
				names1=allnames[f1]
				first=0
			} else {
				names1=paste(names1," + ",allnames[f1],sep="")
			}
		}
		stat=log(sqrt(exp(res1)),2)
		means=apply(stat,2,mean)
		hpds=HPDinterval(stat)
		ps=mcmc.pval(stat,testlim=testlimit,sided=1)
		names.final=names1
		ylabel="deviation, log2(fold change)"
	}
	if (plot==F) {
		table.final=data.frame(cbind("mean"=means,hpds,"pval.mcmc"=ps))
		row.names(table.final)=names.final
		return (table.final)
		stop
	}
	if(is.null(ylimits) && is.null(factors2)) {
		ylimits=c(min(hpds[,1])-0.1,max(hpds[,2])+0.1)
	} else if (is.null(ylimits)) ylimits=c(0,1)
	plot(c(1:ngenes)+jitter,means,xlim=c(0.5,ngenes+0.5),xaxt="n",ylim=ylimits,xlab="",ylab=ylabel,mgp=c(2.3,1,0),...)
	if (hpdtype=="l") {
		lines(hpds[,1],lty=2,...)
		lines(hpds[,2],lty=2,...)
	}	else {
		arrows(c(1:ngenes)+jitter,means,c(1:ngenes)+jitter,hpds[,1],angle=90,code=2,length=0.03,...)
		arrows(c(1:ngenes)+jitter,means,c(1:ngenes)+jitter,hpds[,2],angle=90,code=2,length=0.03,...)
		if (grid==TRUE){
			verts=seq(0.5,ngenes+0.5,1)
			abline(v=verts,lty=3,col="grey60")
		} 		
	} 		
	if (zero==TRUE) abline(h=0,lty=3,col="grey60")
	axis(side=1,labels=genes,at=c(1:ngenes),las=2)
}
