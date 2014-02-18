HPDsummary <-
function(model,data,genes=NA,relative=FALSE,
log.base=2,summ.plot=TRUE,ptype="z",...) {
	mm=model;base=log.base
	gene.results=list()
	trms=sub("gene:","",attr(terms(mm$Fixed$formula),
	"term.labels")[2:length(attr(terms(mm$Fixed$formula),"term.labels"))])
	sols=colnames(mm$Sol)
	if (is.na(genes[1])) { genes=sub("gene","",sols[grep("gene[^:]*$",sols)])}
	facts =list()
	for (t in trms) {
		if (!grepl(":",t)) { facts =append(facts,list(levels(data[,t])))}
	}
	names(facts)=trms[1:length(facts)]
	nfactors=length(facts)
	if (nfactors>2) { 
		stop("not implemented for more than 2 crossed factors")
	}
	gsols=c();fac1=c();fac2=c();samps=c()
	for (gene in genes) {
		fac1g=c();fac2g=c();sampsg=c()
		for (lev1 in 1:length(facts[[1]])) {
			if (nfactors==2) {
				for (lev2 in 1:length(facts[[2]])) {
					if (lev1==1 & lev2==1) { 
						sol=paste("gene",gene,sep="")
						samp=mm$Sol[,sol]*as.numeric(!relative)
						int0=samp
					} else {
						if (lev2==1) { 
							sol=paste("gene",gene,":",
							names(facts)[1],facts[[1]][lev1],sep="")
							samp=int0+mm$Sol[,sol]
							int2=mm$Sol[,sol]
						} else {
							if (lev1==1) { 
								sol=paste("gene",gene,":",
								names(facts)[2],facts[[2]][lev2],sep="") 
								samp=int0+mm$Sol[,sol]
								int1=mm$Sol[,sol]
							} else {
								sol=paste("gene",gene,":",
								names(facts)[1],facts[[1]][lev1],
								":",names(facts)[2],facts[[2]][lev2],sep="") 
								samp=int0+int1+int2+mm$Sol[,sol]
							}
						}
					}
		#			print(paste(lev1,lev2,sol))
					gsols=append(gsols,gene)
					fac1g=append(fac1g,facts[[1]][lev1])
					fac2g=append(fac2g,facts[[2]][lev2])
					sampsg=cbind(sampsg,samp)
				}
			} else {
				if (lev1==1) { 
						sol=paste("gene",gene,sep="") 
						samp=mm$Sol[,sol]*as.numeric(!relative)
						int0=samp
				} else {
					sol=paste("gene",gene,":",
					names(facts)[1],facts[[1]][lev1],sep="")
					samp=int0+mm$Sol[,sol]
				}
				gsols=append(gsols,gene)
				fac1g=append(fac1g,facts[[1]][lev1])
				sampsg=cbind(sampsg,samp)
			}
		}
		fac1=append(fac1,fac1g)
		samps=cbind(samps,sampsg)
		if (nfactors==2) { fac2=append(fac2,fac2g) }
		sampsg=data.frame(sampsg)
		if (nfactors==2) {
			gres=matrix(nrow=length(fac1g),ncol=length(fac1g),
			dimnames=list("pvalue"=paste(names(facts)[1],fac1g,":",
			names(facts)[2],fac2g,sep=""),
			"difference"=paste(names(facts)[1],fac1g,":",
			names(facts)[2],fac2g,sep="")))
		} else { 
			gres=matrix(nrow=length(fac1g),ncol=length(fac1g),
			dimnames=list("pvalue"=fac1g,"difference"=fac1g))
		}
		for (i in 1:(length(fac1g)-1)) {
			for (j in (i+1):length(fac1g)) {
				diff=sampsg[,j]-sampsg[,i]
				gres[j,i]=mcmc.pval(diff,ptype=ptype)
				gres[i,j]=mean(diff)/log(base)
			}	
		}
		gene.results=append(gene.results,list(gres))
		names(gene.results)[length(gene.results)]=gene		
	}
	
	big.summary=data.frame(cbind("gene"=gsols,"f1"=fac1))
	names(big.summary)[2]=names(facts[1])
	if(nfactors==2) {
		big.summary$f2=fac2
		names(big.summary)[3]=names(facts[2])
	}
	samps=apply(samps,2, function(x){ return(x/log(base)) })
	mns=apply(samps,2,mean)
	sds=apply(samps,2,sd)
	lower=apply(samps,2,function(x) { return(quantile(x,0.05)) })
	upper=apply(samps,2,function(x) { return(quantile(x,0.95)) })
	big.summary=cbind(big.summary,"mean"=mns,
	"sd"=sds,"lower"=lower,"upper"=upper)
	
	if(relative) {
		if (nfactors==2) { 
			remov=which(
				big.summary[,2]==facts[[1]][1] & big.summary[,3]==facts[[2]][1]
			)
			big.summary=big.summary[-remov,] 
		} else { 	
			big.summary=big.summary[big.summary[,2]!=facts[[1]][1],] 
			}
	}

	if(summ.plot) { 
		if (!relative) { 
			if (nfactors==2) {
				summaryPlot(big.summary,xgroup=names(facts[1]),
				facet=names(facts[2]),type="line",...) 
			} else {
				summaryPlot(big.summary,xgroup=names(facts[1]),type="line",...) 
			}
		} else {
			if (nfactors==2) {
				summaryPlot(big.summary,xgroup=names(facts[1]),
				facet=names(facts[2]),type="bar",...) 
			} else {
				summaryPlot(big.summary,xgroup=names(facts[1]),type="bar",...) 
			}
		}
	}

	return(list("summary"=big.summary,"geneWise"=gene.results))
}
