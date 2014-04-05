mcmc.qpcr.lognormal <-
function(fixed,random=NULL,data,controls=NULL,normalize=FALSE,include=NULL,m.fix=1.2,v.fix=NULL,geneSpecRes=TRUE,vprior="uninf",...) {
	if (normalize) {
		if (is.null(controls)) { 
			stop("Error: normalize=TRUE is specified wihtout specifying controls")
		}
		ccc=c();saa=c();ccs=c()
		for (sa in unique(data$sample)){
			for (cc in controls) {
				ccs=append(ccs,round(mean(data[data$gene==cc & data$sample==sa,"count"],na.rm=TRUE),0))
				saa=append(saa,sa)
				ccc=append(ccc,cc)
			}
		}
		CTR=data.frame(cbind("count"=ccs,"gene"=ccc,"sample"=saa))
		CTR$count=as.integer(as.character(CTR$count))
		CTR=data.frame(unstack(CTR,count~gene))
		CTR$sample=unique(data$sample)
		if (length(controls)>1) {
			CTR$NORM=round(apply(CTR[,controls],1,mean),3)
		}
		else { CTR$NORM=CTR[,controls[1]] }		
		for (sa in CTR$sample){
			dss=data[data$sample==sa,]
			dss=dss[1,]
			dss$count=CTR[CTR$sample==sa,"NORM"]
			dss$gene="NORM"
			data=rbind(data,dss)
		}
		data$gene=relevel(data$gene,ref="NORM")
	}

	ngenes=length(levels(data[,"gene"]))
	if(vprior=="uninf") { 
		vstr.g1=list(V=1, nu=0)
		vstr.gs=list(V=diag(ngenes), nu=0)
	} else {
		if (vprior=="iw") {
			vstr.g1=list(V=1, nu=0.002, alpha.mu=0, alpha.V=1000)
			vstr.gs=list(V=diag(ngenes), nu=ngenes)
		} else {
			if (vprior=="iw01") {
				vstr.g1=list(V=1, nu=0.002, alpha.mu=0, alpha.V=1000)
				vstr.gs=list(V=diag(ngenes)*0.1, nu=ngenes)
			} else { 
				print("vprior is not recognized, should be uninf, iw, or iw01")
				stop
				}
		}
	}

#checking if there are technical reps
	g1=levels(data[,"gene"])[1]
	ss=data[data[,"gene"]==g1,]
	if (length(levels(ss$sample))<length(ss[,1])) reps=1 else reps=0
	
# assembling the fixed effects formula:
	if (is.null(fixed)) {
		ff="count~0+gene"
	} else {
		ff=gsub('\\s?\\+\\s?',"+gene:",x=fixed,perl=TRUE)
		if(normalize) {
			ff=paste("count~0+gene+",fixed,"+gene:",ff,sep="")
		}
		else {
			ff=paste("count~0+gene+gene:",ff,sep="")
		}
	}	

# assembling the random effects formula:
	rr="~sample"
	for (r in random){
			rr=paste(rr,"+idh(gene):",r,sep="")
	}

	if (is.null(controls)==F){
		if (is.null(include)) {
			include=length(controls)
		} else if (include>length(controls)) {
			stop("'include' cannot exceed the number of controls")
		}
	}	
	Gstruct=list(G1=vstr.g1)
	va.r=diag(ngenes)
	va.r.cov=diag(ngenes)*0
	if (geneSpecRes) {
		Rr=list(V=diag(ngenes), nu=ngenes-1+0.02)
	} else { Rr=list(V=1, nu=0) }
	
	if (is.null(controls)) {
		if(length(random)>0){
			for (ri in 1:length(random)) {
				Gstruct[[paste("G",1+ri,sep="")]]=vstr.gs
			}
		}
		prior=list(  # inverse gamma, no normalizers 
				R=Rr, 
				G=Gstruct
				)
				print(list("PRIOR"=prior,"FIXED"=ff,"RANDOM"=rr))
				if (geneSpecRes) {
					mc=MCMCglmm(
						formula(ff),random=formula(rr),
						rcov=~idh(gene):units,
						data=data,prior=prior,...)
				} else {
					mc=MCMCglmm(
						formula(ff),random=formula(rr),
						data=data,prior=prior,...)
				}
				return(mc)
	} else {
		
# moving controls to the end of the factor list (for variance fixing later)
		controls=controls[length(controls):1]
		relev=c()
		for (g in levels(data[,"gene"])) {
			if (g %in% controls) {
				next
			}
			relev=append(relev,g)
		}
		var.estimate=length(relev)+length(controls)-include
		relev=append(relev,controls)
		data[,"gene"]=factor(data[,"gene"],levels=relev)
		if (normalize) { data$gene=relevel(data$gene,ref="NORM") }
		
		gm1=MCMCglmm(formula(ff),data=data,verbose=F,nitt=100,thin=10,burnin=2)
		fixs=names(posterior.mode(gm1$Sol))
		fl=length(fixs)
		m.fix=m.fix+0.001
		Mu=rep(0,length(fixs)) # zero mean priors for fixed effects
		va.m=diag(length(fixs))*1e+8 # very large prior variances for fixed effects on non-normalizers;
		
		# fishing out the numbers of control-related effects:
		controls=controls[length(controls):1]
		nn=c()
		if (include>0 & !normalize) {
			for (n in controls[1:include]){
				nn=append(nn,grep(n,fixs))
				nn=nn[nn>ngenes]
			}
			va.m[nn,nn]=diag(length(nn))*log(m.fix^2) # specifying prior variance for main effects of normalizers;
		}
				
		if (is.null(v.fix)==FALSE & include>0 & !normalize){
			v.fix=v.fix+0.001
			va.r[(var.estimate+1):ngenes,(var.estimate+1):ngenes]=diag(length((var.estimate+1):ngenes))*log(v.fix^2)

		if (geneSpecRes) { # fixing residual variance
					Rr=list(V=va.r, nu=0,fix=var.estimate+1)
			} else { # not fixing variance components
				if(length(random)>0){
					for (ri in 1:length(random)) {
						Gstruct[[paste("G",1+ri,sep="")]]=vstr.gs
					}
				}
			}
			prior=list( 
				B=list(mu=Mu,V=va.m),
				R=Rr, 
				G=Gstruct
			)
		} else {
			if(length(random)>0){
				for (ri in 1:length(random)) { # not fixing variance components
					Gstruct[[paste("G",1+ri,sep="")]]=vstr.gs
				}
			}
			if (geneSpecRes) { 
				Rr=list(V=diag(ngenes), nu=ngenes-1+0.02)
			} else {
				Rr=list(V=1, nu=0)
			}	 
			prior=list(  
				B=list(mu=Mu,V=va.m),
				R=Rr, 
				G=Gstruct
			)	
		}
		print(list("PRIOR"=prior,"FIXED"=ff,"RANDOM"=rr))
		if (geneSpecRes) {
			mc=MCMCglmm(
				formula(ff),random=formula(rr),
				rcov=~idh(gene):units,
				data=data,prior=prior,...)
		} else {
			mc=MCMCglmm(
				formula(ff),random=formula(rr),
				data=data,prior=prior,...)
		}
		return(mc)
	}
}
