mcmc.qpcr.lognormal <-
function(fixed,random=NULL,data,controls=NULL,include=NULL,m.fix=1.2,v.fix=NULL,genebysample=T,vprior="flat",...) {
	ngenes=length(levels(data[,"gene"]))
	if(vprior=="flat") { 
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
				print("vprior is not recognized, should be flat, iw, or iw01")
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
		ff=paste("count~0+gene+gene:",ff,sep="")
	}	

# assembling the random effects formula:
	rr="~sample"
	if(reps & genebysample)	random=append(random,"sample")
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
	Rr=list(V=diag(ngenes), nu=0)
	
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
mc=MCMCglmm(formula(ff),random=formula(rr),rcov=~idh(gene):units,data=data,prior=prior,...)
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
		
		gm1=MCMCglmm(formula(ff),data=data,verbose=F,nitt=100,thin=10,burnin=2)
		fixs=names(posterior.mode(gm1$Sol))
		fl=length(fixs)
		m.fix=m.fix+0.000001
		Mu=rep(0,length(fixs)) # zero mean priors for fixed effects
		va.m=diag(length(fixs))*1e+8 # very large prior variances for fixed effects on non-normalizers;
		
		# fishing out the numbers of control-related effects:
		controls=controls[length(controls):1]
		nn=c()
		if (include>0) {
			for (n in controls[1:include]){
				nn=append(nn,grep(n,fixs))
				nn=nn[nn>ngenes]
			}
			va.m[nn,nn]=diag(length(nn))*log(m.fix^2) # specifying prior variance for main effects of normalizers;
		}
				
		if (is.null(v.fix)==FALSE & include>0){
			v.fix=v.fix+0.000001
			va.r[(var.estimate+1):ngenes,(var.estimate+1):ngenes]=diag(length((var.estimate+1):ngenes))*log(v.fix^2)

			if (reps==1 & genebysample) { # fixing all but residual variance
				for (ri in 1:length(random)) {
					Gstruct[[paste("G",1+ri,sep="")]]=list(V=va.r, nu=0,fix=var.estimate+1)
				}
			} else if (reps==0 & genebysample) { # fixing residual variance
					Rr=list(
						V=va.r, nu=0,fix=var.estimate+1
						)
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
			prior=list(  
				B=list(mu=Mu,V=va.m),
				R=list(
					V=diag(ngenes), nu=0), 
				G=Gstruct
			)
	
		}
		print(list("PRIOR"=prior,"FIXED"=ff,"RANDOM"=rr))
mc=MCMCglmm(formula(ff),random=formula(rr),rcov=~idh(gene):units,data=data,prior=prior,...)
		return(mc)
	}
}
