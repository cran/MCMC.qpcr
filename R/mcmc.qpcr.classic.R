mcmc.qpcr.classic <-
function(fixed,random=NULL,data,controls,genebysample=T,center=T,...) {
	ngenes=length(levels(data[,"gene"]))

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
	rr="~"
	if(reps & genebysample)	random=append(random,"sample")
	for (r in random){
			if (r==random[1]) rr=paste(rr,"idh(gene):",r,sep="") else rr=paste(rr,"+idh(gene):",r,sep="") 
	}
# moving controls to the end of the factor list (for variance fixing later)
	controls=controls[length(controls):1]
	relev=c()
	for (g in levels(data[,"gene"])) {
		if (g %in% controls) {
			next
		}
		relev=append(relev,g)
	}
	relev=append(relev,controls)
	data[,"gene"]=factor(data[,"gene"],levels=relev)

	ddn=normalize.qpcr(data,controls,center)
	
mc=MCMCglmm(formula(ff),random=formula(rr),rcov=~idh(gene):units,data=ddn,...)
	return(mc)
}
