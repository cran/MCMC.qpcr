padj.hpdsummary=function(hpdsumm,controls=NULL,method="BH") {
	pvs=c()
	pvn=c()
	for (g in names(hpdsumm$geneWise)) {
		if (g %in% controls) {next}
		rr= hpdsumm$geneWise[[g]]
		pvs=append(pvs, rr[lower.tri(rr)])
		pvn=append(pvn,rep(g,length(rr[lower.tri(rr)])))
	}
	pvs=p.adjust(pvs,method=method)
	for (g in names(hpdsumm$geneWise)) {
		if (g %in% controls) {next}
		hpdsumm$geneWise[[g]][lower.tri(hpdsumm$geneWise[[g]])]=pvs[pvn==g]
	}
	return(hpdsumm)
}
