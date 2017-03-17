require(MASS)
require(VariantAnnotation)
require(parallel)
require(HardyWeinberg)
require(DESeq2)
require(limma)
require(OutFLANK)
options(stringsAsFactors=F)
#### this function assumes that you have 1. a vcf of biallelic SNPs produced by freebayes
#### 2. the same type of expression matrix that DESeq uses and 3. a factor vector of populations
#### it uses DESeq to test allelic imbalance and voom/lm to test genotype-expression associations within pops
#### then it combines evidence from these tests into an overall eQTL pvalue
#### for convenience, it also can compute Fst and test for differential expression

make012<-function(mat){
	row_names=rownames(mat)
	col_names=colnames(mat)
	allowed=c('0|0','0/0','0/1','1/0','0|1','1|0','1/1','1|1')
	if(any(!mat%in%allowed)){
		stop(paste('unexpected genotype(s):',paste(mat[!mat%in%allowed],collapse=','),'-- are you sure these are biallelic SNPs?'))
	}
	newmat<-mat
	newmat[mat=='0/0'|mat=='0|0']<-0
	newmat[mat=='0/1'|mat=='1/0'|mat=='0|1'|mat=='1|0']<-1
	newmat[mat=='1/1'|mat=='1|1']<-2
	newmat<-apply(newmat,2,as.numeric)
	rownames(newmat)<-row_names
	colnames(newmat)<-col_names
	return(newmat)
}

makeHierfstat<-function(mat,pops,remove.names=F){
	newmat=mat
	newmat[mat==0]<-11
	newmat[mat==1]<-12
	newmat[mat==2]<-22
	newmat<-as.data.frame(cbind(factor(pops),t(newmat)))
	if(remove.names) {rownames(newmat)=NULL;colnames(newmat)=c('pops',1:(ncol(newmat)-1))}
	return(newmat)
}

DEbin<-function(bin,bins,c.ord,p.ord,fitType=c('parametric','local','mean')){
	numBin= length(which(bins==bin))
	currc<-c.ord[bins==bin,]
	currp<-p.ord[bins==bin,][1,]
	if(numBin<20) fitType='mean'
	nafilt=!is.na(currp)
	currc<-currc[,nafilt]
	currp<-currp[nafilt]
	form=as.formula(~currp)
	currdat=data.frame(currp=factor(currp))
	cat(paste('testing sites with',table(currdat$currp)[1],'heterozygotes...\n'))	
	cds<-DESeqDataSetFromMatrix(currc,currdat,form)
	cds<-DESeq(cds,fitType=fitType)
	res<-results(cds)
	return(res)
}

imbalanceTest<-function(AImat,AIdat,fitType=c('parametric','local','mean')){
	cat('reordering into groups with same number of hets...\n')
	counts=AImat
	predictor=do.call(rbind,rep(list(AIdat),nrow(counts)))
	predictor[is.na(counts)]<-NA
	counts[is.na(predictor)]<-NA
	orders<-t(apply(predictor,1,order))
	c.ord<-t(sapply(1:nrow(orders),function(i) counts[i,orders[i,]]))
	rownames(c.ord)<-rownames(counts)
	p.ord<-t(sapply(1:nrow(orders),function(i) predictor[i,orders[i,]]))
	rownames(p.ord)<-rownames(counts)
	notAllNA<-apply(p.ord,1,function(v) any(!is.na(v)))
	c.ord<-c.ord[notAllNA,]
	p.ord<-p.ord[notAllNA,]
	bins<-apply(p.ord,1,function(v) length(which(!is.na(v))))/2
	bintab<-table(bins)[table(bins)>1]
	cat('running DESeq...\n')
	res<-do.call(rbind,sapply(names(bintab),DEbin,bins=bins,c.ord=c.ord,p.ord=p.ord,fitType=fitType))
	res<-as.data.frame(res)[rownames(counts),]
	rownames(res)<-rownames(counts)
	return(res)
}

associationTest<-function(curr,currexpr,currweights,genos,covariates=NULL,withinPop=T){
	gt=genos[curr,]
	components='gt'
	if(withinPop) components=c('pops','gt')
	if(!is.null(covariates)) components=c('covariates',components)
	form<-as.formula(paste('ex~',paste(components,collapse='+')))
	ex= currexpr[curr,]
	if(is.null(currweights)) {
		lm.out<-summary(lm(form))$coefficients
	} else{
		wt= currweights[curr,]
		lm.out<-summary(lm(form,weights=wt))$coefficients
	}
	if(!'gt'%in%rownames(lm.out)) {
		warning(paste('could not fit',curr,' -- skipping and returning NA...'))
		return(rep(NA,3))
	}
	lm.res<-lm.out['gt',c('t value','Pr(>|t|)')]
	fc=2^(lm.out['gt','Estimate']*2)
	return(c(lm.res,fc))
}

propExplain<-function(curr,currexpr,genos,pops){
	ex=currexpr[curr,]
	gt=genos[curr,]
	aov1<-summary(aov(ex~pops))[[1]]	
	aov2<-summary(aov(ex~gt+pops))[[1]]	
	popDiffExplained=max(1-aov2['pops','Sum Sq']/aov1['pops','Sum Sq'],0)
	return(popDiffExplained)
}

vcf2eqtl<-function(vcf,
expr,
pops=NULL,
minHet=3,
mc.cores=1,
alpha=0.05,
calculateFst=T,
testDE=F,
all3=T,
hweFilter=T,
hweAlpha=0.05,
covariates=NULL,
propExplained=T,
withinPop=T,
transcripts=NULL,
keepSamples=NULL){
	
	if(is.null(pops)){
		calculateFst=F
		testDE=F
		propExplained=F
		withinPop=F
	}
	#extract SNP genotype and SNP depth data from vcf
	globalFst=NULL
	cat('reading vcf...\n')
	currvcf<-suppressWarnings(readVcf(vcf,genome='curr'))
	genoInfo<-geno(currvcf)
	genos<-make012(genoInfo$GT)
	AOs<-apply(genoInfo$AO,2,as.numeric)
	ROs<-apply(genoInfo$RO,2,as.numeric)
	AOs[genos!=1]<-NA
	ROs[genos!=1]<-NA
	AImat<-cbind(ROs,AOs)
	rownames(AImat)<-rownames(genos)
	AIdat<-factor(c(rep('ref',ncol(genos)),rep('alt',ncol(genos))),levels=c('ref','alt'))

	#collect SNP metadata	
	CHROM=as.vector(seqnames(rowRanges(currvcf)))
	POS=as.numeric(as.character(start(ranges(rowRanges(currvcf)))))
	REF=as.character(mcols(rowRanges(currvcf))[,'REF'])
	ALT=as.character(unlist(mcols(rowRanges(currvcf))[,'ALT']))
	AF=unlist(info(currvcf)$AF)
	snpInfo<-cbind(CHROM,POS,REF,ALT,AF)
	rownames(snpInfo)<-rownames(genos)

	#organize and normalize expression data
	if(is.null(transcripts)) transcripts=CHROM
	expr<-expr[rownames(expr)%in%transcripts,]
	voomExpr<-voom(expr)
	rownames(voomExpr$weights)<-rownames(expr)
	currexpr<-as.matrix(voomExpr$E[transcripts,])
	currweights<-as.matrix(voomExpr$weights[transcripts,])
	rownames(currexpr)<-rownames(genos)
	rownames(currweights)<-rownames(genos)
	
	#subset samples if desired
	if(!is.null(keepSamples)){
		cat('subsetting to only keep specified samples and assuming pops correspond to post-filtered samples...\n')
		keep=colnames(genos)%in%keepSamples
		expr<-expr[,keep]
		currexpr<-currexpr[,keep]
		currweights<-currweights[,keep]
		AImat<-AImat[,rep(keep,2)]
		AIdat<-AIdat[keep]
		genos<-genos[,keep]
	}

	#if desired, filter genotypes	
	if(all3){
		cat('filtering for sites with all three genotypes...\n')
		num_gts<-apply(genos,1,function(v) length(table(v)))
		genos<-genos[num_gts==3,]
	}
	if(hweFilter){
		cat('filtering out sites out of HWE in at least one population...\n')
		hwepops=pops
		hwe=rep(1,nrow(genos))
		if(is.null(pops)) hwepops=rep(1,ncol(genos))
		for(pop in unique(hwepops)){
			hwe<-pmin(hwe,apply(genos[,which(pops==pop)],1,function(v) HWExact(table(factor(v,levels=c(0,1,2))),verbose=F)$pval))
		}
		genos<-genos[hwe>hweAlpha,]
	}
	if(!is.null(minHet)){
		cat('filtering out sites with fewer than',minHet,'heterozygotes...\n')
		numHet<-apply(genos,1,function(v) length(which(v==1)))
		genos<-genos[numHet>=minHet,]
	}

	cat('filtering out sites without allele observations\n')
	AImat<-AImat[rownames(genos),]
	imbalanceInfo<-apply(AImat,1,function(v) any(!is.na(v)))
	genos<-genos[rownames(AImat)[imbalanceInfo],]
	if(nrow(genos)==0){
		stop('No sites left after filtering, check to make sure you have a freebayes VCF of biallelic SNPs with allele observations in it')
	}	
	#subset other data sets after filtering genos
	AImat<-AImat[rownames(genos),]
	currexpr<-currexpr[rownames(genos),]
	snpInfo<-snpInfo[rownames(genos),]
	cat(paste(nrow(genos),'sites left after filtering, testing for eQTL status...\n'))

	cat('allelic imbalance test...\n')
	imb.out<-imbalanceTest(AImat=AImat,AIdat=AIdat)
	
	#single-threaded test
	if(mc.cores==1){	
		cat('association test...\n')
		assoc.out<-t(sapply(rownames(genos),associationTest,
		currexpr=currexpr,currweights=currweights,genos=genos,withinPop=withinPop))
	} else{
	#multithreaded test
		cat('mulithreaded association test...\n')
		assoc.out<-do.call(rbind,mclapply(rownames(genos),associationTest,
			currexpr=currexpr,currweights=currweights,genos=genos,
			withinPop=withinPop,mc.cores=mc.cores))
	}

	#collect results and calculate p values using Stouffer's method
	res<-cbind(imb.out[,c('stat','pvalue')],2^imb.out$log2FoldChange,assoc.out)
	colnames(res)<-c('AIz','AIp','AIfc','ASSOCz','ASSOCp','ASSOCfc')	
	rownames(res)<-rownames(genos)
	res<-as.data.frame(res)
	res$z<-(res[,'AIz']+res[,'ASSOCz'])/(2**.5)
	# res$z[is.na(res$AIz)&!is.na(res$ASSOCz)]<-res$ASSOCz[is.na(res$AIz)&!is.na(res$ASSOCz)]
	# res$z[!is.na(res$AIz)&is.na(res$ASSOCz)]<-res$ASSOCz[!is.na(res$AIz)&is.na(res$ASSOCz)]
	res$p<-sapply(res$z,function(val) min(pnorm(val,lower.tail=T),pnorm(val,lower.tail=F))*2)
	res$padj<-p.adjust(res$p,method='BH')
	res$AIpadj<-p.adjust(res$AIp,method='BH')
	res$ASSOCpadj<-p.adjust(res$ASSOCp,method='BH')
	res$eQTL<-res$padj<alpha
	res<-cbind(snpInfo,res)
	res$REF<-as.character(res$REF)
	res$ALT<-as.character(res$ALT)
	res$POS<-as.numeric(as.character(res$POS))
	
	#calculate Fst and call outliers using OutFLANK
	if(calculateFst){
		cat('calculating Fst...\n')
		wc.out<-MakeDiploidFSTMat(t(genos),rownames(genos),pops)
		res$Fst=wc.out$FST
		res$FstNum<-wc.out$T1
		res$FstDen<-wc.out$T2
		fl.out<-OutFLANK(wc.out)$results
		globalFst=fl.out$FSTbar
		res$FstOutlier=fl.out$OutlierFlag
		res$FstOutlierP=fl.out$pvaluesRightTail
	}
	
	#test for differential expression using DESeq2
	DEres =NULL
	if(testDE){
		cat('testing for differential expression with DESeq2...\n')
		currdat<-data.frame(pops=pops)
		cds<-DESeqDataSetFromMatrix(expr,currdat,~pops)
		cds<-DESeq(cds,test='LRT',reduced=~1)
		DEres <-results(cds)
		res$DE<-DEres[res$CHROM,'padj']<alpha
		res$DEp<-DEres[res$CHROM,'pvalue']
		res$DEpadj<-DEres[res$CHROM,'padj']
	}
	
	#calculate reduction in population differentiation after accounting for eQTL
	if(propExplained){
		cat('calculating proportion of population differences explained by eQTLs...\n')
		if(mc.cores==1){
			propE<-sapply(rownames(genos)[which(res$eQTL)],propExplain,currexpr=currexpr,genos=genos,pops=pops)
		} else{
		cat('\tmulithreading...\n')
			propE<-unlist(mclapply(rownames(genos)[which(res$eQTL)],propExplain,
				currexpr=currexpr,genos=genos,pops=pops,mc.cores=mc.cores))
		}
		res$popDiffExplained<-NA
		res$popDiffExplained[which(res$eQTL)]=propE
	}
	
	return(list(res=res,snpContigExpr=currexpr,genos=genos,AImat=AImat,AIdat=AIdat,globalFst=globalFst,pops=pops,DEres= DEres))
}

flattenOutput<-function(output){
	ex<-output$snpContigExpr
	colnames(ex)<-paste(colnames(ex),'_expr',sep='')
	gt<-output$genos
	colnames(gt)<-paste(colnames(gt),'_geno',sep='')
	return(cbind(output$res,ex,gt))
}

plotAssociation<-function(output,curr,titles=T,plotLegend=T,legendPos=NULL,snpLab=NULL,usecol=NULL,geneName=NULL,lty=2,lwd=2,pops=NULL){
	if(is.null(pops)){pops=output$pops}
	if(is.null(pops)) {
		plotLegend=F
		cols=rep('black',ncol(output$genos))
	}
	ex=output$snpContigExpr[curr,]
	gt=output$genos[curr,]
	currsnp<-output$res[curr,]
	if(!is.null(pops)) cols=rainbow(length(levels((pops))),v=0.8)[as.numeric((pops))]
	if(!is.null(usecol)) {cols=usecol; plotLegend=F}
	xlab=paste(currsnp$CHROM,'pos.', currsnp['POS'])
	if(!is.null(snpLab)) xlab=snpLab
	plot(ex~jitter(gt,0.5),col=cols,pch=19,axes=F, ylab='Log2 Normalized Counts', xlab=xlab,xlim=c(-.1,2.1))
	abline(lm(ex~gt),lwd=lwd,lty=lty)
	axis(2)
	legpos='topleft'
	if(currsnp$z<0) legpos='topright'
	if(!is.null(legendPos)) legpos=legendPos
	if(plotLegend) legend(legpos,fill=rainbow(length(levels(pops)),v=0.8),legend=levels(pops),bty='n')
	axis(1,at=c(0,1,2),labels=c(paste(currsnp['REF'],currsnp['REF'],sep=''),
		paste(currsnp['REF'],currsnp['ALT'],sep=''),
		paste(currsnp['ALT'],currsnp['ALT'],sep='')))
	if(titles) title('eQTL Association')
	box()
	if(!is.null(geneName)) mtext(geneName,1,2)
}

plotImbalance<-function(output,curr,titles=T,usecol=NULL,geneName=NULL,pops=NULL){
	if(is.null(pops)){pops=output$pops}
	if(is.null(pops)) {
		cols=rep('black',ncol(output$genos))
	} else{
		cols=rainbow(length(levels(factor(pops))),v=0.8)[as.numeric(factor(pops))]
	}
	currsnp<-output$res[curr,]
	if(!is.null(usecol)) cols= usecol; legend=F
	AIs<-cbind(output$AImat[curr,output$AIdat=='ref'],output$AImat[curr,output$AIdat=='alt'])
	c1=cols[!is.na(AIs[,1])]
	if(length(which(!is.na(AIs[,1])))==0){stop('Cannot plot imbalance without ref and alt read depths, check AO and RO VCF fields')}
	AIs<-na.omit(AIs)
	plot(1,type='n',xlim=c(0.8,2.2),ylim=c(min(AIs),max(AIs)),axes=F,xlab='Allele',ylab='Counts')
	axis(1,at=1:2,labels=c(currsnp['REF'],currsnp['ALT']))
	axis(2)
	sapply(1:nrow(AIs),function(ind) lines(AIs[ind,],type='b',lwd=2,col=c1[ind],pch=19))
	if(titles) title('Allelic Imbalace')
	box()
}

plotDistribution<-function(output,curr,titles=T,usecol=NULL,plotLegend=T,legendPos='topright',xlab='Population',pops=NULL){
	if(is.null(pops)) {
		if(!is.null(output$pops)){
			pops=output$pops
		} else{
		stop('Need populations to display allele distributions between populations...')
		}
	}
	cols=rainbow(length(levels(pops)),v=0.8)[1:length(levels(pops))]
	if(!is.null(usecol)) cols= usecol; legend=F
	agg<-aggregate(output$genos[curr,]~pops,FUN=mean)
	xlim=c(1,length(levels(pops)))
	xlim=mean(xlim)+(xlim-mean(xlim))*1.25
	ylim=c(min(agg[,2]),max(agg[,2]))/2
	ylim=mean(ylim)+(ylim-mean(ylim))*1.25
	plot(agg[,2]/2,col=cols,pch=19,cex=2,xlab=xlab,ylab='Alternate allele freq.',axes=F,xlim=xlim,ylim=ylim)
	axis(1,at=1:length(levels(pops)),levels(pops))
	axis(2)
	lines(agg[,2]/2,lty=2,lwd=2)
	box()
	if(titles) title('Allele distribution')
	if(plotLegend) {legend(legendPos,fill=cols,legend=levels(pops),bty='n')}
}

plotCline<-function(cline,output,curr,titles=T,usecol=NULL,plotLegend=T,legendPos='topright',xlab='Latitude'){
	clinefac<-factor(cline)
	cols=rainbow(length(levels(clinefac)),v=0.8)[1:length(levels(clinefac))]
	if(!is.null(usecol)) cols= usecol; legend=F
	agg<-aggregate(output$genos[curr,]~cline,FUN=mean)
	xlim=c(min(agg[,1]),max(agg[,1]))
	xlim=mean(xlim)+(xlim-mean(xlim))*1.25
	ylim=c(min(agg[,2]),max(agg[,2]))/2
	ylim=mean(ylim)+(ylim-mean(ylim))*1.25
	plot(agg[,1],agg[,2]/2,col=cols,pch=19,cex=2,xlab=xlab,ylab='Alternate allele freq.',axes=F,xlim=xlim,ylim=ylim)
	axis(1)
	axis(2)
	lines(agg[,1],agg[,2]/2,lty=2,lwd=2)
	box()
	if(titles) title('Allele distribution')
}

plotQTL<-function(output,curr,
	titles=T,
	plotLegend=T,
	ourpar=T,
	legendPos='topright',
	snpLab=NULL,
	col=NULL,
	geneName=NULL,
	lty=2,
	lwd=2,
	cline=NULL,
	pops=NULL,
	plotDistribution=T){
	if(ourpar) {
		par(tck=-0.01,mgp=c(1.2,0.2,0),mar=c(3,3,2,1),bty='l')
		if(plotDistribution) {
			par(mfrow=c(1,3))	
		} else{
			par(mfrow=c(1,2))
		}
	}
	plotAssociation(output,curr,titles,plotLegend=F,legendPos,snpLab,col,geneName,lty,lwd,pops=pops)
	plotImbalance(output,curr,titles,col,geneName,pops=pops)
	if(plotDistribution){
		if(is.null(cline)){
			plotDistribution(output,curr,titles,col,plotLegend,legendPos,pops=pops)
		} else{
			plotCline(cline,output,curr,titles,col,plotLegend,legendPos)
		}
	}
}