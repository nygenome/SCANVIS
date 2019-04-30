## Phaedra Agius, April 2019
##
## SCANVIS.annotation - script for preparing a gencode object for use in SCANVIS
##
## eg.ftp.url='ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/'

SCANVIS.annotation<-function(ftp.url){

	out.dir=getwd()

	ftp.files <- ls_url(ftp.url)

	gtf.file <- Filter(function(x) grepl("\\d\\.annotation.gtf", x), ftp.files)
	gtf.url <- paste0(ftp.url, gtf.file)
	fout=gsub('gtf.gz','rda',gtf.file)

	download.file(gtf.url, destfile = file.path(out.dir, gtf.file))
	gencode<-rtracklayer::import.gff(file.path(out.dir, gtf.file), version = "2")
	save(gencode, file = fout)


	#################################################
	##read in gencode data and pull components needed
	tmp=unlist(strsplit(unlist(strsplit(fout,'/')),'\\.'))
	v=tmp[grep('v',tmp)]
	print(paste0('****** Loading up gencode data: gencode version ',v,' ******'))

	gencode=get(load(fout)) 
	x=as.matrix(gencode@ranges)
	x[,2]=x[,1]+x[,2]-1
	colnames(x)=c('start','end')
	xg=cbind(as.vector(gencode@seqnames),
			gencode@elementMetadata@listData$transcript_id,
			gencode@elementMetadata@listData$gene_name,
			gencode@elementMetadata@listData$gene_id,
			gencode@elementMetadata@listData$exon_id,
			gencode@elementMetadata@listData$exon_number,
			as.vector(gencode@elementMetadata@listData$type),
			as.vector(gencode@strand),
			as.vector(gencode@elementMetadata@listData$gene_type))
	colnames(xg)=c('chr','transcript','gene_name','gene_id','exon_id',
		'exon_number','type','strand','gene_type')

	gen.tmp=cbind(x,xg)
	q=union(which(gen.tmp[,'type']=='exon'),which(gen.tmp[,'type']=='gene'))
	gen.tmp=gen.tmp[q,]
	rm(gencode)

	#################################################################
	#get intronic regions that do not overlap with any coding regions
	print('***** Collecting intronic coordinates ... *****')
	gen=NULL
	gen$INTRONS=NULL
	for(chr in unique(gen.tmp[,'chr'])){
		print(chr)
		q=which(gen.tmp[,'chr']==chr)
		v=coverage(IRanges(as.numeric(gen.tmp[q,'start']),as.numeric(gen.tmp[q,'end'])))
		#I=cov2coor(v)
		h=which(v@values==0)
		end.pos=unlist(lapply(h,function(x) sum(v@lengths[seq(1,x,1)])))
		tmp=c(1,v@lengths)
		start.pos=unlist(lapply(h,function(x) sum(tmp[seq(1,x,1)])))
		I=cbind(start.pos,end.pos)

		#these should all have zero intersection with genes
		v2=findOverlaps(IRanges(I[,1],I[,2]),
			IRanges(as.numeric(gen.tmp[q,'start']),as.numeric(gen.tmp[q,'end'])))
		if(length(v2)>0) stop()

		#now let's extract all exonic coordinates
		q2=q[which(gen.tmp[q,'type']=='exon')]
		v2=coverage(IRanges(as.numeric(gen.tmp[q2,'start']),as.numeric(gen.tmp[q2,'end'])))
		#I2=cov2coor(v2)
		h=which(v2@values==0)
		end.pos=unlist(lapply(h,function(x) sum(v2@lengths[seq(1,x,1)])))
		tmp=c(1,v2@lengths)
		start.pos=unlist(lapply(h,function(x) sum(tmp[seq(1,x,1)])))
		I2=cbind(start.pos,end.pos)

		tmp=paste(I[,1],I[,2])
		tmp2=paste(I2[,1],I2[,2])
		I2=I2[which(!is.element(tmp2,tmp)),]

		if(length(I)>2) I=cbind(I,'IG')
		if(length(I)==2) I=c(I,'IG')
		if(length(I2)>2) I2=cbind(I2,'nonIG')
		if(length(I2)==2) I2=c(I2,'nonIG')
		if(length(I)*length(I2)>0){
			out.tmp=cbind(chr,I)
			out.tmp=rbind(out.tmp,cbind(chr,I2))
			gen$INTRONS=rbind(gen$INTRONS,out.tmp)
		}
	}
	print('***** DONE: Collecting intronic coordinates *****')
	colnames(gen$INTRONS)=c('chr','start','end','InterGenic')
	d=as.numeric(gen$INTRONS[,'end'])-as.numeric(gen$INTRONS[,'start'])
	gen$INTRONS=gen$INTRONS[which(d>=2),]

	##################################################
	gen$EXONS=gen.tmp[which(gen.tmp[,'type']=='exon'),]
	gen$GENES=gen.tmp[which(gen.tmp[,'type']=='gene'),]
	##################################################

	##################################################
	#merge genes that have overlapping genomic regions
	gen$GENES.merged=NULL
	CHR=unique(gen$GENES[,'chr'])
	for(chr in CHR){
		q=which(gen$GENES[,'chr']==chr)
		IR=IRanges(as.numeric(gen$GENES[q,'start']),as.numeric(gen$GENES[q,'end']))
		M=reduce(IR) #merge intervals
		v=as.matrix(findOverlaps(M,IR))
		tmp=gen$GENES[q[v[,2]],'gene_name']
		tt=table(v[,1])
		if(max(tt)>1){
			for(j in as.numeric(names(tt)[which(tt>1)])){
				q=which(v[,1]==j)
				tmp[q]=rep(paste(unique(tmp[q]),collapse=','),length(q))
			}
		}
		M=as.matrix(M)
		M=cbind(M[,1],M[,1]+M[,2]-1)
		M=cbind(M[v[,1],],tmp)
		M=unique(cbind(chr,M))
		gen$GENES.merged=rbind(gen$GENES.merged,M)
	}
	colnames(gen$GENES.merged)=c('chr','start','end','gene_name')

	return(gen)
}

ls_url <- function(url) {
	print(url)	
	stopifnot(url.exists(url))
	out <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
	readLines(textConnection(out))
}


