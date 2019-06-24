## Phaedra Agius, June 2019
##
## SCANVIS.annotation - script for preparing a gencode object for use in SCANVIS
##
## Eg.
## ftp.url='ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/'

SCANVIS.annotation<-function(ftp.url){

    ftp.files <- ls_url(ftp.url)

    gtf.file <- Filter(function(x) grepl("\\d\\.annotation.gtf", x), ftp.files)
    gtf.url <- paste0(ftp.url, gtf.file)
    fout=gsub('gtf.gz','rda',gtf.file)

    download.file(gtf.url, destfile = file.path(out.dir, gtf.file))
    gencode<-rtracklayer::import.gff(file.path(out.dir, gtf.file),version="2")

    #################################################
    ##read in gencode data and pull components needed
    tmp=unlist(strsplit(unlist(strsplit(fout,'/')),'\\.'))
    v=tmp[grep('v',tmp)]
    message(paste0('**** Loading up gencode data: gencode version ',v,' ****'))

    x=as.matrix(ranges(gencode))
    x[,2]=x[,1]+x[,2]-1
    colnames(x)=c('start','end')
    tmp=c('transcript_id','gene_name','gene_id','exon_id',
          'exon_number','type','gene_type')
        tmp=elementMetadata(gencode)[,tmp]
        xg=matrix('',nrow(tmp),ncol(tmp))
        for(i in seq(1,ncol(xg),1))
        xg[,i]=as.vector(tmp[,i])
        xg=cbind(as.vector(seqnames(gencode)),xg,as.vector(strand(gencode)))
        colnames(xg)=c('chr','transcript','gene_name','gene_id','exon_id',
               'exon_number','type','gene_type','strand')
    gen.tmp=cbind(x,xg)
    q=union(which(gen.tmp[,'type']=='exon'),which(gen.tmp[,'type']=='gene'))
    gen.tmp=gen.tmp[q,]
    rm(gencode)

    qE=which(gen.tmp[,'type']=='exon')
    qG=which(gen.tmp[,'type']=='gene')
    gen=list(EXONS=gen.tmp[qE,],GENES=gen.tmp[qG,],
        GENES.merged=gen.tmp[qG,c('chr','start','end','gene_name')],
        INTRONS=matrix('',length(qE),4))

    ##################################################
    #merge genes that have overlapping genomic regions
    CHR=unique(gen[['GENES']][,'chr'])
    for(chr in CHR){
        q=which(gen[['GENES']][,'chr']==chr)
        IR=IRanges(as.numeric(gen[['GENES']][q,'start']),
            as.numeric(gen[['GENES']][q,'end']))
        M=reduce(IR) #merge intervals
        v=as.matrix(findOverlaps(M,IR))
        tmp=gen[['GENES']][q[v[,2]],'gene_name']
        tt=table(v[,1])
        if(max(tt)>1){
            for(j in as.numeric(names(tt)[which(tt>1)])){
                h=which(v[,1]==j)
                tmp[h]=rep(paste(unique(tmp[h]),collapse=','),length(h))
            }
        }
        M=as.matrix(M)
        M=cbind(M[,1],M[,1]+M[,2]-1)
        M=cbind(M[v[,1],],tmp)
        M=unique(cbind(chr,M))

        q1=q[1:nrow(M)]
        q0=setdiff(q,q1)
        gen[['GENES.merged']][q1,]=M
        gen[['GENES.merged']][q0,]=''
    }
    h=which(apply(gen[['GENES.merged']]=='',1,sum)!=4)
    gen[['GENES.merged']]=gen[['GENES.merged']][h,]

    #################################################################
    #get intronic regions that do not overlap with any coding regions
    message('***** Collecting intronic coordinates ... *****')
    for(chr in unique(gen.tmp[,'chr'])){
        message(chr)
        q=which(gen.tmp[,'chr']==chr)
        v=coverage(IRanges(as.numeric(gen.tmp[q,'start']),
            as.numeric(gen.tmp[q,'end'])))
        #I=cov2coor(v)
        h=which(v@values==0)
        end.pos=unlist(lapply(h,function(x) sum(v@lengths[seq(1,x,1)])))
        tmp=c(1,v@lengths)
        start.pos=unlist(lapply(h,function(x) sum(tmp[seq(1,x,1)])))
        I=cbind(start.pos,end.pos)

        #these should all have zero intersection with genes
        v2=findOverlaps(IRanges(I[,1],I[,2]),
            IRanges(as.numeric(gen.tmp[q,'start']),
                as.numeric(gen.tmp[q,'end'])))
        if(length(v2)>0) stop()

        #now let's extract all exonic coordinates
        q2=q[which(gen.tmp[q,'type']=='exon')]
        v2=coverage(IRanges(as.numeric(gen.tmp[q2,'start']),
            as.numeric(gen.tmp[q2,'end'])))
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
        if(length(I)>0 & length(I2)>0){
            out.tmp=cbind(chr,I)
            out.tmp=rbind(out.tmp,cbind(chr,I2))
            Q=which(gen[['INTRONS']][,1]=='')[1:nrow(out.tmp)]
            gen[['INTRONS']][Q,]=out.tmp
        }
    }
    gen[['INTRONS']]=gen[['INTRONS']][which(gen[['INTRONS']][,1]!=''),]
    colnames(gen[['INTRONS']])=c('chr','start','end','InterGenic')
    d=as.numeric(gen[['INTRONS']][,'end'])-
      as.numeric(gen[['INTRONS']][,'start'])
    gen[['INTRONS']]=gen[['INTRONS']][which(d>=2),]

    message('***** DONE: Collecting intronic coordinates *****')


    return(gen)
}

ls_url <- function(url) {
    message(url)    
    stopifnot(url.exists(url))
    out <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
    readLines(textConnection(out))
}

