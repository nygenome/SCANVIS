##Phaedra Agius, April 2019
##
##INPUT: scn=urls to scn files output by SCANVIS scan or linkvar functions
##            OR scn can also be a list of matrices
##        method='mean' or 'median'
##        roi=region of interest chr,start,end
##        gen=gencode data
##
##OUTPUT: object with components sj, sj.samples and muts (muts available only if 
##            scn are outputs of SCANVISlinkvar)
##           out$sj=averaged RRS and uniq.reads across samples
##           out$sj.samples=comma sep samples, RRSs and uniq.reads for each sj
##           out$muts=all mutations mapped to sjs in out$sj, with comma sep
##           samples and number of samples containing the mutation
SCANVISmerge<-function(scn,method='mean',roi=NULL,gen=NULL){

    roi.g=NULL
    if(length(roi)>0){
        if(length(grep('chr',roi))==0){
            if(length(gen)==0){
                stop('Please supply gencode details so that the 
                appropriate genomic region for your query genes can be derived')
            }
            roi.g=roi
            roi=gene2roi(roi,gen)
        }
    }
    if(length(roi)==0 & length(scn)>50){
        stop('For more than 50 samples, we recommend merging or agglomerating 
            samples for a query region or chr, 
            otherwise resulting data matrices will be too large')
    }

    #load up data
    usjcc=c('chr','start','end','JuncType','gene_name','FrameStatus')
    roi.nn=NULL
    if(length(roi)>1)
        roi.nn=as.numeric(roi[2:3])
    if(is.list(scn)){
        sj=scn
        names(sj)=names(scn)
        if(length(roi)>0)
            sj=lapply(sj,function(x) x[which(x[,'chr']==roi[1]),])
        if(length(roi)>1){
            sj=lapply(sj,function(x)
                x[which(as.numeric(x[,'start'])>=roi.nn[1]),])
            sj=lapply(sj,function(x)
                x[which(as.numeric(x[,'end'])<=roi.nn[2]),])
        }
    }
    if(is.character(scn)){
        sj=NULL
        for(i in seq(1,length(scn),1)){
            id=scn[i]
            tmp=gsub(' ','',as.matrix(read.delim(scn[i])))
            if(length(roi)>0)
                tmp=tmp[which(tmp[,'chr']==roi[1]),]
            if(length(roi.nn)>0){
                tmp=tmp[which(as.numeric(tmp[,'start'])>=roi.nn[1]),]
                tmp=tmp[which(as.numeric(tmp[,'end'])<=roi.nn[2]),]
            }
            sj[[i]]=tmp
        }
    }
    #adding rownames
    for(i in seq(1,length(sj),1)){
        tmp=sj[[i]][,usjcc]
        h=which(tmp=='')
        if(length(h)>0) tmp[h]='NA'
        tmp=gsub(' ','',tmp)
        rownames(sj[[i]])=apply(tmp,1,function(x) paste(x,collapse=' '))
    }

    #agglomerate data into matrices
    usj=unique(unlist(lapply(sj,function(x) rownames(x))))
    N=length(usj)
    RRS=matrix(0,N,length(scn))
    rownames(RRS)=usj
    colnames(RRS)=names(sj)
    NR=RRS
    #initializing MUTS matrix
    MUTS=NULL
    n1=unlist(lapply(sj,function(x) ncol(x)))
    n2=unlist(lapply(sj,function(x)
        max(which(is.element(colnames(x),c('FrameStatus','covRRS'))))))
    h=which((n1-n2)>0)
    if(length(h)>0){
        mm=lapply(h,function(x) sj[[x]][,(n2[x]+1):n1[x]])
        tmp=unique(unlist(strsplit(unlist(mm),'\\|')))
        tmp=setdiff(tmp,'')
        tmp=tmp[which(!is.na(tmp))]
        MUTS=matrix(0,length(tmp),length(scn))
        rownames(MUTS)=tmp
        colnames(MUTS)=colnames(RRS)
    }

    for(i in seq(1,length(scn),1)){
        sj.tmp=sj[[i]]
        xx=intersect(rownames(RRS),rownames(sj.tmp))
        RRS[xx,i]=as.numeric(sj.tmp[xx,'RRS'])
        NR[xx,i]=as.numeric(sj.tmp[xx,'uniq.reads'])
        if(length(MUTS)>0){
            q=max(which(is.element(colnames(sj.tmp),c('FrameStatus','covRRS'))))
            if(ncol(sj.tmp)>q){
                Q=(q+1):ncol(sj.tmp)
                mut.new=as.vector(sj.tmp[,Q])
                mut.new=mut.new[which(mut.new!='')]
                mut.new=mut.new[which(!is.na(mut.new))]
                mut.new=setdiff(unique(unlist(strsplit(mut.new,'\\|'))),'')
                MUTS[mut.new,i]=1
            }
        }
    }

    ##representative sample
    S0=t(matrix(unlist(strsplit(rownames(RRS),' ')),6,nrow(RRS)))
    if(method=='mean') S0=cbind(S0,apply(RRS,1,mean),apply(NR,1,mean))
    if(method=='median') S0=cbind(S0,apply(RRS,1,median),apply(NR,1,median))
    colnames(S0)=c('chr','start','end','JuncType',
        'gene_name','FrameStatus','RRS','uniq.reads')

    out=NULL
    if(length(roi)>0) out$roi=roi
    if(length(roi.g)>0) out$roi=c(roi,paste(roi.g,collapse=','))
    out$RRS=RRS
    out$NR=NR
    if(length(MUTS)>0) out$MUTS=MUTS
    out$SJ=S0
    rownames(out$SJ)=NULL
    return(out)
}

##get genomic coordinates for gene name
gene2roi<-function(g,gen){

    q=which(is.element(gen$GENES[,'gene_name'],g))
    tmp=intersect(g,gen$GENES[q,'gene_name'])
    if(length(tmp)==0){
        stop('Gene/s not found in the gencode object. 
            Please check your gene list')
    }
    if(length(tmp)!=length(g)){
        print(paste('Gene/s',setdiff(g,tmp),'are not listed in the 
            supplied gencode object and will be excluded.'))
    }
    if(length(unique(gen$GENES[q,'chr']))!=1){
        print(paste('The genes supplied fall into the following 
            chromosomes:',paste(unique(gen$GENES[q,'chr']),collapse=',')))
        stop('Only one chromosome at a time is accepted. 
            Please supply genes in the same chromosome')
    }

    roi=NULL
    tt=sort(table(gen$GENES[q,'chr']))
    q=intersect(q,which(gen$GENES[,'chr']==tail(names(tt),1)))
    if(length(unique(gen$GENES[q,'chr']))!=1){
        stop('Gene list supplied occur in more than one chromosome. 
            Please supply gene/s that occur in one chromosome only')
    }
    roi=range(as.numeric(as.vector(gen$GENES[q,c('start','end')])))
    roi=c(tail(names(tt)[1]),roi)

    return(roi)
}

