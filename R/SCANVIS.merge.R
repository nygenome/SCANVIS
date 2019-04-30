##Phaedra Agius, April 2019
##
##INPUT: scn=urls to scn files output by SCANVIS scan or linkvar functions
##            OR scn can also be a list of matrices
##        method='mean' or 'median'
##        roi=region of interest chr,start,end
##        gen=gencode data
##
##OUTPUT: object with components sj, sj.samples and muts (muts available only if 
##            scn are outputs of SCANVIS.linkvar)
##           out$sj=averaged RRS and uniq.reads across samples
##           out$sj.samples=comma sep samples, RRSs and uniq.reads for each sj
##           out$muts=all mutations mapped to sjs in out$sj, with comma sep
##           samples and number of samples containing the mutation
SCANVIS.merge<-function(scn,method='mean',roi=NULL,gen=NULL){

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

    RRS=matrix(0,2,2)
    rownames(RRS)=c('nada1','nada2')
    colnames(RRS)=rownames(RRS)
    NR=RRS
    MUTS=RRS

    for(i in seq(1,length(scn),1)){
        if(is.list(scn)){
            id=names(scn)[i]
            sj=scn[[i]]
        }
        if(is.character(scn)){
            id=scn[i]
            sj=as.matrix(read.delim(scn[i]))
        }
        if(length(roi)>0) sj=sj[which(sj[,'chr']==roi[1]),]
        if(length(roi)>1){
            q=which(as.numeric(sj[,'start'])>=as.numeric(roi[2]))
            q=q[which(as.numeric(sj[q,'end'])<=as.numeric(roi[3]))]
            sj=sj[q,]
        }
        sj=gsub(' ','',sj)
        rownames(sj)=paste(sj[,'chr'],sj[,'start'],sj[,'end'],sj[,'JuncType'],
            sj[,'gene_name'],sj[,'FrameStatus'])
        
        RRS=cbind(RRS,0)
        colnames(RRS)[ncol(RRS)]='RRS'
        NR=cbind(NR,0)
        colnames(NR)[ncol(NR)]='uniq.reads'
        dd=setdiff(rownames(sj),rownames(RRS))
        if(length(dd)>0){
            tmp=c(rownames(RRS),dd)
            RRS=rbind(RRS,matrix(0,length(dd),ncol(RRS)))
            rownames(RRS)=tmp
            NR=rbind(NR,matrix(0,length(dd),ncol(NR)))
            rownames(NR)=tmp
        }
        xx=intersect(rownames(sj),rownames(RRS))
        RRS[xx,'RRS']=as.numeric(sj[xx,'RRS'])
        NR[xx,'uniq.reads']=as.numeric(sj[xx,'uniq.reads'])

        q=max(which(is.element(colnames(sj),c('FrameStatus','covRRS'))))
        if(ncol(sj)>q){
            Q=(q+1):ncol(sj)
            mut.new=sj[,Q] 
            #MUTS=addmut(MUTS,sj[,Q],id)
            MUTS=cbind(MUTS,0)
            colnames(MUTS)[ncol(MUTS)]=id
            mut.new=as.vector(mut.new)
            mut.new=mut.new[which(mut.new!='')]
            mut.new=mut.new[which(!is.na(mut.new))]
            mm=setdiff(unique(unlist(strsplit(mut.new,'\\|'))),'')
            if(length(mm)>0){
                dd=setdiff(mm,rownames(MUTS))
                if(length(dd)>0){
                    tmp=c(rownames(MUTS),dd)
                    MUTS=rbind(MUTS,matrix(0,length(dd),ncol(MUTS)))
                    rownames(MUTS)=tmp
                }
                xx=intersect(mm,rownames(MUTS))
                MUTS[xx,id]=1
            }
        }
    }

    RRS=RRS[3:nrow(RRS),3:ncol(RRS)]
    NR=NR[3:nrow(NR),3:ncol(NR)]
    if(nrow(MUTS)==2) MUTS=NULL
    if(length(MUTS)>0){
        if(nrow(MUTS)>3)
            MUTS=MUTS[3:nrow(MUTS),3:ncol(MUTS)]
        if(nrow(MUTS)==3){
            tmp=rownames(MUTS)[3]
            MUTS=t(as.matrix(MUTS[3,3:ncol(MUTS)]))
            rownames(MUTS)=tmp
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

# deprecated function
# add2mat<-function(mat,sj.new,id,whattoadd='RRS'){
#     mat=cbind(mat,0)
#     colnames(mat)[ncol(mat)]=id
#     dd=setdiff(rownames(sj.new),rownames(mat))
#     if(length(dd)>0){
#         tmp=c(rownames(mat),dd)
#         mat=rbind(mat,matrix(0,length(dd),ncol(mat)))
#         rownames(mat)=tmp
#     }
#     xx=intersect(rownames(sj.new),rownames(mat))
#     mat[xx,id]=as.numeric(sj.new[xx,whattoadd])
#     return(mat)
# }
