##Phaedra Agius, April 2019
##
###SCANVIS.scan - this function scores and annotates splice junctions
##
#INPUT sj=splice junction matrix with colnames chr,start,end,uniq.reads 
##           OR 
##         url to file
##      gen=output or URL to output from SCANVIS.annotation
##      Rcut = (Default=5) only SJs with >=Rcut reads are processed
##      bam = (Default=NULL) url to bamfile (recommended for RRC scores for NEs)
##      samtools = URL to your samtools executable
##
##OUTPUT sj file from INPUT with additional columns indicating the following:
##           gene_name=names of any genes with intervals that overlap SJ
##            JuncType: annot, exon.skip, alt5, alt3, IsoSwitch, Unknown or NE
##           RRS=Relative Read Support score
##           genomic_interval = interval used to compute RRS score and defined 
##                              s.t. it contains >=1 gene overlapping SJ 
##                              and >=1 annotated SJ
##           RRC=Relative Read Coverage for NEs only, 
##               computed only if bam is supplied
##           For NE events, additional rows are introduced and 
##           the start/end coordinates correspond to two separate SJs
##
##Example: scn=SCANVIS.scan('~/sj.tab',gen25,5)
SCANVIS.scan<-function(sj,gen,Rcut=5,bam=NULL,samtools=NULL){

    if(length(bam)>0 & length(samtools)==0)
        stop('Please specify your path to samtools')

    if(length(gen)==1) gen=get(load(gen))
    gen.EXONS=gen$EXONS
    gen.GENES=gen$GENES
    gen.GENES.merged=gen$GENES.merged
    gen.INTRONS=gen$INTRONS
    rm(gen)
    
    if(length(sj)==1)
        sj=as.matrix(read.delim(sj))
    if(length(intersect(colnames(sj),c('chr','start','end','uniq.reads')))!=4){
        stop('Error! INPUT sj must have at least 4 columns with chr, 
            start, end and uniq.reads as colnames')
    }
    sj=sj[which(as.numeric(sj[,'uniq.reads'])>=Rcut),]
    sj=gsub(' ','',sj)
    sj=sj[which(is.element(sj[,'chr'],unique(gen.GENES[,'chr']))),]
    d=as.numeric(sj[,'end'])-as.numeric(sj[,'start'])
    if(min(abs(d))==0)
        sj=sj[which(abs(d)>0),]
    if(min(d)<0){
        h=which(d<0)
        sj[h,c('start','end')]=sj[h,c('end','start')]
    }
    tmp=intersect(sj[,'chr'],gen.GENES[,'chr'])
    if(length(tmp)==0){
        sj[,'chr']=paste0('chr',sj[,'chr'])
        tmp=intersect(sj[,'chr'],gen.GENES[,'chr'])
    }
    sj=sj[which(is.element(sj[,'chr'],tmp)),] ##only keep chr present in gen
    ##sort sj by start position
    for(i in unique(sj[,'chr'])){
        q=which(sj[,'chr']==i)
        sj[q,]=sj[q[order(as.numeric(sj[q,'start']))],]
    }

    ############################################################################
    ##START: annotate SJs, ASJ and USJ
    sj=gsub(' ','',sj)
    gen=gen.EXONS
    options(scipen=999)
    out=rep(0,nrow(sj))
    q=which(as.numeric(sj[,'end'])<as.numeric(sj[,'start']))
    if(length(q)>0)
            sj[q,c('start','end')]=sj[q,c('end','start')]
    sj.tmp=paste(sj[,'chr'],sj[,'start'],sj[,'end'])
    Q=which(gen[,'type']=='exon')
    q=intersect(Q,which(gen[,'strand']=='+'))
    tx=as.vector(t(cbind(gen[q,'transcript'],gen[q,'transcript'])))
    tmp=as.vector(t(cbind(gen[q,'chr'],gen[q,'chr'])))
    tmp2=as.vector(t(cbind(as.numeric(gen[q,'start'])-1,
        as.numeric(gen[q,'end'])+1)))       
    tx=t(matrix(tx[2:(length(tx)-1)],2,length(q)-1))
    tmp=t(matrix(tmp[2:(length(tmp)-1)],2,length(q)-1))
    tmp2=t(matrix(tmp2[2:(length(tmp2)-1)],2,length(q)-1))
    h=which(tx[,1]==tx[,2])
    tmp=tmp[h,]
    tmp2=tmp2[h,]
    annot.tmp=paste(tmp[,1],tmp2[,1],tmp2[,2])
    z=intersect(sj.tmp,annot.tmp)
    q=which(is.element(sj.tmp,z))
    out[q]=1

    q=intersect(Q,which(gen[,'strand']=='-'))
    tx=t(matrix(rev(as.vector(t(cbind(gen[q,'transcript'],
        gen[q,'transcript'])))),2,length(q)))
    tmp=t(matrix(rev(as.vector(t(cbind(gen[q,'chr'],
        gen[q,'chr'])))),2,length(q)))
    tmp2=t(matrix(rev(as.vector(t(cbind(as.numeric(gen[q,'start']),
        as.numeric(gen[q,'end']))))),2,length(q)))
    tmp3=cbind(tmp2[,2]-1,tmp2[,1]+1)
    z=c(1,length(tmp3))
    tx=matrix(tx[-z],length(q)-1,2)
    tmp=matrix(tmp[-z],length(q)-1,2)
    tmp3=matrix(tmp3[-z],length(q)-1,2)
    annot.tmp=paste(tmp[,1],tmp3[,2],tmp3[,1])
    #annot.tmp=c(annot.tmp,paste(tmp[,1],tmp3[,1],tmp3[,2]))
    z=intersect(sj.tmp,annot.tmp)
    q=which(is.element(sj.tmp,z))
    out[q]=1
    annot.tmp=paste(tmp[,1],tmp3[,1],tmp3[,2])
    z=intersect(sj.tmp,annot.tmp)
    q=which(is.element(sj.tmp,z))
    out[q]=1

    out[which(out==1)]='annot'
    sj=cbind(sj,out)
    colnames(sj)[ncol(sj)]='JuncType'
    gen=NULL

    ############################################################################
    ##START:describe USJs
    print('*** Categorizing unannotated splice junctions ... ***')
    Q=which(sj[,'JuncType']!='annot')
    if(length(Q)>0){
    	sj.orig=sj
    	sj=sj[Q,]
	    out=rep('Unknown',length(Q))
	    genEXONS=gen.EXONS[which(gen.EXONS[,'type']=='exon'),]
	    CHR=unique(sj[,'chr'])
	    for(chr in CHR){
	        q=which(sj[,'chr']==chr)
	        sj.tmp=sj[q,]
	        sj.coor=cbind((as.numeric(sj[q,'start'])-1),
	            (as.numeric(sj[q,'end'])+1))
	        if(length(q)==1)
	            sj.tmp=t(as.matrix(sj.tmp))
	        out.tmp=rep('',nrow(sj.tmp))

	        h=which(genEXONS[,'chr']==chr)
	        gen.tmp=genEXONS[h,]
	        gen.coor=cbind(as.numeric(genEXONS[h,'start']),
                as.numeric(genEXONS[h,'end']))
	        gen.strand=unique(genEXONS[h,c('transcript','strand')])
	        if(length(h)>1){
	            rownames(gen.strand)=gen.strand[,1]
	            gen.strand=gen.strand[,2]
	        }
	        if(length(h)==1) gen.strand=gen.strand[2]

	        B=cbind(1*is.element(sj.coor[,1],gen.coor),
	            1*is.element(sj.coor[,2],gen.coor))
	        if(length(h)>1){
	            tmp=c(gen.tmp[,'transcript'],gen.tmp[,'transcript'])
	            names(tmp)=c(gen.tmp[,'start'],gen.tmp[,'end'])
	        }
	        if(length(h)==1){
	            tmp=rep(genEXONS[h,'transcript'],2)
	            names(tmp)=genEXONS[h,c('start','end')]
	        }
	        z=which(B[,1]==1)
	        B[z,1]=tmp[as.character(sj.coor[z,1])]
	        z=which(B[,2]==1)
	        B[z,2]=tmp[as.character(sj.coor[z,2])]

	        BS=cbind(1*is.element(sj.coor[,1],gen.coor[,1]),
	            1*is.element(sj.coor[,2],gen.coor[,1]))
	        BE=cbind(1*is.element(sj.coor[,1],gen.coor[,2]),
	            1*is.element(sj.coor[,2],gen.coor[,2]))

	        #exon skip or isoform jumping
	        k=intersect(which(B[,1]!=0),which(B[,2]!=0))
	        if(length(k)>0){

	            kk=k[which(B[k,1]==B[k,2])]
	            if(length(kk)>0){
	                BX=B
	                tmp=c(gen.tmp[,'exon_number'],gen.tmp[,'exon_number'])
	                names(tmp)=c(gen.tmp[,'start'],gen.tmp[,'end'])
	                BX[kk,1]=tmp[as.character(sj.coor[kk,1])]
	                BX[kk,2]=tmp[as.character(sj.coor[kk,2])]
	                d=abs(apply(cbind(as.numeric(BX[kk,1]),
	                    as.numeric(BX[kk,2])),1,diff))
	                kkk=kk[which(d>1)]
	                if(length(kkk)>0)
	                    out.tmp[kkk]='exon.skip'
	                kkk=kk[which(d<=1)]
	                if(length(kkk)>0){
	                    if(length(kkk)>1){
	                        kkk.s=kkk[which(apply(BS[kkk,],1,sum)==2)]
	                        if(length(kkk.s)>0)
	                            out.tmp[kkk.s]='exon.skip.55'
	                        kkk.e=kkk[which(apply(BE[kkk,],1,sum)==2)]
	                        if(length(kkk.e)>0)
	                            out.tmp[kkk.e]='exon.skip.33'
	                    }
	                    if(length(kkk)==1){
	                        tmp=c(sum(BS[kkk,]),sum(BE[kkk,]))
	                        if(tmp[1]==2) out.tmp[kkk]='exon.skip.55'
	                        if(tmp[2]==2) out.tmp[kkk]='exon.skip.33'

	                    }
	                }
	            }
	            kk=k[which(B[k,1]!=B[k,2])]
	            if(length(kk)>0)
	                out.tmp[kk]='IsoSwitch'
	        }

	        #alt 3p
	        k=intersect(which(B[,1]!=0),which(B[,2]==0))
	        if(length(k)>0){
	            #paying attention to orientation of the gene
	            kk=k[which(gen.strand[B[k,1]]=='+')]
	            if(length(kk)>0) out.tmp[kk]='alt3p'
	            kk=k[which(gen.strand[B[k,1]]=='-')]
	            if(length(kk)>0) out.tmp[kk]='alt5p'
	        }
	        #alt5p
	        k=intersect(which(B[,1]==0),which(B[,2]!=0))
	        if(length(k)>0){
	            #paying attention to orientation of the gene
	            kk=k[which(gen.strand[B[k,2]]=='+')]
	            if(length(kk)>0) out.tmp[kk]='alt5p'
	            kk=k[which(gen.strand[B[k,2]]=='-')]
	            if(length(kk)>0) out.tmp[kk]='alt3p'
	        }
	        out[q]=out.tmp
	    }
	    h=which(out=='')
	    if(length(h)>0) out[h]='Unknown'
	    sj=sj.orig
	    sj[Q,'JuncType']=out
	    rm(sj.orig)
	}
    print('*** DONE: Categorizing unannotated splice junctions ***')
    ##END:describe USJs
    ############################################################################

    ##END: annotate SJs, ASJ and USJ
    ############################################################################

    ############################################################################
    ##START: compute RRS and assign gene names
    print('*** Computing RRS by median of local ASJ reads ... ***')
    #sj=sj2RRS(sj,gen.GENES)
    RRS=rep(-1,nrow(sj))
    gi=matrix('na',nrow(sj),2)
    colnames(gi)=c('gene_name','genomic_interval')
    q=which(sj[,'JuncType']=='annot')
    u=sort(table(sj[q,'chr']))
    u=u[which(u>0)]
    for(chr in names(u)){
    	#################################################################
    	##START: get genomic intervals o'lap SJs s.t. >=1 gene and >=1 ASJ
        #INTERVALS=getIntervals(sj,gen.GENES,chr)
	    gen.tmp=gen.GENES[which(gen.GENES[,'chr']==chr),]
	    gI=IRanges(as.numeric(gen.tmp[,'start']),as.numeric(gen.tmp[,'end']))

	    ##define intervals on ASJ only
	    q=which(sj[,'chr']==chr)
	    sjI=IRanges(as.numeric(sj[q,'start']),as.numeric(sj[q,'end']))
	    qA=intersect(q,which(sj[,'JuncType']=='annot'))
	    qU=intersect(q,which(sj[,'JuncType']!='annot'))

	    aI=IRanges(as.numeric(sj[qA,'start']),as.numeric(sj[qA,'end']))
	    if(length(aI)==0){
	        ##all SJs are USJs ... happens in smaller chr like chrM
	        INTERVALS=range(IR2Mat(gI))        
	        INTERVALS=reduce(IRanges(INTERVALS[1],INTERVALS[2]))
	    }
	    ##merge annot intervals with gene intervals
	    if(length(aI)>0){
	        v=as.matrix(findOverlaps(aI,gI))
	        gI.tmp=IR2Mat(gI[unique(v[,2]),])
	        INTERVALS=reduce(IRanges(gI.tmp[,1],gI.tmp[,2]))
	    }
	    ##make sure all the unannotated events are included in one interval
	    if(length(qU)>0){
	        uI=IRanges(as.numeric(sj[qU,'start']),as.numeric(sj[qU,'end']))
	        v=as.matrix(findOverlaps(uI,INTERVALS))
	        h=setdiff(seq(1,length(qU),1),v[,1])
	        if(length(h)>0){
	            uI=IR2Mat(uI)
	            INTERVALS=IR2Mat(INTERVALS)
	            for(j in h){
	                d1=cbind(abs(INTERVALS[,1]-uI[j,1]),
	                    abs(INTERVALS[,2]-uI[j,1]))
	                d2=cbind(abs(INTERVALS[,1]-uI[j,2]),
	                    abs(INTERVALS[,2]-uI[j,2]))
	                dd=apply(cbind(d1,d2),1,min)
	                k=which(dd==(min(dd)))[1]
	                INTERVALS[k,]=range(c(INTERVALS[k,],uI[j,]))
	            }
	            if(!is.matrix(INTERVALS)) INTERVALS=t(as.matrix(INTERVALS))
	            INTERVALS=IRanges(INTERVALS[,1],INTERVALS[,2])
	        }
	    }
	    v=as.matrix(findOverlaps(sjI,INTERVALS))
	    INTERVALS=IR2Mat(INTERVALS)
	    if(!is.matrix(INTERVALS)) INTERVALS=t(as.matrix(INTERVALS))
	    if(length(unique(v[,1]))<nrow(v)){
	        tt=table(v[,1])
	        tt=tt[which(tt>1)]
	        for(j in as.numeric(names(tt))){
	            I.tmp=range(INTERVALS[v[which(v[,1]==j),2]])
	            INTERVALS=rbind(INTERVALS,I.tmp)
	            v[which(v[,1]==j),2]=nrow(INTERVALS)
	        }
	        v=unique(v)
	    }
	    INTERVALS=INTERVALS[v[,2],]
	    if(!is.matrix(INTERVALS)) INTERVALS=t(as.matrix(INTERVALS))
	    INTERVALS=IRanges(INTERVALS[,1],INTERVALS[,2])

    	##END: get genomic intervals o'lap SJs s.t. >=1 gene and >=1 ASJ
    	###############################################################

        IIu=unique(IR2Mat(INTERVALS))
        IIr=reduce(IRanges(IIu[,1],IIu[,2]))
        if(nrow(as.matrix(IIr))<nrow(IIu)){
            v=as.matrix(findOverlaps(INTERVALS,IIr))
            if(nrow(v)!=nrow(as.matrix(INTERVALS))) 
                stop('Error with row mismatch in sj2RRS')
            INTERVALS=IIr[v[,2],]
        }

        II=IR2Mat(INTERVALS)
        gi.tmp=paste0(chr,':',II[,1],'-',II[,2])
        q=which(sj[,'chr']==chr)
        gi[q,'genomic_interval']=gi.tmp
        ugi=unique(gi.tmp)
        Q=lapply(ugi,function(x) which(gi.tmp==x))        
        qA=which(sj[q,'JuncType']=='annot')
        Q=lapply(Q,function(x) intersect(x,qA))
        nA=unlist(lapply(Q,function(x) length(x)))
        if(sum(nA==0)>0){ 
            stop('Error with allocating intervals 
                encompassing at least one annotated SJ')
        }

        mm=unlist(lapply(Q,
            function(x) median(as.numeric(sj[x,'uniq.reads']))))
        names(mm)=ugi
        mm=mm[gi.tmp]
        RRS.tmp=as.numeric(sj[q,'uniq.reads'])
        RRS[q]=RRS.tmp/(RRS.tmp+mm)

        #get gene names
        k=which(gen.GENES[,'chr']==chr)
        gentmp=gen.GENES[k,c('start','end')]
        gentmp=cbind(as.numeric(gentmp[,1]),as.numeric(gentmp[,2]))
        v=as.matrix(findOverlaps(INTERVALS,IRanges(gentmp[,1],gentmp[,2])))
        v1=v[,1]
        v[,2]=gen.GENES[k[v[,2]],'gene_name']
        gtmp=unlist(lapply(unique(v1),function(x) 
            paste(sort(unique(v[which(v1==x),2])),collapse=',')))
        gi[q,'gene_name']=gtmp
    }
    tmp=colnames(sj)
    sj=cbind(sj,RRS,gi[,'genomic_interval'])
    colnames(sj)=c(tmp,'RRS','genomic_interval')
    print('*** DONE: Computing RRS by median of local ASJ reads ***')
    ##END: compute RRS and assign gene names
    ############################################################################

    ############################################################################
    #START: assigning gene names
    print('*** Designating gene names to sj coordinates ... ***')
    gn=rep('intergenic',nrow(sj))
    ##merge gen.GENES first
    CHR=unique(sj[,'chr'])
    for(chr in CHR){
        I1=which(sj[,'chr']==chr)
        IR1=IRanges(as.numeric(sj[I1,'start']),as.numeric(sj[I1,'end']))
        I2=which(gen.GENES[,'chr']==chr)
        IR2=IRanges(as.numeric(gen.GENES[I2,'start']),
            as.numeric(gen.GENES[I2,'end']))
        v=as.matrix(findOverlaps(IR1,IR2))
        tmp=gen.GENES[I2[v[,2]],'gene_name']
        tt=table(v[,1])
        if(max(tt)>0){
            n=as.numeric(names(tt)[which(tt>1)])
            for(j in n){
                q=which(v[,1]==j)
                tmp[q]=paste(tmp[q],collapse=',')
            }
        }
        tmp=unique(cbind(v[,1],tmp))
        gn[I1[as.numeric(tmp[,1])]]=tmp[,2]
    }
    sj=cbind(sj,gn)
    colnames(sj)[ncol(sj)]='gene_name'
    q=which(sj[,'gene_name']=='')
    if(length(q)>0)
        sj[q,'gene_name']='NONE'
    #for multiple genes, sort in lexicographic order
    q=grep(',',sj[,'gene_name'])
    if(length(q)>0)
        sj[q,'gene_name']=unlist(lapply(sj[q,'gene_name'],
            function(x) paste(sort(unique(unlist(strsplit(x,',')))),
                collapse=',')))

    print('*** DONE: Designating gene names to sj coordinates ***')
    ##END: assigning gene names
    ############################################################################

    ############################################################################
    ##START: Querying inFrameStatus
    print('*** Querying inFrameStatus ... ***')
    sj=cbind(sj,rep('NA',nrow(sj)))
    colnames(sj)[ncol(sj)]='FrameStatus'
    QFS=which(sj[,'JuncType']!='annot')
    if(length(QFS)>0){
    	sj.tmp=sj[QFS,]
	    FS=rep('Unknown',nrow(sj.tmp))
	    for(chr in unique(sj.tmp[,'chr'])){
	        q=which(sj.tmp[,'chr']==chr)
	        q=intersect(q,which(sj.tmp[,'JuncType']!='NE'))
	        q=intersect(q,which(sj.tmp[,'JuncType']!=''))
	        q=intersect(q,union(grep('exon.skip',sj.tmp[,'JuncType']),
	            grep('alt',sj.tmp[,'JuncType'])))
	        if(length(q)>0){
	            gen.chr=gen.EXONS[which(gen.EXONS[,'chr']==chr),]
	            gen.chr.se=cbind(as.numeric(gen.chr[,'start']),
	                as.numeric(gen.chr[,'end']))
	            
	            for(j in q){
	                v=as.numeric(sj.tmp[j,'start']):as.numeric(sj.tmp[j,'end'])
	                sj.se=c(v[1]-1,tail(v,1)+1)
	                h=which(gen.chr.se==sj.se[1],arr.ind=TRUE)[,1]
	                h=c(h,which(gen.chr.se==sj.se[2],arr.ind=TRUE)[,1])
	                TX=unique(gen.chr[h,'transcript'])
	                out.tmp=NULL
	                for(tx in TX){
	                    d=NULL
	                    k=which(gen.chr[,'transcript']==tx)
	                    h.tx=h[which(gen.chr[h,'transcript']==tx)]
	                    #exon skipping events
	                    if(length(grep('exon.skip',sj.tmp[j,'JuncType']))>0){
	                        tmp=unlist(lapply(k,
	                            function(x) gen.chr.se[x,1]:gen.chr.se[x,2]))
	                        d=length(intersect(tmp,v))
	                    }
	                    #alt5/3p events
	                    if(length(grep('alt',sj.tmp[j,'JuncType']))>0){
	                        m=setdiff(sj.se,gen.chr.se[k,])
	                        if(length(m)>1)
	                          stop('Both ends of an alt5/3p event do not align')
	                        inExon=FALSE
	                        z=which(gen.chr.se[k,1]<m)    
	                        if(length(z)>0)
	                            z=intersect(z,which(gen.chr.se[k,2]>m))
	                        if(length(z)>0)
	                            inExon=TRUE
	                        if(inExon){
	                            tmp=unlist(lapply(k,function(x) 
	                                gen.chr.se[x,1]:gen.chr.se[x,2]))
	                            d=length(intersect(tmp,v))
	                        }
	                        if(!inExon){ #adding nt
	                            d0=cbind(gen.chr.se[k,1]-m,gen.chr.se[k,2]-m)
	                            da=apply(abs(d0),1,min)
	                            d=min(da)-1
	                        }
	                    }

	                    tt=floor(d/3)==(d/3)
	                    if(tt) out.tmp=rbind(out.tmp,c(tx,'INframe'))
	                    if(!tt) out.tmp=rbind(out.tmp,c(tx,'OUTframe'))
	                }

	                if(length(out.tmp)==2) FS[j]=out.tmp[2]
	                if(length(out.tmp)>2){
	                    u=unique(out.tmp[,2])
	                    if(length(u)==1) FS[j]=u
	                    if(length(u)>1){
	                        FS[j]=paste(unlist(lapply(seq(1,nrow(out.tmp),1),
	                            function(x) paste(out.tmp[x,1],
	                                out.tmp[x,2],sep=':'))),collapse=',')
	                    }
	                }
	            }
	        }
	    }
        sj[QFS,'FrameStatus']=FS
        print('*** DONE: Querying FrameStatus ***')
    }
    ##END: Querying inFrameStatus
    ############################################################################

    ############################################################################
    ##START: getting/scoring NEs by finding USJs coinciding in intronic regions
    print('*** Collecting and scoring all potential NEs ... ***')
    NE=NULL    
    USJ=gsub(' ','',sj[which(sj[,'JuncType']!='annot'),])
    CHR=intersect(unique(USJ[,'chr']),gen.INTRONS[,'chr'])
    if(length(CHR)==0) USJ[,'chr']=paste0('chr',USJ[,'chr'])
    CHR=intersect(unique(USJ[,'chr']),gen.INTRONS[,'chr'])
    colIDs=c('uniq.reads',colnames(USJ)[grep('RRS',colnames(USJ))])
    QE=which(as.numeric(USJ[,'start'])==as.numeric(USJ[,'end']))
    for(chr in CHR){
        #print(chr)
        INtmp=gen.INTRONS[which(gen.INTRONS[,'chr']==chr),]
        q=setdiff(which(USJ[,'chr']==chr),QE)
        sj.coor=cbind(as.numeric(USJ[q,'start']),as.numeric(USJ[q,'end']))
        NREADS=as.numeric(USJ[q,'uniq.reads'])
        RRS=as.numeric(USJ[q,'RRS'])

        h=which(sj.coor[,2]-sj.coor[,1]<0)
        if(length(h)>0) sj.coor=sj.coor[h,]
        if(!is.matrix(sj.coor)) sj.coor=t(as.matrix(sj.coor))
        I=unique(sj.coor[,1])
        I=cbind(I,unlist(lapply(I,function(x) 
            sum(NREADS[which(sj.coor[,1]==x)]))))
        I=cbind(I,unlist(lapply(I[,1],function(x) 
            mean(RRS[which(sj.coor[,1]==x)]))))

        J=unique(sj.coor[,2])
        J=cbind(J,unlist(lapply(J,function(x) 
            sum(NREADS[which(sj.coor[,2]==x)]))))
        J=cbind(J,unlist(lapply(J[,1],function(x) 
            mean(RRS[which(sj.coor[,2]==x)]))))
        J=J[which(J[,2]>0),]

        if(length(I)>0 & length(J)>0){
            if(!is.matrix(J) & length(J)==3) J=t(as.matrix(J))                        
            for(k in seq(1,nrow(I),1)){
                #print(c(k,nrow(I)))
                j=J
                j=j[which(j[,1]<I[k,1]),] #end_junct plus start_junc
                if(!is.matrix(j) & length(j)==3) j=t(as.matrix(j))                        
                if(length(j)>0){
                    #find all intronic ranges containing end_junc (i)
                    z=intersect(which(as.numeric(INtmp[,'start'])<I[k,1]),
                        which(as.numeric(INtmp[,'end'])>I[k,1]))
                    if(length(z)>0){
                        INtmp.z=cbind(as.numeric(INtmp[z,'start']),
                            as.numeric(INtmp[z,'end']))
                        if(length(j)>3){
                           zz=unique(unlist(lapply(seq(1,length(z),1),
                            function(x) intersect(which(j[,1]>INtmp.z[x,1]),
                                which(j[,1]<INtmp.z[x,2])))))
                       }
                        if(length(j)==3){
                            zz=NULL
                            if(length(intersect(which((j[1]-INtmp.z[,1])>0),
                            	which((INtmp.z[,2]-j[1])>0)))>0)
                                zz=1
                        }
                        if(length(zz)>0){
                            if(length(zz)==1){
                             NE=rbind(NE,c(chr,I[k,1],j[zz,1],
                                 mean(c(I[k,2],j[zz,2])),
                                 mean(c(I[k,3],j[zz,3]))))
                            }
                            if(length(zz)>1){
                             NE=rbind(NE,cbind(chr,I[k,1],
                                 j[zz,1],apply(cbind(I[k,2],j[zz,2]),1,mean),
                                 apply(cbind(I[k,3],j[zz,3]),1,mean)))
                            }
                        }
                    }
                }
            }
        }
    }
    if(length(NE)>0){
        if(!is.matrix(NE)) NE=t(as.matrix(NE))
        NE[,2:3]=NE[,c(3,2)]
        if(length(NE)>4)
            colnames(NE)=c('chr','start','end','uniq.reads','RRS')

        ##print('*** Annotating by gene name ***')
        g=rep('',nrow(NE))
        for(i in unique(NE[,'chr'])){
            q=which(NE[,'chr']==i)
            IR1=IRanges(as.numeric(NE[q,'start']),as.numeric(NE[q,'end']))
            h=which(gen.GENES[,'chr']==i)
            IR2=IRanges(as.numeric(gen.GENES[h,'start']),
                as.numeric(gen.GENES[h,'end']))
            v=as.matrix(findOverlaps(IR1,IR2))
            v[,2]=gen.GENES[h[v[,2]],'gene_name']
            if(length(v)>0){
                tt=table(v[,1])
                if(max(tt)>1){
                    for(j in names(tt)[which(tt>1)]){
                        k=which(v[,1]==j)
                        u=unique(v[k,2])
                        v[k,2]=paste(u,collapse=',')
                    }
                    v=unique(v)
                }
                g[q[as.numeric(v[,1])]]=v[,2]
            }
        }
        NE=cbind(NE,g)
        colnames(NE)[ncol(NE)]='gene_name'
        NE=cbind(NE,0)
        colnames(NE)[ncol(NE)]='annotated'
        NE=cbind(NE,'NE')
        colnames(NE)[ncol(NE)]='JuncType'
        h=which(g=='')
        if(length(h)>0){
            q=which(USJ[,'gene_name']=='')
            if(length(q)>0) USJ[q,'gene_name']='NA'
            NE[h,'JuncType']='NE_IG'
            rownames(USJ)=paste(USJ[,'chr'],USJ[,'end'])
            tmp=paste(NE[h,'chr'],NE[h,'start'])
            g=USJ[tmp,'gene_name']
            rownames(USJ)=paste(USJ[,'chr'],USJ[,'start'])
            tmp=paste(NE[h,'chr'],NE[h,'end'])
            g=paste(g,USJ[tmp,'gene_name'],sep='---')
            NE[h,'gene_name']=g
        }
    }
    ##return only NEs with at least 3bp
    if(length(NE)>0){
    	NE=gsub(' ','',NE)
        dd=as.numeric(NE[,'end'])-as.numeric(NE[,'start'])
        q=which(dd>=3)
        ##filter for any NEs that are actually Unknown SJs
        tmp=apply(USJ[,c(1,2,3)],1,function(x) paste(x,collapse=' '))
        if(length(NE)>8)
            tmp2=apply(NE[,c(1,2,3)],1,function(x) paste(x,collapse=' '))
        if(length(NE)==8) 
            tmp2=paste(NE[c(1,2,3)],collapse=' ')
        xx=intersect(tmp2,tmp)
        if(length(xx)>0) q=intersect(q,which(!is.element(tmp2,xx)))
        NE=NE[q,]
    }
    ############################################################################
    #START: assigning RRS scores to NEs
    if(length(NE)>0){
        if(!is.matrix(NE)) NE=t(as.matrix(NE))
        print('*** Computing RRCs using bam file ***')
        if(length(bam)>0){
            BED=NE[,c('chr','start','end')]
            p=0.2
		    if(!is.matrix(BED)) BED=t(as.matrix(BED))
		    n=as.numeric(BED[,3])-as.numeric(BED[,2])+1
		    n=ceiling(n*(p/2))
		    bed5=cbind(BED[,1],as.numeric(BED[,2])-n,as.numeric(BED[,2])-1)
		    bed3=cbind(BED[,1],as.numeric(BED[,3])+1,as.numeric(BED[,3])+n)

		    BED=rbind(BED,bed5,bed3)
		    id=unlist(strsplit(tail(unlist(strsplit(bam,'/')),1),'\\.'))[1]
		    fbam=paste0('tmp_',id)
		    write.table(bam,fbam,sep='\t',quote=FALSE,
		        row.names=FALSE,col.names=FALSE)
		    sam.bed=paste0(BED[,1],':',BED[,2],'-',BED[,3])
			nn=NULL
			cmd=paste("-r sam.bed -f",fbam)
            cmd=paste(cmd,"| awk '{ sum += $3 } END { print sum }'")
			cmd=paste(samtools,'depth',cmd)
			for(i in seq(1,length(sam.bed),1)){
		        fout=paste0('tmp_',id,'__',i)
		        cmd.tmp=gsub('sam.bed',sam.bed[i],cmd)
		        cmd.tmp=paste(cmd.tmp,'>',fout)
		        system(cmd.tmp)
			    if(file.info(fout)$size>1)
		            nn=c(nn,as.matrix(read.delim(fout,header=FALSE)))
		        system(paste('rm',fout))
			}
		    system(paste('rm',fbam))
            nn=as.numeric(nn)
		    ii=1
		    covB=nn[ii:nrow(BED)]
		    ii=nrow(BED)+1
		    cov5p=nn[ii:(ii+nrow(bed5)-1)]
		    ii=ii+nrow(bed5)
		    cov3p=nn[ii:length(nn)]
		    
		    RRC=covB/(as.numeric(BED[,3])-as.numeric(BED[,2]))    
		    RRC=cbind(RRC,(covB/(covB+cov5p+cov3p)))
		    NE=cbind(NE,RRC[,2])
		    colnames(NE)[ncol(NE)]='RRC'
        }
    }
    #END: assigning RRS scores to NEs
    ############################################################################
    if(length(NE)>0){
        tmp=matrix(NA,nrow(NE),ncol(sj))
        colnames(tmp)=colnames(sj)
        cc=c('chr','start','end','uniq.reads','gene_name','RRS')
        tmp[,cc]=NE[,cc]
        tmp[,'JuncType']='NE'
        if(length(bam)==0) sj=rbind(sj,tmp)
        if(length(bam)>0){
            sj=rbind(cbind(sj,0),cbind(tmp,NE[,'RRC']))
            colnames(sj)[ncol(sj)]='RRC'
        }    
    }
    print('*** DONE: Collecting and scoring all potential NEs ... ***')
    ##END: getting/scoring NEs by finding USJs coinciding in intronic regions
    ############################################################################

    return(sj)
}

#function for converting IRanges intervals to matrix
IR2Mat<-function(I){
    tmp=as.matrix(I)
    tmp[,2]=tmp[,2]+tmp[,1]-1
    return(tmp)
}



