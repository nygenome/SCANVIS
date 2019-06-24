##Phaedra Agius, April 2019
##
##INPUT: roi=region of interest either 
###          a 3 bit vector (chr,start,end) 
###          OR a single gene name 
###          OR a character vector of gene names
###      gen=gencode object as output by SCANVISannotation.R script
###      scn=one or more urls to scn files (if multiple then these will be merged) 
###              OR a scn file in matrix format OR output from SCANVISmerge
###      SJ.special=default is NULL, otherwise specify SJs of interest as a 3 
##        col matrix with chr,start,end which will be shown as cyan-colored arcs
##       TITLE=default is '', otherwise specify figure title
##       bam=default is NULL, otherwise specify url to bam for read profile plot
##       samtools=default is NULL, otherwise path to your samtools, 
##                MUST be specified if bam is not NULL (default=NULL)
##       full.annot=TRUE for full annotation details, 
##                  FALSE for condensed format (default=FALSE)
##       USJ=default is "NR", this parameter indicates whether to print 
##           NR=num reads OR "RRS"=Relative Read Support at top of USJ arcs
##
##OUTPUT: an object containing SJs plotted, any mutations plotted and 
##       coordinates of the region plotted

SCANVISvisual<-function(roi,gen,scn,SJ.special=NULL,TITLE=NULL,bam=NULL,
    samtools=NULL,full.annot=FALSE,USJ="NR"){

    if(length(bam)>0 & length(samtools)==0){
        stop('For read coverage plots with bam, 
            you must supply your url to samtools')
    }

    roi.g=''
    if(length(grep('chr',roi))==0){
        roi.g=roi
        roi=gene2roi(roi,gen)
    }
    chr=roi[1]
    XLIM=as.numeric(roi[2:3])
    print(paste('*** Sashimi plot  spanning',diff(XLIM),
        'base pairs is being generated'))

    muts=NULL
    if(is.matrix(scn)) sj=scn
    if(!is.matrix(scn)){
        if(is.character(scn) | is.list(scn)){
            tmp=names(scn)
            if(length(tmp)==0 | sum(tmp=='SJ')!=1)
                scn=SCANVISmerge(scn,'mean',roi,gen)
            sj=scn$SJ
            if(length(grep('MUTS',names(scn)))>0){
                muts=cbind(rownames(scn$MUTS),apply(scn$MUTS,1,sum))
                if(!is.matrix(muts)) muts=t(as.matrix(muts))
            }
        }
    }
    sj=gsub(' ','',sj)

    ##keeping only SJ within the roi
    qs=intersect(which(as.numeric(sj[,'start'])>=as.numeric(roi[2])),
        which(as.numeric(sj[,'start'])<=as.numeric(roi[3])))
    qe=intersect(which(as.numeric(sj[,'end'])>=as.numeric(roi[2])),
        which(as.numeric(sj[,'end'])<=as.numeric(roi[3])))
    q=intersect(qs,qe)
    if(length(q)==0){
        stop('There are no SJs in the roi specified. 
        Consider expanding your genomic coordinates')
    }
    if(length(q)==0){
        print(paste('There are no splice junctions in the region',
            paste(roi.g,collapse=','),paste0(roi[1],':',roi[2],'-',roi[3])))
        print('Please consider expanding your region of interest')
        stop()
    }
    if(length(q)>0){
        sj=sj[q,]
        if(length(q)==1) sj=t(as.matrix(sj))
       ##expand XLIM to allow for any borderline Sjs
        XLIM=range(c(XLIM,as.numeric(as.vector(sj[,c('start','end')])))) 
        if(length(muts)==0){
            tmp=max(which(is.element(colnames(sj),c('FrameStatus','RRC'))))
            if(ncol(sj)>tmp){
                tmp2=unlist(strsplit(unlist(strsplit(sj[,(tmp+1):
                    ncol(sj)],':')),';'))
                tmp2=as.numeric(tmp2[(grep('chr',tmp2)+1)])
                XLIM=range(c(XLIM,tmp2)) #expand XLIM to inc all mapped varts
            }
        }
        if(length(muts)>0){
            tmp2=lapply(muts[,1],function(x) 
                unlist(strsplit(unlist(strsplit(unlist(strsplit(x,';'))[1],':'))
                    [2],'-')))
            tmp2=as.numeric(tmp2[(grep('chr',tmp2)+1)])
            XLIM=range(c(XLIM,tmp2))
        }
    }
    XLIM[1]=XLIM[1]-1
    XLIM[2]=XLIM[2]+1

    ############################################################################
    #START: plot exons
    IR.roi=IRanges(as.numeric(roi[2]),as.numeric(roi[3]))
    h=which(gen$EXONS[,'chr']==roi[1])
    tmp=gen$EXONS[h,]
    IR.gen=IRanges(as.numeric(tmp[,'start']),as.numeric(tmp[,'end']))
    v=as.matrix(findOverlaps(IR.roi,IR.gen))
    if(length(v)==0){
    	stop('The region of interest specified is an 
        intergenic space - no plot can be generated for these coordinates')
    }
    gen$EXONS=gen$EXONS[h[unique(v[,2])],]

    #making note of genes with >1 transcript
    tmp=unique(gen$EXONS[,c('gene_name','transcript')])
    tt=table(tmp[,1])
    TX=unique(gen$EXONS[which(is.element(gen$EXONS[,'gene_name'],
        names(tt)[which(tt>1)])),'transcript'])
    h=which(!duplicated(gen$EXONS[,'transcript']))
    gen$EXONS.utx=gen$EXONS[h,]

    #leaving space for legends
    XLIM=range(c(as.numeric(roi[2:3]),XLIM))
    d=diff(XLIM)
    XLIM[1]=XLIM[1]-(floor(0.02*d))
    XLIM[2]=XLIM[2]+(floor(0.1*d))
    
    #setting the plot
    G=unique(gen$EXONS[,'gene_name'])
    YLIM=c(-0.5,1.2)    
    if(length(bam)>0) YLIM[1]=YLIM[1]-0.2 #space for read coverage
    if(length(G)>1) YLIM[1]=YLIM[1]-0.2 #space for all genes
    if(full.annot) YLIM[1]=YLIM[1]-0.8 #space for all isoforms
    YLIM[1]=YLIM[1]-0.05 #just adding a bit more space
    if(length(TITLE)>0){
        plot(c(0,0),xlim=XLIM,ylim=YLIM,xaxt='n',yaxt='n',
            xlab='',ylab='',bty='n',type='l',col='white',main=TITLE)
    }
    if(length(TITLE)==0){
        plot(c(0,0),xlim=XLIM,ylim=YLIM,xaxt='n',yaxt='n',
            xlab='',ylab='',bty='n',type='l',col='white')
    }

    #plotting all exons 
    bh=0.008
    yy=(-0.05)
    u=unique(gen$EXONS[,c('start','end')])
    if(!is.matrix(u)) u=t(as.matrix(u))
    u=cbind(as.numeric(u[,1]),as.numeric(u[,2]))
    kul='darkblue'
    for(j in seq(1,nrow(u),1))
        rect(u[j,1],(yy+bh),u[j,2],(yy-bh),col=kul,border=kul)
    lines(range(u),c(yy,yy),lwd=1,col=kul)

    yy=-0.1    
    if(length(bam)>0) yy=-0.35

    #plot all GENES
    dyy=0.15/length(G)
    for(g in G){
        q=which(gen$EXONS[,'gene_name']==g)
        STRAND='->'
        if(unique(gen$EXONS[q,'strand'])=='-') STRAND='<-'
        u=unique(gen$EXONS[q,c('start','end')])
        if(!is.matrix(u)) u=t(as.matrix(u))
        u=cbind(as.numeric(u[,1]),as.numeric(u[,2]))
        if(length(G)>1){
        	#plotExonSet(u,yy,bh,'blue')
		    kul='blue'
		    for(j in seq(1,nrow(u),1))
		        rect(u[j,1],(yy+bh),u[j,2],(yy-bh),col=kul,border=kul)
		    lines(range(u),c(yy,yy),lwd=1,col=kul)
        }
        ulim=range(u)
        text(max(ulim)+(0.01*diff(ulim)),yy,paste0(STRAND,g),
            col='blue',cex=0.6,adj=0)
        yy=yy-dyy
    }

    ##for full annotation, plot every transcript
    if(full.annot){
        yy=yy-0.05
        dyy=0.8/length(TX)
        for(tx in TX){
            ##plot exons and a line to string all exons for one gene
            q=which(gen$EXONS[,'transcript']==tx)            
            u=unique(gen$EXONS[q,c('start','end')])
            if(!is.matrix(u)) u=t(as.matrix(u))
            u=cbind(as.numeric(u[,1]),as.numeric(u[,2]))
            #plotExonSet(u,yy,bh,'cornflowerblue')
            kul='cornflowerblue'
            for(j in seq(1,nrow(u),1))
                rect(u[j,1],(yy+bh),u[j,2],(yy-bh),col=kul,border=kul)
            lines(range(u),c(yy,yy),lwd=1,col=kul)

            ulim=range(u)
            id=unique(gen$EXONS[q,c('gene_name','transcript')])
            text((max(ulim)+(0.005*diff(ulim))),yy,id[1],
                col='cornflowerblue',cex=0.4,adj=0)
            text((max(ulim)-(0.005*diff(ulim))),yy,id[2],
                col='cornflowerblue',cex=0.4,adj=0)
            yy=yy-dyy
        }
    }
    roi=c(roi[1],range(as.numeric(as.vector(gen$EXONS[,c('start','end')]))))
    #roi=c(roi[1],range(c(as.numeric(gen$EXONS[,'start']),as.numeric(gen$EXONS[,'end']))))
    #END: plot exons
    ############################################################################

    ############################################################################
    #START: plot scn
    IRroi=IRanges(as.numeric(roi[2]),as.numeric(roi[3]))
    sj=sj[which(sj[,'chr']==roi[1]),]
    IRsj=IRanges(as.numeric(sj[,'start']),as.numeric(sj[,'end']))
    v=as.matrix(findOverlaps(IRsj,IRroi))
    sj=sj[unique(v[,1]),]
    #plot.scn(sj,SJ.special,USJ)
	scn=sj
    QQ=NULL
    if(length(SJ.special)>0){
        SJ.special=apply(SJ.special,1,function(x) paste(x,collapse=' '))
        QQ=which(is.element(apply(scn[,c('chr','start','end')],1,
            function(x) paste(x,collapse=' ')),SJ.special))
    }

    ##defining line widths according to RRS
    tmp=as.numeric(scn[,'RRS'])
    q0=which(tmp==min(tmp))[1]
    aa=abs(tmp-median(tmp))
    q1=which(aa==(min(aa)))[1]
    q2=which(tmp==max(tmp))[1]
    leg.lwd=cbind(round(as.numeric(scn[c(q0,q1,q2),'uniq.reads']),2),
        round(as.numeric(scn[c(q0,q1,q2),'RRS']),4))
    LWD=rep(0.25,length(tmp))
    zz=seq(0.05,1,0.05)
    zz=cbind(zz-0.05,zz)
    lwd0=0.25
    for(i in seq(1,nrow(zz),1)){
        q=intersect(which(tmp>zz[i,1]),which(tmp<=zz[i,2]))
        if(length(q)>0) LWD[q]=lwd0
        lwd0=lwd0+0.25
    }

    ##defining arc height according to uniq.reads
    AH=as.numeric(scn[,'uniq.reads'])
    AH=log2(AH)/max(log2(AH))
    AH=AH+0.1

    ##defining line type according to frame status
    LTY=rep(1,nrow(scn))
    q=grep('ENST',scn[,'FrameStatus']); if(length(q)>0) LTY[q]=2
    q=grep('OUT',scn[,'FrameStatus']); if(length(q)>0) LTY[q]=3

    q=grep('exon.skip',scn[,'JuncType'])
    if(length(q)>0) scn[q,'JuncType']='exon.skip'
    q=grep('alt',scn[,'JuncType'])
    if(length(q)>0) scn[q,'JuncType']='alt'
    ##setting junction colors ... blue for ASJs, variety for USJs
    kul=c('red','magenta','orange','mediumpurple2','green','azure3')
    names(kul)=c('exon.skip','alt','IsoSwitch','Unknown','NE','annot')

    ##defining plot limits
    y0=0
    NR=as.numeric()
    coor=cbind(as.numeric(scn[,'start']),as.numeric(scn[,'end']))
    coor=cbind(coor,abs(coor[,2]-coor[,1]))
    coor=cbind(coor,coor[,1]+(0.5*coor[,3]))

    scn[,'uniq.reads']=round(as.numeric(scn[,'uniq.reads']))
    I=grep('annot',scn[,'JuncType'],invert=TRUE)
    I=setdiff(I,grep('NE',scn[,'JuncType']))
    I=I[order(as.numeric(scn[I,'uniq.reads']),decreasing=TRUE)]
    Q=setdiff(c(grep('annot',scn[,'JuncType']),I),QQ)
    for(i in Q){
        if(scn[i,'JuncType']!='annot'){
            if(USJ=='NR'){
                text(coor[i,4],(AH[i]+0.03),scn[i,'uniq.reads'],
                col=kul[scn[i,'JuncType']],cex=0.3)
            }
            if(USJ=='RRS'){
                text(coor[i,4],(AH[i]+0.03),
                round(as.numeric(scn[i,'RRS']),2),
                col=kul[scn[i,'JuncType']],cex=0.3)
            }
        }
        draw.ellipse(coor[i,4],y0,coor[i,3]/2,AH[i],
            angle=0,segment=c(0,180),arc.only=TRUE,
            lwd=LWD[i],border=kul[scn[i,'JuncType']],lty=LTY[i])
    }
    if(length(QQ)>0){
        for(i in QQ){
            draw.ellipse(coor[i,4],y0,coor[i,3]/2,AH[i],
                angle=0,segment=c(0,180),arc.only=TRUE,
                lwd=LWD[i],border='cyan4',lty=LTY[i])
            if(USJ=='NR'){
                text(coor[i,4],(AH[i]+0.03),
                    scn[i,'uniq.reads'],col='cyan4',cex=0.3)
            }
            if(USJ=='RRS'){
                text(coor[i,4],(AH[i]+0.03),
                    round(as.numeric(scn[i,'RRS']),2),col='cyan4',cex=0.3)
            }
        }
    }

    tmp=range(coor[,c(1,2)])
    tmp2=0.05*diff(tmp)
    tmp=tmp+c(-tmp2,tmp2)
    lines(tmp,c(-0.002,-0.002),col='white',lwd=max(LWD)+2)

    for(i in grep('NE',scn[,'JuncType'])){
        jj=abs(jitter(0.01,100))
        lines(c(coor[i,1],coor[i,2]),c(-jj,-jj),col='green',lwd=1)
        ##lines(c(coor[i,1],coor[i,2]),c(-0.01,-0.01),col='green',lwd=0.5)
    }

    ##legend for USJ
    LEG=rep(1,6)
    LEG=cbind(LEG,c('exon.skip','alt5/3p',
        'IsoSwitch','Unknown','NovelExon','annot'))
    LEG=cbind(LEG,c('red','magenta','orange',
        'mediumpurple2','green','azure3'))
    legend('topright',LEG[,2],lwd=2,lty=as.numeric(LEG[,1]),
        col=LEG[,3],cex=0.4,bty='n')

    ##legend for LWD and LTY
    LEG=c('InFrame','UnknownFrame','OutFrame','',
        'Min.reads/RRS:','Med.reads/RRS:','Max.reads/RRS: ')
    LEG[5]=paste(LEG[5],paste0(leg.lwd[1,],collapse=','))
    LEG[6]=paste(LEG[6],paste0(leg.lwd[2,],collapse=','))
    LEG[7]=paste(LEG[7],paste0(leg.lwd[3,],collapse=','))
    LEG=cbind(c(1,2,3,1,1,1,1),LEG)
    LEG=cbind(LEG,c(1,1,1,1,LWD[c(q0,q1,q2)]))
    LEG=cbind(LEG,c(rep('grey',3),'white',rep('grey',3)))
    LEG=gsub(' ','',LEG)
    legend('topleft',LEG[c(1,2,3),2],lwd=LEG[c(1,2,3),3],
        lty=as.numeric(LEG[c(1,2,3),1]),col=LEG[c(1,2,3),4],cex=0.4,bty='n')
    legend('top',LEG[c(5,6,7),2],lwd=LEG[c(5,6,7),3],
        lty=as.numeric(LEG[c(5,6,7),1]),col=LEG[c(5,6,7),4],cex=0.4,bty='n')

    #END: plot scn
    ############################################################################

    ############################################################################
    #START: if bam url provided, plot read profile
    if(length(bam)!=0){
        #plot.bam(roi,bam,samtools)
		#plot.bam<-function(roi,bam,samtools){
	    bed=paste0(roi[1],':',roi[2],'-',roi[3])
	    write.table(bam,'tmpGRP.bam',sep='\t',quote=FALSE,
	        row.names=FALSE,col.names=FALSE)
	    #print('*** Getting read depth from bam files ... ***')
	    system(paste(samtools,'depth -a -r',bed,'-f tmpGRP.bam > tmpGRP'))
	    tmp=as.matrix(read.delim('tmpGRP',header=FALSE))
	    system('rm tmpGRP*')
	    dat=cbind(as.numeric(tmp[,2]),as.numeric(tmp[,3]))
	    rm(tmp)
	    n=nrow(dat)
	    if(n>20000){
	        # dat0=dat
	        # d=floor(n/20000)
	        # q=seq(1,n,d)
	        # z=tail(q,1)
	        # q=head(q,-1)
	        # dat=dat0[q,]
	        # for(j in seq(2,d,1))
	        #     dat[,2]=dat[,2]+dat0[q+j-1,2]
	        # dat[,2]=dat[,2]/d
	        # I=seq(z,nrow(dat0),1)
	        # if(length(I)>1) dat=rbind(dat,c(dat0[I[1],1],mean(dat0[I,2])))
	        # if(length(I)==1) dat=rbind(dat,dat0[I,])
	    }
	    dat[,2]=dat[,2]/max(dat[,2])
	    dat[,2]=dat[,2]*0.2
	    #lines(dat[,1],dat[,2]+0.1,col='gray18',lwd=0.5)
	    lines(dat[,1],-dat[,2]-0.1,col='gray18',lwd=0.5)
	    #print('*** DONE: Getting read depth from bam files ***')
	}
    #END: if bam url provided, plot read profile
    ############################################################################

    ############################################################################
    ##START: plot any mutations
    mut4plot=NULL
    if(length(muts)==0){
        q=grep('RRC',colnames(sj))
        if(length(q)==0) q=grep('FrameStatus',colnames(sj))
        n=ncol(sj)-q
        if(n>0){
            for(j in (q+1):ncol(sj)){
                h=which(sj[,j]!='')
                if(length(h)>0)
                    #plot_mut(sj[h,j])
                    mut4plot=c(mut4plot,sj[h,j])
            }
        }
        if(length(mut4plot)>1) mut4plot=cbind(mut4plot,1)
        if(length(mut4plot)==1) mut4plot=t(as.matrix(c(mut4plot,1)))
    }
    if(length(muts)>0){
        print('plotting mutations ...')
        #tmp=lapply(muts[,1],function(x) unlist(strsplit(unlist(strsplit(unlist(strsplit(x,';'))[1],':'))[2],'-')))
        #qq=which(unlist(lapply(tmp,function(x) sum(as.numeric(x)>=as.numeric(roi[2]))*sum(as.numeric(x)<=as.numeric(roi[3]))))>0)
        tmp=unlist(lapply(muts[,1],function(x) unlist(strsplit(x,';'))[1]))
        tmp=unlist(lapply(tmp,function(x) unlist(strsplit(x,':'))[2]))
        tmp=lapply(tmp,function(x) as.numeric(unlist(strsplit(x,'-'))))
        rn=as.numeric(roi[2:3])
        qq=which(unlist(lapply(tmp,function(x) sum(x>=rn[1])*sum(x<=rn[2])))>0)
        # if(length(qq)>1) plot_mut(muts[qq,c(1,ncol(muts))])
        # if(length(qq)==1) plot_mut(t(as.matrix(muts[qq,c(1,ncol(muts))])))
        if(length(qq)>1) mut4plot=rbind(mut4plot,muts[qq,c(1,ncol(muts))])
        if(length(qq)==1) mut4plot=rbind(mut4plot,muts[qq,c(1,ncol(muts))])
    }
    if(length(mut4plot)>0){
    	if(!is.matrix(mut4plot)) mut4plot=t(as.matrix(mut4plot))
		##   A typical mut4plot entry in col1 might look like this: 
		##   "chr1:953778 rs13303056 G->C|chr1:953779 rs13302945 A->C"
        ym0=0.005
        ymt0=-0.015
        if(!is.matrix(mut4plot))
                mut4plot=cbind(unique(unlist(strsplit(mut4plot,'\\|'))),1)
        N=as.numeric(mut4plot[,2])
        for(i in seq(1,nrow(mut4plot),1)){
	        d=tail(unlist(strsplit(mut4plot[i,1],';')),1)
	        tmp=unlist(strsplit(mut4plot[i,1],';'))[1]
	        tmp=unlist(strsplit(gsub('-',' ',gsub(':',' ',tmp)),' '))
	        tmp=c(tmp,d)
	        d=tail(tmp,1)
	        tmp=head(tmp,-1)
	        m=as.numeric(tmp[2:length(tmp)])
	        if(N[i]<2){
	                m=as.numeric(tmp[2:length(tmp)])
	                if(length(m)==1)
	                        points(m,ym0,col='black',pch=19,cex=0.5)
	                if(length(m)==2)
	                        lines(m,c(ym0,ym0),col='black',lwd=1)
	        }
	        if(N[i]>1){
	                ym0.tmp=ym0
	                dy=min(1/N[i],0.04)
	                cx=0.5
	                for(j in seq(1,N[i],1)){
	                        if(j==N[i]) cx=0.8
	                        if(length(m)==1) points(m,ym0.tmp,col='black',pch=19,cex=cx)
	                        if(length(m)>1) lines(m,c(ym0.tmp,ym0.tmp),col='black',lwd=0.7)
	                        #if(j==N[i]) text(m,ym0.tmp,N[i],col='white',cex=0.3)
	                        ym0.tmp=ym0.tmp+dy
	                }
	                if(length(m)==1){
	                        lines(c(m,m),c(ym0,ym0.tmp-dy),col='black')
	                        text(m,ym0.tmp-dy,N[i],col='white',cex=0.3)
	                }
	        }
	        if(length(m)>1) m=m[1]+(diff(m)/2)
	        text(m,ymt0,d,cex=0.4)                  
	    }
    }
    ##END: plot any mutations
    ############################################################################

    out=NULL
    out$roi=roi
    if(roi.g[1]!='') out$roi=c(out$roi,paste(roi.g,collapse=','))
    out$sj=sj
    if(length(muts)>0) out$muts=muts
    return(out)
}

## plotExonSet<-function(coor,yy,bh,kul='blue'){
##     for(j in seq(1,nrow(coor),1))
##         rect(coor[j,1],(yy+bh),coor[j,2],(yy-bh),col=kul,border=kul)
##     LIM=range(coor)
##     lines(LIM,c(yy,yy),lwd=1,col=kul)
## }

##get genomic coordinates for one or more gene names in same chr
gene2roi<-function(g,gen){
	GENg=gen$GENES[,'gene_name']
	if(length(intersect(g,GENg))==0)
	    GENg=gen$GENES[,'gene_id']
    q=which(is.element(GENg,g))
    if(length(q)==0){
        stop('Gene/s not found in the gencode object. 
        	Please check your gene list')
    }
    tmp=intersect(g,GENg[q])
    if(length(tmp)!=length(g)){
        print(paste('Gene/s',setdiff(g,tmp),'are not listed in the 
            supplied gencode object and will be excluded.'))
    }
    if(length(unique(gen$GENES[q,'chr']))!=1){
        print(paste('The genes supplied fall into the following 
            chromosomes:',paste(unique(gen$GENES[q,'chr']),collapse=',')))
        stop('The function gene2roi only accepts genes in one chromosome.
            Please adjust your input gene name/s')
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


