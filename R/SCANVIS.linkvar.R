##Phaedra Agius, April 2019
##
##This script maps variants (input in bed format) to SJs (as output by 
##SCANVIS.scan function) by overlapping variant coordinates with the 
##coordinates of the gene/s that overlap the SJ
##A column is tagged onto the SJs matrix indicating any present variants. 
##Script may be run any number of times for multiple variant sets 
##(eg. variants from different callers)
##
#INPUT: scn=output of SCANVIS.scan
##   bed=variants in bed format (chr,start,end plus a single description 
##       column which can contain REF->ALT and any other details in text)
##       NOTE: fourth column of bed file should have a descriptive name 
##           eg. ssSNP for splice site SNPs
##           These ids will be useful when visualizing with SJs
##           Bed input for structural variants:
##           Eg. A str variant chr1:10001-10010 chr2:2100-2500 rs12345 
##               would be represented in the bed file as
##               chr1 10001 10010 rs12345
##               chr2 2100 2500 rs12345
##   gen=gencode data as output by SCANVIS.annotation, 
##       user must supply this if variants are to be mapped
##       to gene coordinates that overlap SJs. If not supplied, variants are 
##       mapped to the genomic intervals defined by SCANVIS for RRS scoring 
##   p=padding bp length where SJ intervals get extended just a little further 
##     to include nearby variants (default: p=0, no padding)
##OUTPUT: scn as input with one additional column indicating any variants 
##       present in the genes that overlap the SJ

SCANVIS.linkvar<-function(scn,bed,gen,p=0){

    scn=gsub(' ','',scn)
    bed=gsub(' ','',bed)

    if(ncol(bed)>4){
        bed[,4]=unlist(lapply(bed[,4:ncol(bed)],
            function(x) paste(x,collapse=',')))
        bed=bed[,seq(1,4,1)]
    }
    
    if(!is.matrix(scn)) scn=t(as.matrix(scn))
    out=rep('',nrow(scn))

    CHR=intersect(bed[,'chr'],scn[,'chr'])
    if(length(CHR)==0)
        stop('chr formats in scn & bed do not match ... no communal chrs found')
    
    SNP=is.logical(all.equal(bed[,3],bed[,2]))
    GENES=gen$GENES

    for(chr in rev(CHR)){
        print(paste('Mapping variants to SJs for',chr))
        B=which(bed[,'chr']==chr)
        bed.tmp=unique(bed[B,])
        if(length(B)==1) bed.tmp=t(as.matrix(bed[B,]))
        bed.tmp.w=as.numeric(bed.tmp[,3])-as.numeric(bed.tmp[,2])
        bed.id=paste0(bed.tmp[,1],':',bed.tmp[,2],'-',
            bed.tmp[,3],';',bed.tmp[,4])
        h=which(bed.tmp.w==1)
        if(length(h)>0)
            bed.id[h]=paste0(bed.tmp[h,1],':',bed.tmp[h,2],';',bed.tmp[h,4])
        IRmut=IRanges(as.numeric(bed.tmp[,'start']),as.numeric(bed.tmp[,'end']))

        Q=which(GENES[,'chr']==chr)
        IRgen=IRanges(as.numeric(GENES[Q,'start']),as.numeric(GENES[Q,'end']))
        v1=as.matrix(findOverlaps(IRmut,IRgen))
        if(p>0){
            IRmut2=IRanges(as.numeric(bed.tmp[,'start'])-p,
                as.numeric(bed.tmp[,'end'])+p)
            v1=as.matrix(findOverlaps(IRmut2,IRgen))
        }

        q=which(scn[,'chr']==chr)
        IRscn=IRanges(as.numeric(scn[q,'start']),as.numeric(scn[q,'end']))
        v2=as.matrix(findOverlaps(IRscn,IRgen))

        xx=intersect(v1[,2],v2[,2])
        if(length(xx)>0){
            v2=v2[which(is.element(v2[,2],xx)),]
            if(length(v2)>0){
                if(!is.matrix(v2)) v2=t(as.matrix(v2))            
                u=unique(v2[,1])
                tmp=unlist(lapply(u,
                    function(x) paste(unique(bed.id[v1[which(is.element(v1[,2],
                        v2[which(v2[,1]==x),2])),1]]),collapse='|')))
                out[q[u]]=tmp
            }
        }        
    }
    scn=cbind(scn,out)
    colnames(scn)[ncol(scn)]=colnames(bed)[4]
    return(scn)
}







