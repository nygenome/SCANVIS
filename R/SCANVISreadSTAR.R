##INPUT: url to star sj file
SCANVISread_STAR<-function(sj_file){
    out=as.matrix(read.delim(sj_file,header=FALSE))
    colnames(out)=c("chr","start","end","strand","intron.motif","is.annotated",
        "uniq.reads","multi.reads","max.overhang")
    out=gsub(' ','',out)
    return(out)
}
