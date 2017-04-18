"zinv"<-function(y){
    y<-as.matrix(y)
    n<-nrow(y)
    yi<-y
    abc1<-.Fortran("inv1",as.double(y),as.integer(n),yi=as.double(yi))
    yi<-matrix(abc1$yi,nrow=n,byrow=T)
    list(yi=yi)
}



