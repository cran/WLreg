"wrlogistic"<-function(aindex,z,b0=rep(0,ncol(z)),tol=1.0e-04,maxiter=20){
    z<-as.matrix(z)
    n<-nrow(z)
    p<-ncol(z)
    b1<-Ubeta<-Wtotal<-Ltotal<-matrix(b0,nrow=p,ncol=1)
    VU<-Vbeta<-Imatrix<-matrix(0,nrow=p,ncol=p)
    err<-1.0
    iter<-1
    abc<-.Fortran("wrlogistica",as.integer(n),as.integer(p),as.double(z),as.integer(aindex),as.double(b0),
                                as.double(tol), as.integer(maxiter),b1=as.double(b1),Ubeta=as.double(Ubeta),
                                VU=as.double(VU),Vbeta=as.double(Vbeta),Imatrix=as.double(Imatrix),
                                Wtotal=as.double(Wtotal),Ltotal=as.double(Ltotal),err=as.double(err),iter=as.integer(iter))
    VU<-matrix(abc$VU,ncol=p,byrow=T)
    Vbeta<-matrix(abc$Vbeta,ncol=p,byrow=T)
    Imatrix<-matrix(abc$Imatrix,ncol=p,byrow=T)
    Waldtest<-abc$b1/sqrt(diag(Vbeta)/n)
    pWaldtest<-(1-pnorm(abs(Waldtest)))*2
    list(b=abc$b1,Ubeta=abc$Ubeta,VU=VU,Vbeta=Vbeta,Wald=Waldtest,pvalue=pWaldtest,Imatrix=Imatrix,Wtotal=abc$Wtotal,Ltotal=abc$Ltotal,
         err=abc$err,iter=abc$iter)
}


