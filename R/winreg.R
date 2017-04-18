winreg<-function(y1,y2,d1,d2,z){
  z<-as.matrix(z)
  nc<-ncol(z)
  fit2<-coxph(Surv(y2,d2) ~ z)
  rfit2<-matrix(residuals(fit2,type="score"),byrow=F,ncol=nc)
  fit1<-coxph(Surv(y1,d1) ~ z)
  rfit1<-matrix(residuals(fit1,type="score"),byrow=F,ncol=nc)
  sigma2<-solve(t(rfit2)%*%rfit2)
  sigma1<-solve(t(rfit1)%*%rfit1)
  beta1<-matrix(fit1$coef,byrow=F,ncol=1)
  beta2<-matrix(fit2$coef,byrow=F,ncol=1)
  beta<-beta1+beta2
  sigmax<-matrix(0,ncol=2*nc,nrow=2*nc)
  sigmax[1:nc,1:nc]<-sigma1
  sigmax[(nc+1):(2*nc),(nc+1):(2*nc)]<-sigma2
  sigmax[1:nc,(nc+1):(2*nc)]<-sigma1%*%t(rfit1)%*%rfit2%*%sigma2
  sigmax[(nc+1):(2*nc),1:nc]<-sigma2%*%t(rfit2)%*%rfit1%*%sigma1
  dx<-cbind(diag(nc),diag(nc))
  sigma<-dx%*%sigmax%*%t(dx)
  tb1<-beta1/sqrt(diag(sigma1));pb1<-2*(1-pnorm(abs(tb1)))
  tb2<-beta2/sqrt(diag(sigma2));pb2<-2*(1-pnorm(abs(tb2)))
  tb<-beta/sqrt(diag(sigma));pb<-2*(1-pnorm(abs(tb)))
  list(beta1=beta1,sigma1=matrix(sigma1,byrow=F,ncol=nc),tb1=tb1,pb1=pb1,
       beta2=beta2,sigma2=matrix(sigma2,byrow=F,ncol=nc),tb2=tb2,pb2=pb2,
       beta=beta,  sigma=matrix(sigma,byrow=F,ncol=nc),tb=tb,pb=pb)
}



