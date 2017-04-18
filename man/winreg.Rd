\name{winreg}
\alias{winreg}
\title{Double Cox regression for win product}
\description{Use two Cox regression models (one for the terminal event and the other for the non-trminal event) to model the win product adjusting for covariates}
\usage{
winreg(y1,y2,d1,d2,z)
}
\arguments{
\item{y1}{ a numeric vector of event times denoting the minimum of event times \eqn{T_1}, \eqn{T_2} and censoring time \eqn{C}, where the endpoint \eqn{T_2}, corresponding to the terminal event, is considered of higher clinical importance than the endpoint \eqn{T_1}, corresponding to the non-terminal event. Note that the terminal event may censor the non-terminal event, resulting in informative censoring.}
\item{y2}{ a numeric vector of event times denoting the minimum of event time \eqn{T_2} and censoring time \eqn{C}. Clearly, y2 is not smaller than y1. }
\item{d1}{ a numeric vector of event indicators with 1 denoting the non-terminal event is observed and 0 else. }
\item{d2}{ a numeric vector of event indicators with 1 denoting the terminal event is observed and 0 else. }
\item{z}{ a numeric matrix of covariates.}
}
\details{This function uses two Cox regression models (one for the terminal event and the other for the non-trminal event) to model the win product adjusting for covariates.}
\value{
\item{beta1}{Estimated regression parameter based on the non-terminal event times \code{y1}, \eqn{\exp}(\code{beta1}) is the adjusted hazard ratio}
\item{sigma1}{Estimated variance of \code{beta1} using the residual method instead of the inverse of Fisher information}
\item{tb1}{Wald test statistics based on \code{beta1} and \code{sigma1}}
\item{pb1}{Two-sided p-values of the Wald test statistics \code{tb1}}
\item{beta2}{Estimated regression parameter based on the terminal event times \code{y2}, \eqn{\exp}(\code{beta2}) is the adjusted hazard ratio}
\item{sigma2}{Estimated variance of \code{beta2} using the residual method instead of the inverse of Fisher information}
\item{tb2}{Wald test statistics based on \code{beta2} and \code{sigma2}}
\item{pb2}{Two-sided p-values of the Wald test statistics \code{tb2}}
\item{beta}{\code{beta1}+\code{beta2},\eqn{\exp}(-\code{beta}) is the adjusted win product}
\item{sigma}{Estimated variance of \code{beta} using the residual method}
\item{tb}{Wald test statistics based on \code{beta} and \code{sigma}}
\item{pb}{Two-sided p-values of the Wald test statistics \code{tb}}
}
\references{
Pocock S.J., Ariti C.A., Collier T. J. and Wang D. 2012. The win ratio: a new approach to the analysis of composite endpoints in clinical trials based on clinical priorities. European Heart Journal, 33, 176-182.

Luo X., Tian H., Mohanty S. and Tsai W.-Y. 2015. An alternative approach to confidence interval estimation for the win ratio statistic. Biometrics, 71, 139-145.

Luo X., Qiu J., Bai S. and Tian H. 2017. Weighted win loss approach for analyzing prioritized outcomes. Statistics in Medicine, to appear.
}
\author{Xiaodong Luo}
\seealso{\code{\link{wrlogistic}}}
\examples{
###Generate data
n<-300
rho<-0.5
b2<-c(1.0,-1.0)
b1<-c(0.5,-0.9)
bc<-c(1.0,0.5)
lambda10<-0.1;lambda20<-0.08;lambdac0<-0.09
lam1<-rep(0,n);lam2<-rep(0,n);lamc<-rep(0,n)
z1<-rep(0,n)
z1[1:(n/2)]<-1
z2<-rnorm(n)
z<-cbind(z1,z2)

lam1<-lam2<-lamc<-rep(0,n)
for (i in 1:n){
    lam1[i]<-lambda10*exp(-sum(z[i,]*b1))
    lam2[i]<-lambda20*exp(-sum(z[i,]*b2))
    lamc[i]<-lambdac0*exp(-sum(z[i,]*bc))
}
tem<-matrix(0,ncol=3,nrow=n)

y2y<-matrix(0,nrow=n,ncol=3)
y2y[,1]<-rnorm(n);y2y[,3]<-rnorm(n)
y2y[,2]<-rho*y2y[,1]+sqrt(1-rho^2)*y2y[,3]
tem[,1]<--log(1-pnorm(y2y[,1]))/lam1
tem[,2]<--log(1-pnorm(y2y[,2]))/lam2
tem[,3]<--log(1-runif(n))/lamc

y1<-apply(tem,1,min)
y2<-apply(tem[,2:3],1,min)
d1<-as.numeric(tem[,1]<=y1)
d2<-as.numeric(tem[,2]<=y2)

y<-cbind(y1,y2,d1,d2)
z<-as.matrix(z)
aa<-winreg(y1,y2,d1,d2,z)
aa
}
\keyword{win ratio}
\keyword{win loss statistics}
\keyword{win product}
\keyword{Cox regression}
\keyword{Win loss regression}






