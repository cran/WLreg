\name{wrlogistic}
\alias{wrlogistic}
\title{Logistic regression for win ratio}
\description{Use a logistic regression model to model win ratio adjusting for covariates with the user-supplied comparison results}
\usage{
wrlogistic(aindex,z,b0=rep(0,ncol(z)),tol=1.0e-04,maxiter=20)
}
\arguments{
\item{aindex}{a vector that collects the pairwise comparison results. Suppose there are a total of \eqn{n} subjects in the study, there are \eqn{n(n-1)/2} elements in \code{aindex}. The \eqn{(i-1)*(i-2)/2+j}-th element, denoted by \eqn{C_{ij}}, is the comparison result between subject \eqn{i} and subject \eqn{j}, where \eqn{i=2,\ldots,n} and \eqn{j=1,\ldots,i-1}. The element \eqn{C_{ij}} is equal to \eqn{1} if subject \eqn{i} wins over subject \eqn{j} on the most important outcome, \eqn{C_{ij}} is equal to \eqn{-1} if subject \eqn{i} loses against subject \eqn{j} on the most important outcome; \eqn{C_{ij}} is equal to \eqn{2} if subject \eqn{i} wins over subject \eqn{j} on the second most important outcome after tie on the most important outcome, \eqn{C_{ij}} is equal to \eqn{-2} if subject \eqn{i} loses against subject \eqn{j} on the second most important outcome after tie on the most important outcome; and so forth until all the outcomes have been used for comparison; then \eqn{C_{ij}} is equal to \eqn{0} if an ultimate tie is resulted.}
\item{z}{a matrix of covariates}
\item{b0}{the initial value of the regression parameter}
\item{tol}{error tolerence}
\item{maxiter}{maximum number of iterations}
}
\details{This function uses a logistic regression model to model win ratio adjusting for covaraites. This function uses the pairwise comparision result supplied by the user which hopefully will speed up the program.}
\value{
\item{b}{Estimated regression parameter, \eqn{\exp}(\code{b}) is the adjusted win ratio}
\item{Ubeta}{The score function}
\item{Vbeta}{The estimated varaince of \eqn{\sqrt{n}\times}\code{b}}
\item{Wald}{Wald test statistics for the estimated parameter \code{b}}
\item{pvalue}{Two-sided p-values of the Wald statistics}
\item{Imatrix}{The information matrix}
\item{Wtotal}{Total wins}
\item{Ltotal}{Total losses}
\item{err}{err at convergence}
\item{iter}{number of iterations performed before covergence}
}
\references{
Pocock S.J., Ariti C.A., Collier T. J. and Wang D. 2012. The win ratio: a new approach to the analysis of composite endpoints in clinical trials based on clinical priorities. European Heart Journal, 33, 176-182.

Luo X., Tian H., Mohanty S. and Tsai W.-Y. 2015. An alternative approach to confidence interval estimation for the win ratio statistic. Biometrics, 71, 139-145.

Luo X., Qiu J., Bai S. and Tian H. 2017. Weighted win loss approach for analyzing prioritized outcomes. Statistics in Medicine, to appear.
}
\author{Xiaodong Luo}
\seealso{\code{\link{winreg}}}
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
###################

#####Define the comparison function
comp<-function(y,x){
  y1i<-y[1];y2i<-y[2];d1i<-y[3];d2i<-y[4]
  y1j<-x[1];y2j<-x[2];d1j<-x[3];d2j<-x[4] 
  w2<-0;w1<-0;l2<-0;l1<-0
  
  if (d2j==1 & y2i>=y2j) w2<-1
  else if (d2i==1 & y2j>=y2i) l2<-1
  
  if (w2==0 & l2==0 & d1j==1 & y1i>=y1j) w1<-1
  else if (w2==0 & l2==0 & d1i==1 & y1j>=y1i) l1<-1
  
  comp<-0
  if (w2==1) comp<-1 
  else if (l2==1) comp<-(-1)
  else if (w1==1) comp<-2
  else if (l1==1) comp<-(-2)
  
  comp
}
bin<-rep(0,n*(n-1)/2)
for (i in 2:n)for (j in 1:(i-1))bin[(i-1)*(i-2)/2+j]<-comp(y[i,],y[j,])
###Use the win loss indicator matrix to calculate the general win loss statistics
bb2<-wrlogistic(bin,z,b0=rep(0,ncol(z)),tol=1.0e-04,maxiter=20)
bb2


####Calculate the win, loss, tie result using Fortran loops to speed up the process
####Using the "inline" package to convert the code into Fortran

#install.packages("inline") #Install the package "inline''
library("inline") ###Load the package "inline"

########################################################
# The use of ``inline'' needs ``rtools'' and ``gcc''
# in the PATH environment of R.
# The following code will put these two into
# the PATH for the current R session ONLY.
########################################################

#rtools <- "C:\\Rtools\\bin"
#gcc <- "C:\\Rtools\\gcc-4.6.3\\bin"
#path <- strsplit(Sys.getenv("PATH"), ";")[[1]]
#new_path <- c(rtools, gcc, path)
#new_path <- new_path[!duplicated(tolower(new_path))]
#Sys.setenv(PATH = paste(new_path, collapse = ";"))

codex4 <- "
integer::i,j,indexij,d1i,d2i,d1j,d2j,w2,w1,l2,l1
double precision::y1i,y2i,y1j,y2j
do i=2,n,1
   y1i=y(i,1);y2i=y(i,2);d1i=dnint(y(i,3));d2i=dnint(y(i,4))
   do j=1,(i-1),1
      y1j=y(j,1);y2j=y(j,2);d1j=dnint(y(j,3));d2j=dnint(y(j,4)) 
      indexij=(i-1)*(i-2)/2+j
      w2=0;w1=0;l2=0;l1=0
      if (d2j==1 .and. y2i>=y2j) then
         w2=1
      else if (d2i==1 .and. y2j>=y2i) then
         l2=1
      else if (d1j==1 .and. y1i>=y1j) then
         w1=1
      else if (d1i==1 .and. y1j>=y1i) then
         l1=1
      end if
      aindex(indexij)=0
      if (w2==1) then 
         aindex(indexij)=1
      else if (l2==1) then
         aindex(indexij)=-1
      else if (w2==0 .and. l2==0 .and. w1==1) then
         aindex(indexij)=2
      else if (w2==0 .and. l2==0 .and. l1==1) then
         aindex(indexij)=-2
      end if
   end do
end do
"
###Convert the above code into Fortran
cubefnx4<-cfunction(sig = signature(n="integer", p="integer", y="numeric", aindex="integer"), 
            implicit = "none",dim = c("", "", "(n,p)", "(n*(n-1)/2)"), codex4, language="F95")

###Use the converted code to calculate the win, loss and tie indicators
options(object.size=1.0E+10)
ain<-cubefnx4(length(y[,1]),length(y[1,]), y, rep(0,n*(n-1)/2))$aindex

####Perform the logistic regression
aa2<-wrlogistic(ain,z,b0=rep(0,ncol(z)),tol=1.0e-04,maxiter=20)
aa2

}
\keyword{weighted}
\keyword{win ratio}
\keyword{pairwise comparison}
\keyword{logistic regression}
\keyword{conditional logistic regression}





