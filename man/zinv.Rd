\name{zinv}
\alias{zinv}
\title{Inverse matrix}
\description{This will calculate the inverse matrix by Gauss elimination method}
\usage{
zinv(y)
}
\arguments{
  \item{y}{a sqaure matrix}
}
\details{Inverse matrix}
\value{\item{yi}{the inverse of \code{y}}
}
\author{Xiaodong Luo}
\note{This provides the inverse matrix using Gauss elimination method, this program performs satisfactorily when the size of the matrix is less than 50}
\examples{
y<-matrix(c(1,2,0,1),ncol=2,nrow=2)
zinv(y)
}
