\name{M-Quantile SAE}
\alias{mq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
M-Quantile SAE
}
\description{
A alternative method for small area estimation that robust for data with outlier and doesn't meet normal distribution.
}
\usage{
mq(y, x, regioncode.s, x.outs, regioncode.r, m, p, SE = 0, ydir = 0, tol.value = 1e-04, maxit.value = 100, k.value = 1.345)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  Response variable for observed unit.
}
  \item{x}{
  Auxiliary variable for observed unit.
}
  \item{regioncode.s}{
  Region code for observed unit.
}
  \item{x.outs}{
  Auxiliary variable for unobserved unit.
}
  \item{regioncode.r}{
  Region code for unobserved unit.
}
  \item{m}{
  Total area (area level)
}
  \item{p}{
  Size for x+1 (include intercept).
}
  \item{SE}{
  It is optional. If your data have a different level between observed unit and unobserved unit then you fill SE with standard error values for observed unit.
}
  \item{ydir}{
  It is optional. If your data have a different level between observed unit and unobserved unit then you fill ydir with response variable values for observed unit.
  }
  \item{tol.value}{
  Covergence tolerance limit for M-Quantile model. Default 0,0001.
}
  \item{maxit.value}{
  Maximum iteration for iterative reweighted least square. Default 100.
}
  \item{k.value}{
  Tuning constant for estimate Huber Proposal 2 scale. Default 1,345.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
hasilmq <- mq(data1$y, data1$x, data1$reg.s, data2$x, data2$reg.r, 30, 2)
}
}

