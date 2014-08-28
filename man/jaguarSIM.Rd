\name{jaguarSIM}
\alias{jaguarSIM}
\docType{package}
\title{
Run null or power simulations
}
\description{
Function to run power/null simulations by simulating one gene and one SNP at a time. The objective of these simulations is two pronged - 1) Check for the type I error control, and 2) Compare two different null hypotheses where one's called a global null and other is local null. Under the global null hypotheses, we fit a linear mixed-effects model where we assume that there is no main genotypic effect and group-specific variability in the data. This leads to a 2 degree-of-freedom joint score test. Under the local null, we fit a linear mixed-effects model where we assume only the absence of group-specific variability. This leads to a 1 degree-of-freedom variance component score test.
}
\usage{
jaguarSIM(nobs = 500, k = 4, tau = 1, eps = 1,gamma = 0,bta = 0,maf = 0.10)
}
\arguments{
        \item{nobs}{The number of observations in each group}
        \item{k}{The total number of groups}
        \item{tau}{Variance component of the subject-specific random effect}
        \item{eps}{Variance component of the residual error}
        \item{gamma}{Variance component of group-specific interaction with the genotype as a random effect}
        \item{bta}{Main genotypic effect as a fixed-effect}
        \item{maf}{Minor allele frequency}
}
\value{
        \item{}{A numeric vector consisting of four different p-values. U0_bta and U0_gamma are the individual components of our joint score test statistic, U_psi under the global null. U_gamma shows the variance component score test result under the local null.} 
}
\details{
This function currently implements only balanced designs with equal number of observations in each group. A linear mixed effects model as described in the manuscript is fit to run the analyses. Please refer to our manuscript for more details. 
}
\author{
Chaitanya R. Acharya, Andrew S. Allen
Maintainer: Chaitanya Acharya<c.acharya@duke.edu>
}
\references{
Chaitanya R. Acharya, Kouros Owzar and Andrew S. Allen; Joint analysis of genotype and tissue-specific variability using a novel score test to map eQTL (Manuscript submitted)
}
\seealso{
  \code{\link{calcThreshold},\link{ComputeLinMax},\link{ProcessJaguarResults},\link{SliceGeneData},\link{jaguar}}
}
\examples{
## An example to perform some null simulations
## NOTE: 10 sims are definitely not enough. Please try between 1000-10000.

nsim=10; alpha=0.05;
test = do.call("rbind",rlply(nsim,.progress="time",jaguarSIM(nobs=10,k=4)));
null.sim = apply(test,2,function(x) sum(x<=alpha)/nsim);
}
\keyword{score test}
\keyword{simulations}
\keyword{GWAS}
\keyword{eQTL}