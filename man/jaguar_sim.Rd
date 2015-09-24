\name{jaguar_sim}
\alias{jaguar_sim}
\docType{package}
\title{
Run null or power simulations
}
\description{
Function to run power/null simulations by simulating one gene and one SNP at a time. The objective of these simulations is two pronged - 1) Check for the type I error control for the joint score test statistic, and 2) Compare two different null hypotheses where one's called a global null (bta=0 and PVEg =0) and other is local null (PVEg=0). Under the global null hypotheses, we fit a model where we assume that there is no main genotypic effect and group-specific variability in the data. Under the local null, we fit a model where we assume only the absence of group-specific variability. This is essentially a variance component score test.
}
\usage{
jaguar_sim(nobs = 500, k = 5, tau = 1, eps = 1, PVEg = 0, bta = 0, maf = 0.10)
}
\arguments{
        \item{nobs}{The number of observations in each group}
        \item{k}{The total number of groups}
        \item{tau}{Variance component of the subject-specific random effect}
        \item{eps}{Variance component of the residual error}
        \item{PVEg}{Proportion of variance explained by gamma}
        \item{bta}{Additive genotypic effect as a fixed-effect}
        \item{maf}{Minor allele frequency}
}
\value{
        \item{}{A numeric vector consisting of two different p-values, "VCScoreTest" and "JointScoreTest" with the former indicating the p-value from the variance component score test and the latter indicating the p-value from the joint score test.} 
}
\details{
This function currently implements only balanced designs with equal number of observations in each group. For each individual, we model the potential genetic association between a target SNP and the expression level a target gene (in multiple tissues) at a single locus using the following mixed effects model (i = individual; t = tissue) - 
\deqn{
y_{i,t}=\alpha_t + g_i \beta_i + b_t g_i + u_i + \epsilon_{i,t}
}
where \eqn{y_{i,t}} is a t-dimensional vector of gene expression data for individual i, \eqn{g_i} is the scalar value of genotype in allele-dosage format, \eqn{b_t} is a t-dimensional tissue-specific random effect where \eqn{b \sim N(0,\gamma)}, \eqn{u_i} is the scalar value representing individual-specific random effect where \eqn{u \sim(0,\tau)}.
}

\author{
Chaitanya R. Acharya, Andrew S. Allen
Maintainer: Chaitanya Acharya<c.acharya@duke.edu>
}
\references{
Chaitanya R. Acharya, Kouros Owzar, Janice M. McCarthy and Andrew S. Allen; Exploiting expression patterns across multiple tissues to map expression quantitative trait loci (Manuscript submitted)

Chaitanya R. Acharya and Andrew S. Allen; JAGUAR: An R-package to Implement Joint Analysis of Genotype and Group-Specific Variability Using a Novel Score Test to Map eQTL (Manuscript submitted)
}
\seealso{
  \code{\link{jaguar_process},\link{jaguar_slice},\link{jaguar_gwa},\link{jaguar_cis},\link{jaguar_plotqtl}}
}
\examples{
## An example to perform some null simulations
## NOTE: 10 sims are obviously not enough. Please try between 1000-10000.

nsim=10; alpha=0.05;
test = do.call("rbind",rlply(nsim,.progress="none",jaguar_sim(nobs=100,k=4)));
null.sim = apply(test,2,function(x) sum(x<=alpha)/nsim);
}
\keyword{score test}
\keyword{simulations}
\keyword{GWAS}
\keyword{eQTL}
