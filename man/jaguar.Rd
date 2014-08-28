\name{jaguar}
\alias{jaguar}
\docType{package}
\title{
Computes joint score test statistic
}
\description{
Computes joint score test statistic to map group-specific expression quantitative trait loci  (eQTL) that tests for the shifts in gene expression patterns due to genotype and variability among tissues in a mixed effects model framework. 
}
\usage{
jaguar(geneexp,geno,ngroups)
}
\arguments{
        \item{geneexp}{ A matrix of gene expression data with gene on rows and patient samples on columns. Missing values not allowed. There has to be equal number of samples in each group.}
        \item{geno}{A matrix of genotype data recoded as single allele dosage number (i.e. 0, 1 or 2) with rows representing SNPs and columns representing samples }
        \item{ngroups}{An integer representing the number of groups in the data}
}
\value{
        \item{}{A matrix of score test statistics with rows representing genes and columns representing SNPs}
}
\details{
\tabular{ll}{
Package: \tab JAGUAR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-08-28\cr
License: \tab GPL-3\cr
}
}
\author{
Chaitanya R. Acharya, Andrew S. Allen
Maintainer: Chaitanya Acharya<c.acharya@duke.edu>
}
\references{
Chaitanya R. Acharya, Kouros Owzar and Andrew S. Allen; Joint analysis of genotype and tissue-specific variability using a novel score test to map eQTL (Manuscript submitted)
}
\seealso{
  \code{\link{calcThreshold},\link{ComputeLinMax},\link{ProcessJaguarResults},\link{SliceGeneData},\link{jaguarSIM}}
}
\keyword{score test statistic }
\examples{
# Load the example data
data(example.data);

# Set the parameters
GeneExp = as.matrix(example.data$Gene); 
Geno = as.matrix(example.data$Geno); 
k = example.data$ngroups;

# Run JAGUAR
out = jaguar(GeneExp,Geno,k);
dim(out);
}
\keyword{score test}
\keyword{genotype}
\keyword{eQTL}
\keyword{interaction}