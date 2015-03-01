\name{jaguar}
\alias{jaguar}
\docType{package}
\title{
Runs the joint score test statistic approach
}
\description{
Computes joint score test statistic to map group-specific expression quantitative trait loci  (eQTL) that tests for the shifts in gene expression patterns due to genotype and variability among tissues in a mixed effects model framework. 
}
\usage{
jaguar(geneexp,geno,ngroups,write=FALSE)
}
\arguments{
        \item{geneexp}{ A matrix of gene expression data with gene on rows and patient samples on columns. Missing values not allowed. There has to be equal number of samples in each group.}
        \item{geno}{A matrix of genotype data recoded as single allele dosage number (i.e. 0, 1 or 2) with rows representing SNPs and columns representing samples }
        \item{ngroups}{An integer representing the number of groups in the data}
        \item{write}{Boolean value indicating whether the results should be outputted into a tab-delimited text file}
}
\value{
        \item{}{A matrix of raw unadjusted p-values with rows representing genes and columns representing SNPs}
}

\author{
Chaitanya R. Acharya, Andrew S. Allen
Maintainer: Chaitanya Acharya<c.acharya@duke.edu>
}
\references{
Chaitanya R. Acharya, Kouros Owzar and Andrew S. Allen; Joint analysis of genotype and tissue-specific variability using a novel score test to map eQTL (Manuscript submitted)
}
\seealso{
  \code{\link{calcThreshold},\link{lin},\link{ProcessJaguarResults},\link{SliceGeneData},\link{jaguarSIM}}
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
out = jaguar(GeneExp,Geno,k,write=FALSE);
dim(out);
}
\keyword{score test}
\keyword{genotype}
\keyword{eQTL}
\keyword{interaction}