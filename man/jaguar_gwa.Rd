\name{jaguar_gwa}
\alias{jaguar_gwa}
\docType{package}
\title{
Perform genome-wide analysis
}
\description{
Computes p-value from our joint score test in a cis framework to map group-specific expression quantitative trait loci (eQTL) that tests for the shifts in gene expression patterns due to genotype and variability among tissues in a mixed effects model framework. A gene-level p-value is computed using a permutation-resampling scheme in order to investigate if a gene has at least one eQTL across all the groups. 
}
\usage{
jaguar_gwa(geneexp,genomat,write=FALSE)
}
\arguments{
        \item{geneexp}{ A matrix of gene expression data with gene on rows and patient samples on columns. There has to be equal number of samples in each group. Samples (columns in the geneexp object) with missing gene expression values for any group/tissue MUST be included in the data}
        \item{genomat}{A matrix of genotype data recoded as single allele dosage number (i.e. 0, 1 or 2) with rows representing SNPs and columns representing samples}
        \item{write}{Boolean value indicating whether the results should be outputted into a tab-delimited text file. The default value is FALSE.}
}
\value{
        \item{}{A matrix of raw unadjusted p-values with rows representing genes and columns representing SNPs}
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
  \code{\link{jaguar_process},\link{jaguar_slice},\link{jaguar_sim},\link{jaguar_cis},\link{jaguar_plotqtl}}
}
\keyword{score test statistic }
\examples{
# Load the example data
data(jaguar_example);

# Run a cis analysis with no permutations
Gene_Mat = as.matrix(jaguar_example$GENE_EXP[1:10,])
Geno_Mat = as.matrix(jaguar_example$GENO_MAT)
out = jaguar_gwa(Gene_Mat,Geno_Mat);
dim(out)

}
\keyword{joint score test}
\keyword{genotype}
\keyword{eQTL}
\keyword{interaction}
