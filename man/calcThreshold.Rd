\name{calcThreshold}
\alias{calcThreshold}
\docType{package}
\title{
Calculates various threshold values to control for the family-wise error rate (FWER).
}
\description{
Function to calculate thresholds values. Only Bonferroni and Danyu Lin's efficient Monte Carlo method are allowed.
}
\usage{
calcThreshold(nsnps,ngenes,method=c("bonferroni","lin"),path,alpha=0.05)
}
\arguments{
        \item{nsnps}{An integer indicating the total number of SNPs in the study}
        \item{ngenes}{An integer indicating the total number of genes in the study}
        \item{path}{Path to all the directories containing mcMIN.txt files. This is essential when gene expression data is sliced into multiple fragments.}
        \item{method}{Choose method to adjust the FWER. Only Bonferroni and Lin's method are provided for now.}
        \item{alpha}{The nominal p-value indicating the false positive rate}
}
\value{
        \item{}{Numeric threshold value}
}
\author{
Chaitanya R. Acharya, Andrew S. Allen
Maintainer: Chaitanya Acharya<c.acharya@duke.edu>
}
\references{
An efficient Monte Carlo approach to assessing the statistical significance in genomic studies. Lin, D.Y. Bioinformatics. 21(6) 2005. 

Joint analysis of genotype and tissue-specific variability using a novel score test to map eQTL. Chaitanya R. Acharya, Kouros Owzar and Andrew S. Allen (Manuscript submitted)
}
\seealso{
  \code{\link{jaguar},\link{ProcessJaguarResults},\link{SliceGeneData},\link{jaguar}}
}
\examples{
## Example
# 
# Load the example data (not an ideal dataset)
data(example.data);

# Set the parameter values
GeneExp = as.matrix(example.data$Gene); 
GenoMat = as.matrix(example.data$Geno);
ngenes = nrow(GeneExp)
nsnps = nrow(GenoMat)

# Calculate the threshold values
calcThreshold(ngenes,nsnps,method="bonferroni",alpha=0.05);

}
\keyword{Monte Carlo}
\keyword{Lin}
\keyword{Bonferroni}