\name{ProcessJaguarResults}
\alias{ProcessJaguarResults}
\docType{package}
\title{
Obtain significant gene-SNP pairs based on a predetermined threshold value
}
\description{
Function that processes results from running jaguar and outputs gene-SNP pairs deemed significant by using a predetermined threshold value. It also has an option to print QQ-plot of the p-values from the analysis. 
}
\usage{
ProcessJaguarResults(jaguar.out,threshold,plot=FALSE)
}
\arguments{
        \item{jaguar.out}{A Matrix of joint score test statistic values with genes on rows and SNPs on columns}
        \item{threshold}{An integer representing a threshold value to call for significance}
        \item{plot}{Takes a Boolean value. If 'TRUE', prints a QQ-plot of the p-values from the analysis. In the interests of time and memory management, if there are more than 500,000 gene-SNP pairs in the analysis, only randomly selected 500,000 gene-SNP pairs will be plotted}
}
\value{
        \item{}{A matrix containing three columns -- 1) Genes, 2) SNPs and 3) P-value from the joint score test approach}
}
\author{
Chaitanya R. Acharya, Andrew S. Allen
Maintainer: Chaitanya Acharya<c.acharya@duke.edu>
}
\references{
Joint analysis of genotype and tissue-specific variability using a novel score test to map eQTL. Chaitanya R. Acharya, Kouros Owzar and Andrew S. Allen (Manuscript submitted)
}
\seealso{
  \code{\link{calcThreshold},\link{lin},\link{jaguar},\link{SliceGeneData},\link{jaguarSIM}}
}
\examples{
## Example
# 
# Load the example data
data(example.data);
#
# Set the parameters
GeneExp = as.matrix(example.data$Gene); 
Geno = as.matrix(example.data$Geno);
ngroups=example.data$ngroups;
#
# Run JAGUAR
jag.out = jaguar(GeneExp,Geno,ngroups,write=FALSE);
#
# Process JAGUAR output with a predetermined threshold value
# computed from ComputeLinThreshold()
#
result = ProcessJaguarResults(jag.out,0.50);
dim(result);
head(result);
}
\keyword{score test}
\keyword{GWAS}
\keyword{eQTL}