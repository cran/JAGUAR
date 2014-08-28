\name{calcThreshold}
\alias{calcThreshold}
\docType{package}
\title{
Calculates various threshold values to control for family-wise error rate (FWER).
}
\description{
Function to calculate Bonferroni threshold value and Lin threshold value from the results of Lin's Monte Carlo method. These threshold values control for FWER at the 0.05 level. 
}
\usage{
calcThreshold(file.list,nsnps,ngenes,alpha=0.05)
}
\arguments{
        \item{file.list}{A list of vectors obtained by running ComputeLinMax() on several partitions of the gene expression data set.}
        \item{nsnps}{An integer indicating the total number of SNPs in the study}
        \item{ngenes}{An integer indicating the total number of genes in the study}
        \item{alpha}{The nominal p-value indicating the false positive rate}
}
\value{
        \item{}{A list object containing the Bonferroni and Lin threshold values}
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
  \code{\link{ComputeLinMax},\link{jaguar},\link{ProcessJaguarResults},\link{SliceGeneData},\link{jaguarSIM}}
}
\examples{
## Example
# 
# Load the example data (not an ideal dataset)
data(example.data);

# Set the parameter values
GeneExp = as.matrix(example.data$Gene); 
Geno = as.matrix(example.data$Geno);
ngroups=example.data$ngroups;

# Run the function
out.vec = ComputeLinMax(GeneExp,Geno,ngroups);

# Calculate the threshold values
thresholds = calcThreshold(list(out.vec),15,10);

}
\keyword{Monte Carlo}
\keyword{Lin}
\keyword{Bonferroni}