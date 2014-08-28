\name{ComputeLinMax}
\alias{ComputeLinMax}
\docType{package}
\title{
Computes maximum of the score test statistic using an efficient Monte Carlo approach.
}
\description{
Computes the weighted sum of independent standard normal random variables by multiplying the score test statistic with multiple realizations of normal random samples and taking the max over all the genes. Using an efficient Monte Carlo approach, we approximate the joint distribution of the novel joint score test statistic and then use the Monte Carlo distribution to evaluate the error rates of the statistic. This approach of multiple hypotheses correction in the case of genomic studies has been shown to be computationally inexpensive when compared with standard permutation based approaches.
}
\usage{
ComputeLinMax(geneexp,geno,ngroups,mc.real=5000,parallel=FALSE,ncores=2,snp.slice=10000)
}
\arguments{
        \item{geneexp}{A matrix of gene expression data with gene on rows and patient samples on columns. Missing values not allowed. There has to be equal number of samples in each group.}
        \item{geno}{A matrix of genotype data recoded as single allele dosage number (i.e. 0, 1 or 2) with rows representing SNPs and columns representing samples }
        \item{ngroups}{An integer representing the number of groups in the data}
        \item{mc.real}{An integer representing the number of Monte Carlo realizations of the score test statistic. Higher the number, higher the accuracy and longer the computation time}
        \item{parallel}{A boolean value indicating whether the analysis should be parallelized}
        \item{ncores}{An integer representing the number of cores used for parallel execution. The number must be at least 2. This option is required only when the analysis is parallelized}
        \item{snp.slice}{An integer representing the partition size of the genotype/SNP data. This number basically represents the number of SNPs analyzed at a time. This option is used only when the parallel execution of the analysis is activated.}
}
\value{
        \item{}{A vector of length equal to the number of Monte Carlo realizations containing maximum of the weighted sum of independent standard normal random variables over all the genes}
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
  \code{\link{calcThreshold},\link{jaguar},\link{ProcessJaguarResults},\link{SliceGeneData},\link{jaguarSIM}}
}
\examples{
# Example
#
# Load the example data
data(example.data)

# Set the parameters
GeneExp = as.matrix(example.data$Gene); 
Geno = as.matrix(example.data$Geno); 
k = example.data$ngroups;

#Compute the maximum of the score test statistic over all the genes
#   and SNPs with no parallel execution. We recommend at least 5000 mc realizations
lin.stat.out = ComputeLinMax(GeneExp,Geno,k,mc.real=100)
length(lin.stat.out)
}
\keyword{Lin}
\keyword{Monte Carlo}
\keyword{Joint distribution}