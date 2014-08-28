\name{example.data}
\docType{data}
\alias{example.data}
\title{An example gene expression and genotype data}
\description{
This is a list object containing part simulated data set of 15 SNP values in allele dose format and 10 gene expression values for 10 patients in 4 groups.
}
\format{List containing gene expression data as a matrix with genes on rows and samples in columns, genotype data in allele dosage format with SNPs on rows and samples in columns, the number of groups in the study and a vector of sample ids.}
\value{
        \item{Gene}{A matrix of gene expression data with 10 genes and 10 samples in four groups (so a total of 40 samples)}
        \item{Geno}{A matrix of genotype data with 15 SNPs and 10 samples with SNPs in allele dosage format i.e. 0, 1 or 2}
        \item{ngroups}{Number of groups in the study}
        \item{Samples}{A vector of sample ids in the study}
}
\keyword{datasets}
