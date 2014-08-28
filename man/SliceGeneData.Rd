\name{SliceGeneData}
\alias{SliceGeneData}
\docType{package}
\title{
Slice gene expression data into multiple partitions
}
\description{
Function to 1) create sub-directories, 2) slice gene expression data into partitions of predetermined size, and 3) sliced gene expression partitions are deposited into each sub-directory
}
\usage{
SliceGeneData(geneexp,size,path=getwd())
}
\seealso{
  \code{\link{calcThreshold},\link{ComputeLinMax},\link{jaguar},\link{ProcessJaguarResults},\link{jaguarSIM}}
}
\arguments{
        \item{geneexp}{ A matrix of gene expression data with gene on rows and patient samples on columns. Missing values not allowed. For now, there has to be an equal number of samples in each group.}
        \item{size}{Integer indicating the size of each slice of gene expression data.}
        \item{path}{Location for the sub-directories. Please give the full path. Default is set to the current directory.}
}
\author{
Chaitanya R. Acharya, Andrew S. Allen
Maintainer: Chaitanya Acharya<c.acharya@duke.edu>
}
\references{
Chaitanya R. Acharya, Kouros Owzar and Andrew S. Allen; Joint analysis of genotype and tissue-specific variability using a novel score test to map eQTL (Manuscript submitted)
}
\keyword{ package }
\examples{
# Set the size of the partition
# size = 100; ## Indicates the number of genes in each partitioned gene exp data
#
# Assuming that the path is the default getwd()
# SliceGeneData(geneexp,100)
}