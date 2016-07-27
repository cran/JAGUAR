\name{jaguar_slice}
\alias{jaguar_slice}
\docType{package}
\title{
Slice gene expression data into multiple partitions
}
\description{
Function to 1) create sub-directories, 2) slice gene expression data into partitions of predetermined size, and 3) sliced gene expression partitions are deposited into each sub-directory
}
\usage{
jaguar_slice(geneexp,size,path=getwd())
}
\arguments{
        \item{geneexp}{ A matrix of gene expression data with gene on rows and patient samples on columns. There has to be equal number of samples in each group. Samples (columns in the geneexp object) with missing gene expression values for any group/tissue MUST be included in the data}
        \item{size}{Integer indicating the size of each slice of gene expression data.}
        \item{path}{Location for the sub-directories. Please give the full path. Default is set to the current directory.}
}
\author{
Chaitanya R. Acharya
Maintainer: Chaitanya Acharya<c.acharya@duke.edu>
}
\references{

Chaitanya R. Acharya, Kouros Owzar, Janice M. McCarthy and Andrew S. Allen; Exploiting expression patterns across multiple tissues to map expression quantitative trait loci. BMC Bioinformatics (2016) 17:257 DOI 10.1186/s12859-016-1123-5

}
\keyword{ partition }
\keyword{ slice }
\keyword{ gene expression }
\keyword{ sub-directory }
\seealso{
  \code{\link{jaguar_gwa},\link{jaguar_process},\link{jaguar_sim},\link{jaguar_plotqtl},\link{jaguar_cis}}
}
\examples{
# Set the size of the partition
# size = 100; ## Indicates the number of genes in each partitioned gene exp data
#
# Assuming that the path is the default getwd()
# jaguar_slice(geneexp,size)
}
