\name{jaguar_process}
\alias{jaguar_process}
\docType{package}
\title{
Obtain significant gene-SNP pairs based on a predetermined threshold value
}
\description{
Function that processes results from running a genome-wide analysis of jaguar and outputs gene-SNP pairs deemed significant by using a predetermined threshold value. It also has an option to print QQ-plot of the p-values from the analysis. 
}
\usage{
jaguar_process(jaguar.out,threshold,plot=FALSE)
}
\arguments{
        \item{jaguar.out}{A Matrix of joint score test p-values with genes on rows and SNPs on columns}
        \item{threshold}{An numeric value representing a threshold value to call for significance}
        \item{plot}{Takes a Boolean value. If 'TRUE', prints a QQ-plot of the p-values from the analysis. In the interests of time and memory management, if there are more than 500,000 gene-SNP pairs in the analysis, only randomly selected 500,000 gene-SNP pairs will be plotted}
}
\value{
        \item{}{A matrix containing three columns -- 1) Genes, 2) SNPs and 3) P-value from the joint score test approach}
}
\author{
Chaitanya R. Acharya
Maintainer: Chaitanya Acharya<c.acharya@duke.edu>
}
\references{
Chaitanya R. Acharya, Kouros Owzar, Janice M. McCarthy and Andrew S. Allen; Exploiting expression patterns across multiple tissues to map expression quantitative trait loci. BMC Bioinformatics (2016) 17:257 DOI 10.1186/s12859-016-1123-5
}
\seealso{
  \code{\link{jaguar_gwa},\link{jaguar_slice},\link{jaguar_sim},\link{jaguar_plotqtl},\link{jaguar_cis}}
}
\examples{
## Example
# 
# Load the example data
data(jaguar_example);

# Genome-wide analysis
Gene_Mat = as.matrix(jaguar_example$GENE_EXP[1:10,])
Geno_Mat = as.matrix(jaguar_example$GENO_MAT)
jag.out = jaguar_gwa(Gene_Mat,Geno_Mat);
dim(jag.out);

# Process results based on a predetermined threshold
result = jaguar_process(jag.out,0.05);
dim(result);
}
\keyword{score test}
\keyword{GWAS}
\keyword{eQTL}
