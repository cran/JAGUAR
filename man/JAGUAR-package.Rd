\name{JAGUAR-package}
\alias{JAGUAR-package}
\docType{package}
\title{
Joint analysis of genotype and group-specific variability using a novel score test to map eQTL
}
\description{
The aim of the package is allow users to apply a novel score test method developed to map eQTL in the presence of multiple correlated groups (for example, tissues) from the same individual. We plan to do this by jointly analyzing all the groups by simultaneously measuring the total shift in the gene expression data due to genotypes and group-specific interaction of the genotypes with the gene expression data. Here is an example of a workflow.

1. We assume that the gene expression data and the genotype data are appropriately preprocessed. Usually, gene expression datasets are long and skinny, i.e. \code{p >> n}. We recommend to partition this gene expression data to run simultaneous analyses on all the partitions to save time. This can be performed using \code{\link{SliceGeneData}}

2. Run \code{\link{jaguar}} on each gene expression data partition to obtain a matrix of joint score test p-values with genes on rows and SNPs on columns.

3. A threshold p-value needs to be established to correct for multiplicity.

4. A threshold value can be computed from the above results using \code{\link{calcThreshold}} by combining results from all gene expression data analyses. We provide two methods to adjust for multiplicity. One is the standard Bonferroni correction and other's Danyu Lin's efficient Monte Carlo method. It is important to note that the results from all the partitions of the gene expression data are required to compute the threshold value by Lin's method.

5. This threshold value can be used to call the significant gene-SNP pairs in the analysis by running \code{\link{ProcessJaguarResults}} on the output from running {\code{\link{jaguar}}}.

6. Power or null simulations can be run using \code{\link{jaguarSIM}} by simulating one gene and one SNP at a time. 
}
\details{
\tabular{ll}{
Package: \tab JAGUAR\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2014-03-01\cr
License: \tab GPL-2\cr
}
}
\author{
Chaitanya R. Acharya, Andrew S. Allen
Maintainer: Chaitanya Acharya<c.acharya@duke.edu>
}
\references{
An efficient Monte Carlo approach to assessing the statistical significance in genomic studies. Lin, D.Y. Bioinformatics. 21(6) 2005. 

Joint analysis of genotype and tissue-specific variability using a novel score test to map eQTL. Chaitanya R. Acharya, Kouros Owzar and Andrew S. Allen (Manuscript submitted)
}