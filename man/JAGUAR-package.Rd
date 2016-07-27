\name{JAGUAR-package}
\alias{JAGUAR-package}
\docType{package}
\title{
Joint analysis of genotype and group-specific variability using a novel score test to map eQTL
}
\description{
The aim of the package is allow users to apply a novel score test method developed to map eQTL in the presence of multiple correlated groups (for example, tissues) from the same individual. We plan to do this by jointly analyzing all the groups by simultaneously measuring the total shift in the gene expression data due to genotypes and group-specific interaction of the genotypes with the gene expression data. Here is an example of a workflow.

1. We assume that the gene expression data and the genotype data are appropriately preprocessed. Usually, gene expression datasets are long and skinny, i.e. \code{p >> n}. We recommend to partition this gene expression data to run simultaneous analyses on all the partitions to save time. This can be performed using \code{\link{jaguar_slice}}

2. If performing a genome-wide analysis, run \code{\link{jaguar_gwa}} on each gene expression data partition to obtain a matrix of joint score test p-values with genes on rows and SNPs on columns. If performing a cis analysis, run \code{\link{jaguar_cis}} on each gene expression data partition. 

3. Permutation resampling can be performed while running cis analysis and gene-level p-values can be obtained. We do not recommend permutations for genome-wide analysis due to the computational burden.

4. After running a genome-wide analysis, \code{\link{jaguar_process}} function can be used to identify significant gene-SNP pairs based on a predetermined or user-defined threshold value.

5. Power or null simulations can be run using \code{\link{jaguar_sim}} by simulating one gene-SNP pair at a time. 
}
\details{
\tabular{ll}{
Package: \tab JAGUAR\cr
Type: \tab Package\cr
Version: \tab 3.0.1\cr
Date: \tab 2016-07-11\cr
License: \tab GPL-2\cr
}
}
\author{
Chaitanya R. Acharya
Maintainer: Chaitanya Acharya<c.acharya@duke.edu>
}
\references{
Chaitanya R. Acharya, Kouros Owzar, Janice M. McCarthy and Andrew S. Allen; Exploiting expression patterns across multiple tissues to map expression quantitative trait loci. BMC Bioinformatics (2016) 17:257 DOI 10.1186/s12859-016-1123-5
}
