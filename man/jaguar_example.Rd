\name{jaguar_example}
\docType{data}
\alias{jaguar_example}
\title{Simulated multi-tissue eQTL data}
\description{
This is a list object containing a simulated eQTL data 
}
\format{List containing gene expression data as a matrix with genes on rows and samples in columns, genotype data in allele dosage format with SNPs on rows and samples in columns, gene and SNP information in BED format.}
\value{
        \item{GENE_EXP}{A matrix of gene expression data with 100 genes and 100 individuals in five groups (so a total of 500 samples)}
        \item{GENO_MAT}{A matrix of genotype data with 1,036 SNPs in 100 individuals with SNPs in allele dosage format i.e. 0, 1 or 2}
        \item{GENE_BED}{Gene location information in BED file format}
        \item{SNP_BED}{SNP location information in SNP file format}
}
\keyword{eQTL}
\keyword{eQTL dataset}
