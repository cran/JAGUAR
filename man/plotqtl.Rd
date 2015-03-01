\name{plotqtl}
\alias{plotqtl}
\title{
Plot of eQTL results
}
\description{
Scatter plot displaying eQTL results with transcript location on the y-axis and SNP location on the x-axis. This plot is an implementation of ePlot function from Wei Sun's eMap R-package. 
}
\usage{
plotqtl(geneID,snpID,gene.chr,gene.pos,snp.chr,snp.pos,scores,chroms)
}
\arguments{
  \item{geneID}{A vector indicating the genes to be mapped}
  \item{snpID}{A vector indicating the SNPs to be mapped}
  \item{gene.chr}{A vector indicating the chromosomal location of the genes to be mapped}
  \item{gene.pos}{A vector indicating the start site of all the genes on the Gene Chip}
  \item{snp.chr}{A vector indicating the chromosomal location of the SNPs to be mapped}
  \item{snp.pos}{A vector indicating the chromosomal location of all the SNPs on the SNP Chip}
  \item{scores}{A vector of p-values of each Gene-SNP pair}
  \item{chroms}{A vector indicating the number of chromosomes to me mapped. Usually, it is 1 to 22 (excluding X and Y chromosomes)}
}
\author{
Chaitanya R. Acharya, Andrew S. Allen
Maintainer: Chaitanya Acharya<c.acharya@duke.edu>
}
\references{
eQTL analysis by Linear Model \url{http://www.bios.unc.edu/~weisun/software/eMap.pdf}

Joint analysis of genotype and tissue-specific variability using a novel score test to map eQTL. Chaitanya R. Acharya, Kouros Owzar and Andrew S. Allen (Manuscript submitted)
}
\seealso{
  \code{\link{jaguar},\link{ProcessJaguarResults},\link{SliceGeneData},\link{jaguar}}
}

\examples{
## NOT RUN
### Read the annotation file of the Gene Chip
#genes = read.table("gene_annotation.txt",header=T,check.names=F)
#eChr = genes$Chromosome
#ePos = genes$StartSite
#
### Read the annotation file of the SNP Chip
#snps = read.table("snp_annotation.txt",header=F)
#mChr = snps$Chr
#mPos = snps$Pos
#
### Read the significant Gene-SNP pairs that are needed to be mapped
#out = ProcessJaguarResults(jaguar.out,threshold=0.05)
#
#geneID = match(out$Probes,genes$Probe_Id)
#markerID = match(out$SNPs,snps$SNP_Id)
#scores = out$P.value
#chroms=1:22
#
#plotqtl(geneID,snpID,gene.chr,gene.pos,snp.chr,snp.pos,scores,chroms)

}
\keyword{Lin}
\keyword{eQTL}
\keyword{plot}
\keyword{Linear Model}
