\name{jaguar_cis}
\alias{jaguar_cis}
\docType{package}
\title{
Perform cis-eQTL analysis
}
\description{
Computes p-value from our joint score test in a cis framework to map group-specific expression quantitative trait loci (eQTL) that tests for the shifts in gene expression patterns due to genotype and variability among tissues in a mixed effects model framework. A gene-level p-value is computed using a permutation-resampling scheme in order to investigate if a gene has at least one eQTL across all the groups. 
}
\usage{
jaguar_cis(geneexp,genomat,snp.bed,gene.bed,cisDist=100000,nperm=10000,seed=100)
}
\arguments{
        \item{geneexp}{ A matrix of gene expression data with gene on rows and patient samples on columns. There has to be equal number of samples in each group. Samples (columns in the geneexp object) with missing gene expression values for any group/tissue MUST be included in the data}
        \item{genomat}{A matrix of genotype data recoded as single allele dosage number (i.e. 0, 1 or 2) with rows representing SNPs and columns representing samples}
        \item{snp.bed}{BED file format of SNP description. For more information, see \url{http://genome.ucsc.edu/FAQ/FAQformat.html#format1} or the example data}
        \item{gene.bed}{BED file format of Gene description. For more information, see \url{http://genome.ucsc.edu/FAQ/FAQformat.html#format1} or the example data}
        \item{cisDist}{cis distance is defined as the maximum absolute distance between the gene and a SNP. Default value is 100Kb}
        \item{nperm}{Number of permutations. Default value is 10,000. Note that if it is 0, no permutations are performed}
        \item{seed}{Seed value for permutations}
}
\value{
	If permutations are performed, a list object with the following components is returned -

        \item{results}{A matrix containing gene-level permutation p-values} 
	\item{cisSNPs}{A list object containing the cis-SNPs corresponding to the genes in results}

	If permutations are not performed -
	\item{}{A matrix containing raw unadjusted p-values for all gene-cisSNP pairs}
}


\author{
Chaitanya R. Acharya
Maintainer: Chaitanya Acharya<c.acharya@duke.edu>
}
\references{
Chaitanya R. Acharya, Kouros Owzar, Janice M. McCarthy and Andrew S. Allen; Exploiting expression patterns across multiple tissues to map expression quantitative trait loci. BMC Bioinformatics (2016) 17:257 DOI 10.1186/s12859-016-1123-5
}
\seealso{
  \code{\link{jaguar_process},\link{jaguar_slice},\link{jaguar_sim},\link{jaguar_gwa},\link{jaguar_plotqtl}}
}
\keyword{score test statistic }
\examples{
# Load the example data
data(jaguar_example);
Gene = jaguar_example$GENE_EXP
SNP = jaguar_example$GENO_MAT
gene_loc = jaguar_example$GENE_BED
snp_loc = jaguar_example$SNP_BED
# Run a cis analysis with no permutations
out = jaguar_cis(Gene,SNP,snp_loc,gene_loc,nperm=0);
length(out)
}
\keyword{joint score test}
\keyword{genotype}
\keyword{eQTL}
\keyword{interaction}
