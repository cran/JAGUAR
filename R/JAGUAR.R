## JAGUAR
## Function to compute the score test statistic on gene expression and genotype data with a specified number of groups
## Function looks for matching patient data (same number of patient samples) in all the groups
## There should not be any missing gene expression or genotype data

"jaguar" <- function(geneexp,geno,ngroups){
  
  if(mean(unique(colnames(geneexp)) %in% colnames(geno)) != 1) stop("Sample names in gene expression and genotype data do not match")
  cat("Running JAGUAR... \n")
  geneID = rownames(geneexp);
  snpID = rownames(geno);
  nobs = ncol(geno);
  ind = as.factor(rep(1:nobs,ngroups));
  tissue = as.factor(rep(1:ngroups,each=nobs));
  
  geneexp = suppressWarnings(data.frame(t(geneexp),check.names=F));
  temp = llply(geneexp,as.numeric);
  GeneOBJ = llply(temp,function(x){
    fit = lmer(as.numeric(x)~tissue+(1|ind),REML=F)
    est.eps = attr(VarCorr(fit),"sc")^2;
    est.tau = VarCorr(fit)[[1]][1];
    Y.hat = as.matrix(resid(fit))
    Y.new = matrix(Y.hat,nrow=nobs,ncol=ngroups)
    return(list("Eps"=est.eps,"Tau"=est.tau,"Y"=Y.new));
  })
  Eps = laply(GeneOBJ,function(x){x$Eps})
  Tau = laply(GeneOBJ,function(x){x$Tau})
  Y = llply(GeneOBJ,function(x){x$Y})
  
  U = GENEapply(as.matrix(geno),Y,Eps,Tau,ngroups)
  rownames(U)=geneID
  colnames(U)=snpID
  return(U)
}

### SliceGeneData
### Function to slice the gene expression data
### Creates subdirectoriesto store partitioned gene expression data

"SliceGeneData" <- function(geneexp,size,path=getwd()){
  
  ngenes = size;
  samples = colnames(geneexp);
  tot.size = dim(geneexp)[1];
  part = round(tot.size/ngenes);
  cat( "Slicing the gene expression data in to ", part ," partitions", "\n" );
  output = matrix(1:((part-1)*ngenes),nrow=part-1,byrow=TRUE);
  dataSplit = split(output, 1:nrow(output));
  dataSplit[[length(dataSplit)+1]]= (length(dataSplit)*ngenes+1):tot.size;
  for(i in 1:length(dataSplit)){
    dir.create(paste(path,"/dir",i,sep=""));
    setwd(paste(path,"/dir",i,sep=""));
    print(getwd())
    write.table(geneexp[dataSplit[[i]],],"GeneExp_matrix.txt",sep="\t",col.names=NA,quote=F);
    setwd("../");
  }
  cat( "Please check for the subdirectories under the listed path: ",path, "\n" );
}

## ComputeLinMax
## Function to calculate maximum of the score test statistic over all the genes for all the samples in the given gene expression matrix
## Writes the final result into a file "mcMAX.txt"

"ComputeLinMax" <- function(geneexp,geno,ngroups,mc.real=5000,parallel=FALSE,ncores=2,snp.slice=10000){

  if(mean(unique(colnames(geneexp)) %in% colnames(geno)) != 1) stop("Sample names in gene expression and genotype data do not match")
  
  nobs = ncol(geno);
  ind = as.factor(rep(1:nobs,ngroups));
  tissue = as.factor(rep(1:ngroups,each=nobs));
  set.seed(12345)
  G = matrix(rnorm(mc.real*nobs),nrow=mc.real,ncol=nobs)
  geneexp = suppressWarnings(data.frame(t(geneexp),check.names=F));
  temp = llply(geneexp,as.numeric);
  GeneOBJ = llply(temp,function(x){
    fit = lmer(as.numeric(x)~0+tissue+(1|ind),REML=F)
    est.mu = as.numeric(fixef(fit));
    est.eps = attr(VarCorr(fit),"sc")^2;
    est.tau = VarCorr(fit)[[1]][1];
    A = model.matrix(fit);
    Y.hat = as.matrix(resid(fit))
    Y.new = matrix(Y.hat,nrow=nobs,ncol=ngroups)
    return(list("Eps"=est.eps,"Tau"=est.tau,"Y"=Y.new));
  })
  Eps = laply(GeneOBJ,function(x){x$Eps})
  Tau = laply(GeneOBJ,function(x){x$Tau})
  Y = llply(GeneOBJ,function(x){x$Y})
  if(parallel){
    cat("Parallelizing the analysis over all the SNPs \n")
    if(ncores<2) stop("When parallelizing the analysis, at least 2 cores are necessary")
    registerDoParallel(ncores)
    size.genoData = dim(geno)[1];
    part = round(size.genoData/snp.slice);
    output = matrix(1:((part-1)*snp.slice),nrow=part-1,byrow=TRUE);
    dataSplit = split(output, 1:nrow(output));
    dataSplit[[length(dataSplit)+1]]= (length(dataSplit)*snp.slice+1):size.genoData;
    GenoMat.list = vector(mode="list",length=length(dataSplit));
    for(i in 1:length(dataSplit)){
      GenoMat.list[[i]]=geno[dataSplit[[i]],]
    }
    max.gene = apply(ldply(1:length(GenoMat.list),.parallel=TRUE,function(i)getMAX(G,Y,as.matrix(GenoMat.list[[i]]),Eps,Tau,ngroups)),2,max)
  }else{
    cat("Computing Lin threshold with no parallel execution \n")
    max.gene = getMAX(G,Y,geno,Eps,Tau,ngroups);
  }
  return(max.gene);

}

## calcThreshold
## Function to calculate the so called "Lin MC" and Bonferroni thresholds
## If subdirectories were created on partitioned gene expression data to run the LinMax function on each partition, 
##    this function will grab the results from such analyses from each subdirectory

"calcThreshold" <- function(file.list,nsnps,ngenes,alpha=0.05){
  
  if(class(file.list)!="list") stop("Wrong argument type passed!")
  if(is.null(nsnps) || is.null(ngenes)) stop("Please provide the total number of snps and genes in the analysis")
  mc_MAT = do.call("cbind",file.list);
  maxT = apply(mc_MAT,1,max);
  threshold_lin = as.numeric(quantile(as.numeric(maxT),probs=1-alpha));
  threshold_bon = qchisq(alpha/(nsnps*ngenes),df=2,lower.tail=F);
  return(list("linstat"=threshold_lin,"bonferroni"=threshold_bon));

}

## ProcessJaguarResults
## Function to process results from running JAGUAR based on a predetermined threshold value
## Also plots a QQ-plot of the score test pvalues by randomly sampling 

"ProcessJaguarResults" <- function(jaguar.out,threshold,plot=F){

  if(class(jaguar.out)!="matrix") stop("Wrong argument type passed!");
  if(plot){
    U.vec = as.vector(jaguar.out);
    if(length(U.vec)>50000){
      U.vec.temp = U.vec[sample(1:length(U.vec),50000,replace=F)];
    }else{
      U.vec.temp = U.vec[sample(1:length(U.vec),length(U.vec),replace=F)];
    }
    pval = rep(0,length(U.vec.temp));
    for(i in 1:length(pval)){
      pval[i] = 1-pchisq(U.vec.temp[i],df=2);
    }
    theoretical.quantiles = (1:length(pval))/(1+length(pval)); data.quantiles = sort(pval);
    plot(-log10(theoretical.quantiles),-log10(data.quantiles),main="QQ-plot of the score test pvals",xlab="Theoritical Quantiles",ylab="Data Quantiles");
    abline(0,1,col="red")
  }
  w = which(jaguar.out>threshold,arr.ind=T);
  val = jaguar.out[jaguar.out>threshold];
  genes = rownames(jaguar.out)[w[,1]];
  snps = colnames(jaguar.out)[w[,2]];
  eQTL = cbind(genes,snps,val);
  colnames(eQTL)[1:3] = c("Genes","SNPs","ScoreTestStat");
  return(eQTL);
  
}

## jaguarSIM
## Function to run some simulations

jaguarSIM = function(nobs = 500, k = 4, tau = 1, eps = 1,gamma = 0,bta = 0,maf = 0.10){
  mu = rep(0.5,k);
  repeat{
    snp = rbinom(nobs,2,maf)
    if(sum(snp)>1){ break }
  }
  Y = do.call("rbind",lapply(1:nobs,function(i) mu + bta*snp[i] + rnorm(1,0,tau) + rnorm(k,0,gamma)*snp[i] + rnorm(k,0,eps)))
  data = data.frame("IND"=as.factor(rep(1:nobs,each=k)),"Gene"=as.vector(t(Y)),"Geno" = rep(snp,each=k),"Tissue"=as.factor(rep(1:k,nobs)))
  
  # Fit the model at global null (H0: gamma=0; bta=0)
  fit0 = lmer(Gene~0+Tissue+(1|IND),data,REML=F)
  est0.eps = attr(VarCorr(fit0),"sc")^2;
  est0.tau = VarCorr(fit0)[[1]][1];
  Y0.new = matrix(resid(fit0),nrow=nobs,ncol=k,byrow=T)
  
  # Fit the model at local null (H0: gamma=0)
  # Implements variance component score test
  fit = lmer(Gene~0+Geno+Tissue+(1|IND),data,REML=F)
  est.eps = attr(VarCorr(fit),"sc")^2;
  est.tau = VarCorr(fit)[[1]][1];
  Y.new = matrix(resid(fit),nrow=nobs,ncol=k,byrow=T)
  
  # Run the score test
  st = sim(Eps=est0.eps,Tau=est0.tau,k=k, Eps_g=est.eps, Tau_g=est.tau, Y=Y0.new, Y_g=Y.new, snp=snp)
  return(st)
}