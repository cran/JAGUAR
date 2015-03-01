## JAGUAR
## Function to compute the score test statistic on gene expression and genotype data with a specified number of groups
## Function looks for matching patient data (same number of patient samples) in all the groups
## Missing values are not allowed

"jaguar" <- function(geneexp,geno,ngroups,write=FALSE){
  
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
    fit = lmer(x ~ 0+tissue+(1|ind),REML=F)
    x = matrix(x,nobs,ngroups)
    est.eps = attr(VarCorr(fit),"sc")^2;
    est.tau = VarCorr(fit)[[1]][1];
    Ynew = apply(x,1,function(xx) xx-colMeans(x));
    return(list("Eps"=est.eps,"Tau"=est.tau,"Y"=Ynew));
  })
  Eps = laply(GeneOBJ,function(x){x$Eps})
  Tau = laply(GeneOBJ,function(x){x$Tau})
  Y = llply(GeneOBJ,function(x){x$Y})
  
  U = GENEapply(as.matrix(geno),Y,Eps,Tau,ngroups)
  rownames(U)=geneID
  colnames(U)=snpID
  if(write) write.table(U,"jaguar_out.txt",sep="\t",col.names=NA,quote=F)
  return(U)
}

### SliceGeneData
### Function to slice the gene expression data
### Creates subdirectories to store partitioned gene expression data

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
    dir.create(paste(path,"dir",i,sep=""));
    setwd(paste(path,"dir",i,sep=""));
    print(getwd())
    write.table(geneexp[dataSplit[[i]],],"GeneExp_matrix.txt",sep="\t",col.names=NA,quote=F);
    setwd("../");
  }
  cat( "Please check for the subdirectories under the listed path: ",path, "\n" );
}

## LIN
## Function to identify a new p-value threshold to adjust for multiplicity using Danyu Lin's efficient Monte Carlo approach
## More specifically, it calculates the minimum p-value for all the genes over many realizations of the joint score test statistic
## and writes the final result into a file "mcMIN.txt" at the base directory

"lin" <- function(geneexp,geno,ngroups,mc.real=5000,parallel=FALSE,ncores=2,snp.slice=10000,write=TRUE){
  if(mean(unique(colnames(geneexp)) %in% colnames(geno)) != 1) stop("Sample names in gene expression and genotype data do not match")	
  nobs = ncol(geno);
  ind = as.factor(rep(1:nobs,ngroups));
  tissue = as.factor(rep(1:ngroups,each=nobs));
  set.seed(12345)
  G	= matrix(rnorm(mc.real*nobs),nrow=mc.real,ncol=nobs)
  geneexp = suppressWarnings(data.frame(t(geneexp),check.names=F));
  temp = llply(geneexp,as.numeric);
  GeneOBJ = llply(temp,function(x){
    fit = lmer(x ~ 0+tissue+(1|ind),REML=F)
    x = matrix(x,nobs,ngroups)
    est.eps = attr(VarCorr(fit),"sc")^2;
    est.tau = VarCorr(fit)[[1]][1];
    Ynew = apply(x,1,function(xx) xx-colMeans(x));
    return(list("Eps"=est.eps,"Tau"=est.tau,"Y"=Ynew));
  })
  Eps = laply(GeneOBJ,function(x){x$Eps})
  Tau = laply(GeneOBJ,function(x){x$Tau})
  Y = llply(GeneOBJ,function(x){x$Y})
  if(parallel){
    cat("Parallelizing the analysis over all the SNPs \n")
    if(ncores<2) stop("When parallelizing the analysis, at least 2 cores are necessary!")
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
    minP.gene = apply(ldply(1:length(GenoMat.list),.parallel=TRUE,function(i) getMinP(G,Y,as.matrix(GenoMat.list[[i]]),Eps,Tau,ngroups,TRUE)),2,max)
  }else{
    minP.gene = getMinP(G,Y,geno,Eps,Tau,ngroups,TRUE);
  }
  if(write) write.table(minP.gene,"mcMIN.txt",sep="\t",col.names=F,row.names=F,quote=F)
  return(minP.gene)
}

## calcThreshold
## Function to calculate the Lin Monte Carlo adjusted and Bonferroni thresholds
## If subdirectories were created on partitioned gene expression data to run the 'lin' function on each partition, 
##    this function will grab the "mcMIN.txt" from such analyses from each subdirectory

"calcThreshold" <- function(nsnps,ngenes,method=c("bonferroni","lin"),path,alpha=0.05){
  
  if(method=="bonferroni"){
    print("Calculating threshold by Bonferroni method")
    if(is.null(nsnps) || is.null(ngenes)) stop("Please provide the total number of snps and genes in the analysis")
    threshold = alpha/(nsnps*ngenes);
  }else{
    print("Calculating threshold by Lin's efficient MC method")
    dirs = dir(path=path)[file.info(dir(path=path,full.names=T))$isdir]
    dirs = dirs[grep("dir",dirs)]
    minMAT = vector(mode="list",length=length(dirs))
    for(i in 1:length(minMAT)){
      setwd(paste(path,"/dir",i,sep=""));
    	minMAT[[i]] = read.delim("mcMIN.txt",header=F)
    	setwd("../")
    }
    lin_MAT = do.call("cbind",minMAT);
    minPval = apply(lin_MAT,1,min);
    threshold = as.numeric(quantile(as.numeric(minPval),probs=alpha));
  }
  return(threshold);
}

## ProcessJaguarResults
## Function to process results from running JAGUAR based on a predetermined threshold value
## Also plots a QQ-plot of the score test pvalues by randomly sampling 

"ProcessJaguarResults" <- function(jaguar.out,threshold,plot=FALSE){

  if(class(jaguar.out)!="matrix") stop("Wrong argument type passed!");
  if(plot){
    U.vec = as.vector(jaguar.out);
    if(length(U.vec)>50000){
      U.vec.temp = U.vec[sample(1:length(U.vec),50000,replace=F)];
    }else{
      U.vec.temp = U.vec[sample(1:length(U.vec),length(U.vec),replace=F)];
    }
    theoretical.quantiles = (1:length(U.vec.temp))/(1+length(U.vec.temp)); 
    data.quantiles = sort(U.vec.temp);
    plot(-log10(theoretical.quantiles),-log10(data.quantiles),main="QQ-plot of the score test pvals",xlab="Theoritical Quantiles",ylab="Data Quantiles");
    abline(0,1,col="red")
  }
  w = which(jaguar.out<threshold,arr.ind=T);
  val = jaguar.out[jaguar.out<threshold];
  genes = rownames(jaguar.out)[w[,1]];
  snps = colnames(jaguar.out)[w[,2]];
  eQTL = cbind(genes,snps,val);
  colnames(eQTL)[1:3] = c("Genes","SNPs","P-value");
  return(eQTL);
}

## jaguarSIM
## Function to run some simulations
## Tests "Global" H_0: \beta=0; \gamma=0 and the "Local" H_0: \gamma=0

"jaguarSIM" = function(nobs = 500, k = 5, tau = 1, eps = 1,PVEg = 0,bta = 0,maf = 0.10){
  gamma = ((eps+tau)*PVEg) / (100-PVEg);
  mu = rep(0.5,k);
  repeat{
    snp = rbinom(nobs,2,maf)
    if(sum(snp)>1){ break }
  }
  bkg=rnorm(k,0,gamma)
  Y = do.call("rbind",lapply(1:nobs,function(i) mu + bta*snp[i] + rnorm(1,0,tau) + bkg*snp[i] + rnorm(k,0,eps)))
  data = data.frame("IND"=as.factor(rep(1:nobs,each=k)),"Gene"=as.vector(t(Y)),"Geno" = rep(snp,each=k),"Tissue"=as.factor(rep(1:k,nobs)))
  
  # Fit the model at global null
  fit0 = lmer(Gene~0+Tissue+(1|IND),data,REML=F)
  est0.eps = attr(VarCorr(fit0),"sc")^2; est0.tau = VarCorr(fit0)[[1]][1];
  Y0new = apply(Y,1,function(x) x-colMeans(Y))
  joint_ST = scoreTest(est0.eps,est0.tau,k,Y0new,snp)
  
  # Fit the model at local null for U_gamma
  fit = lmer(Gene~0+Geno+Tissue+(1|IND),data,REML=F)
  est.mu = as.numeric(fixef(fit)); est.eps = attr(VarCorr(fit),"sc")^2;
  est.tau = VarCorr(fit)[[1]][1];	J = model.matrix(fit);
  Ynew = matrix(data$Gene-(J%*%est.mu),k,nobs)
  gamma_ST = gamma_test(est.eps,est.tau,k,Ynew,snp)
  
  return(c("GammaScoreTest"=gamma_ST,"JointScoreTest"=joint_ST))
}

## plotqtl
## Plotting function to visualze eQTL results
## This is a slightly modified plotting function originally written by Wei Sun as a part of eMap R-Package

plotqtl = function(geneID,snpID,gene.chr,gene.pos,snp.chr,snp.pos,scores,chroms){
  
  if(length(geneID) == 0){
    stop("length(geneID)=0\n")
  }
  if(length(geneID) != length(snpID)){
    stop("length(geneID) != length(snpID)\n")
  }
  if(length(geneID) != length(scores)){
    stop("length(geneID) != length(scores)\n")
  }
  valid.chroms = c(1:90, "X", "Y")
  wrongChr = chroms[!(chroms %in% valid.chroms)]
  if(length(wrongChr)>0){
    stop(wrongChr, " are not valid chromosome labels\n")
  }
  
  todrop = which(!(gene.chr %in% chroms))
  gene.chr[todrop] = NA
  gene.pos[todrop] = NA
  gene.chr[gene.chr=="X"] = 99
  gene.chr[gene.chr=="Y"] = 100
  gene.chr = as.integer(gene.chr)
  gene.pos = as.numeric(gene.pos)
  
  todrop = which(!(snp.chr %in% chroms))
  snp.chr[todrop] = NA
  snp.pos[todrop] = NA
  snp.chr[snp.chr=="X"] = 99
  snp.chr[snp.chr=="Y"] = 100
  snp.chr = as.integer(snp.chr)
  snp.pos = as.numeric(snp.pos)
  
  chrs = union(unique(gene.chr), unique(snp.chr))
  chrs = sort(as.numeric(chrs))
  num.chrs = chrs[chrs <= 90]
  chr.chrs = c(num.chrs, "X", "Y")
  max.num = max(num.chrs)
  if(max.num > length(num.chrs)){
    str = "there is no SNP/gene location information for some chromosomes\n"
    stop(str)
  }
  gene.chr[gene.chr==99]  = max.num+1
  gene.chr[gene.chr==100] = max.num+2
  snp.chr[snp.chr==99]  = max.num+1
  snp.chr[snp.chr==100] = max.num+2
  chreMax = tapply(gene.pos, gene.chr, max, na.rm=TRUE, simplify = FALSE)
  chrmMax = tapply(gene.pos, gene.chr, max, na.rm=TRUE, simplify = FALSE)
  chrMax = numeric(length(chrs))
  for(i in 1:length(chrs)){
    ch = as.character(i)
    chrMax[i] = max(chreMax[[ch]], chrmMax[[ch]])
  }
  temp = data.frame(Gene_ID=geneID, Marker_ID=snpID)
  nChr = length(chrMax)
  chrLen = c(0, cumsum(chrMax))
  ep = gene.pos + chrLen[gene.chr]
  mp = snp.pos + chrLen[snp.chr]
  ymax = chrLen[nChr+1]
  bdr1 = -0.016*ymax
  bdr2 = -0.006*ymax
  
  par(mar=c(3,4,0,0))
  plot(c(bdr1,ymax*1.05), c(bdr1,ymax*1.05), type="n", xlab="",ylab="", main="", xaxt="n", yaxt="n", bty="n")
  mtext("SNP Location", side=1, line=1)
  mtext("Transcript Location", side=2, line=1)
  gpos = ep[temp$Gene_ID]
  mpos = mp[temp$Marker_ID]
  points(mpos, gpos, col="black",pch=20, cex=1)
  
  if(nChr>1){
    nchr.plot = floor(nChr/2)
    rect(rep(bdr1,nchr.plot), chrLen[seq(1,nChr,by=2)],rep(bdr2,nchr.plot), chrLen[seq(2,nChr+1,by=2)],border=NA, col="orange")
    rect(chrLen[seq(1,nChr,by=2)], rep(bdr1,nchr.plot),chrLen[seq(2,nChr+1,by=2)], rep(bdr2,nchr.plot),border=NA, col="orange")
  }
  kp = seq(1,nChr,by=2)
  ats = 0.5*(chrLen[-1] + chrLen[-length(chrLen)])
  mtext(chr.chrs[kp], at=ats[kp], side=1, line=-0.5, cex=0.8)
  mtext(chr.chrs[kp], at=ats[kp], side=2, line=-0.5, cex=0.8)
  lines(c(0,ymax*1.01), rep(ymax*1.01,2), lty=2)
  lines(rep(ymax*1.01,2), c(0,ymax*1.01), lty=2)  
}
