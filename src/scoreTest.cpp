/*
 Chaitanya Acharya
 
 ---Functions---
 Calculate the combined score test statistic given a list of Y matrices, and a genotype matrix
 1) scoreTest() -> Function to calculate score test statistic given a vector of genotypes
 2) snpOUT()    -> Function to calculate score test statistic given a matrix of genotypes for all samples
 3) GENEout()   -> Function to calculate score test statistic given a matrix of genotypes and a matrix of gene expression data
 
 */

#include <Rcpp.h>
using namespace Rcpp;
#include<ctime>
#include <iostream>
#include<vector>
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
#include <stdint.h>

// [[Rcpp::export]]
double scoreTest(double Eps, double Tau, double k, NumericMatrix Y, NumericVector snp){	
	double nobs = Y.nrow();
	NumericVector Yhat(k);
	NumericMatrix psiMAT(nobs,2);
	
	double mG = mean(snp);
	double mG2 = mean(snp*snp);
	
	for(double l=0.0; l<nobs; l++){
		Yhat = Y(l,_);
		double G = snp[l];
		double sumY = std::accumulate(Yhat.begin(), Yhat.end(), 0.0);
		double sumY2 = std::inner_product(Yhat.begin(), Yhat.end(), Yhat.begin(),0.0);
		double sumProd = (sumY*sumY) - sumY2;
		double sbeta =  ((G-mG) * sumY)/(Eps + k * Tau);
		double dlde = (1/k) * ( (((k-1)*sumY2 - sumProd - (k*(k-1))*Eps) / (Eps*Eps)) + ((sumY*sumY)/((Eps+k*Tau)*(Eps+k*Tau))) - (k/(Eps+k*Tau)) );
		double sgamma = 0.5 * (((G*G)-mG2)*dlde); 		
		
		NumericVector psi(2);
		psi[0] = sbeta;
		psi[1] = sgamma;
		psiMAT(l,_)=psi;
	}
	
	double meanBETA = mean(psiMAT(_,0));
	double meanGAM = mean(psiMAT(_,1));
	NumericVector cv(nobs);
	for(int i=0;i<nobs;i++){
		cv[i] = (psiMAT(i,0)-meanBETA) * (psiMAT(i,1)-meanGAM);
	}
	double cvSUM = 	std::accumulate(cv.begin(), cv.end(), 0.0);
	double cvVAL = cvSUM/(nobs-1);	
 
	NumericVector psi_bta(nobs);
	NumericVector psi_gamma(nobs);
	psi_bta = psiMAT(_,0);
	psi_gamma = psiMAT(_,1);	
	double sum_psi_bta = std::accumulate(psi_bta.begin(), psi_bta.end(),0.0);
	double sum_psi_gamma = std::accumulate(psi_gamma.begin(), psi_gamma.end(),0.0);
	double var_psi_bta = var(psi_bta);
	double var_psi_gamma = var(psi_gamma);
	double U_bta = (1/(nobs*var_psi_bta)) * (sum_psi_bta*sum_psi_bta);
	double U_gamma = (1/(nobs*var_psi_gamma)) * (sum_psi_gamma*sum_psi_gamma);
	double U_psi_ALL = ( -2*cvVAL*sum_psi_bta*sum_psi_gamma + (sum_psi_gamma*sum_psi_gamma*var_psi_bta) + (sum_psi_bta*sum_psi_bta*var_psi_gamma) ) / ( (var_psi_bta*var_psi_gamma)-(cvVAL*cvVAL) );
	double U_psi = U_psi_ALL/nobs;
	
	return U_psi;
}

// [[Rcpp::export]]
NumericVector snpOUT(NumericMatrix geno, double Eps, double Tau, double k, NumericMatrix Y){
	int nsnps = geno.nrow();
	NumericVector U(nsnps);
	for(int i=0; i<nsnps; i++){
		U[i] = scoreTest(Eps,Tau,k,Y,geno(i,_));
	}
	return wrap(U);
}

// [[Rcpp::export]]
NumericMatrix GENEapply(NumericMatrix geno, List Y, NumericVector Eps, NumericVector Tau, double k){
	
	int ngenes = Y.size();
	int nsnps = geno.nrow();
	NumericMatrix GENEout(ngenes,nsnps);
	for(int i =0; i<ngenes; i++){
		//cout<<"Processing gene "<<i+1<<endl;
		Rcpp::NumericMatrix Ymat = Y[i];
		double est_eps = Eps[i];
		double est_tau = Tau[i];
		GENEout(i,_) = snpOUT(geno, est_eps, est_tau, k, Ymat);
	}
	return GENEout;
}
