/*
 Author: Chaitanya Acharya
 Date: June 17, 2014
 
C++ Function to calculate the Lin threshold
NOTE: cout is disabled for compilation
 */

#include <Rcpp.h>
using namespace Rcpp;
/*#include<ctime>
#include <iostream>
#include<vector>
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
#include <stdint.h>
*/

// [[Rcpp::export]]
NumericVector mc_scoreTest(double Eps, double Tau, double k, NumericMatrix Y, NumericVector snp, NumericMatrix mcMAT){
  
	double nobs = Y.nrow();
	NumericVector Yhat(k);
	NumericMatrix psiMAT(nobs,2);
	int rsim = mcMAT.nrow();
	
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
	NumericVector psi_bta(nobs), psi_bta_new(nobs);
	NumericVector psi_gamma(nobs), psi_gamma_new(nobs);
	psi_bta = psiMAT(_,0);
	psi_gamma = psiMAT(_,1);
	NumericVector U_psi_new(rsim);	
	for(int r=0; r<rsim; r++){
		psi_bta_new = psi_bta * mcMAT(r,_);
		psi_gamma_new = psi_gamma * mcMAT(r,_);
		double sum_psi_bta_new = std::accumulate(psi_bta_new.begin(), psi_bta_new.end(),0.0);
		double sum_psi_gamma_new = std::accumulate(psi_gamma_new.begin(), psi_gamma_new.end(),0.0);
		double var_psi_bta_new = var(psi_bta_new);
		double var_psi_gamma_new = var(psi_gamma_new);
		double U_bta_new = (1/(nobs*var_psi_bta_new)) * (sum_psi_bta_new*sum_psi_bta_new);
		double U_gamma_new = (1/(nobs*var_psi_gamma_new)) * (sum_psi_gamma_new*sum_psi_gamma_new);
		U_psi_new[r] = U_bta_new + U_gamma_new;
	}
	return wrap(U_psi_new);
}

// [[Rcpp::export]]
NumericVector getMAX(NumericMatrix mcMAT, List Y, NumericMatrix geno, NumericVector Eps, NumericVector Tau, double k){
	
	//bool verbose = v;
	int ngenes = Y.size();
	int nsnps = geno.nrow();
	int rsim = mcMAT.nrow();
	int nobs = geno.ncol();
	//cout<<"Num of genes: "<<ngenes<<" Num of snps: "<<nsnps<<endl;
	NumericMatrix maxV(ngenes,rsim);
	
	for(int j=0; j<ngenes; j++){
	/*	if(verbose){
			cout<<"Processing gene "<<j+1<<endl;
		} */
		NumericMatrix Ymat = Y[j];
		double est_eps = Eps[j];
		double est_tau = Tau[j];
		NumericMatrix GenoMAT(nsnps,rsim);

		for(int i=0; i<nsnps; i++){
			NumericVector snp(nobs);
			snp = geno(i,_);
			GenoMAT(i,_) = mc_scoreTest(est_eps,est_tau,k,Ymat,snp,mcMAT);
		}
		int ncol = GenoMAT.ncol();
		NumericVector out(ncol);
		for (int col = 0; col < ncol; col++){
			out[col]=Rcpp::max(GenoMAT(_, col)); 
		}
		maxV(j,_) = out;
	}
	int maxV_cols = maxV.ncol();
	NumericVector geneMax(maxV_cols);
	for(int i=0; i<maxV_cols; i++){
		geneMax[i] = Rcpp::max(maxV(_,i));
	}
	return wrap(geneMax);
}
