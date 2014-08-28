/*
 Author: Chaitanya Acharya
 Date: June 17, 2014
 
 Calculates U_psi using var-cov matrix of U_bta and U_gam
 Also calculates U_gam under local null
 
 */

#include <Rcpp.h>
using namespace Rcpp;

/*#include<ctime>
#include<vector>
#include <algorithm>

using namespace std;
*/

// [[Rcpp::export]]
NumericVector sim(double Eps, double Tau, double Eps_g, double Tau_g, double k, NumericMatrix Y, NumericMatrix Y_g, NumericVector snp){

	double nobs = Y.nrow();
	NumericVector Yhat(k);
	NumericVector Yghat(k);
	NumericMatrix psiMAT(nobs,3);
	double mG = mean(snp);
	double mG2 = mean(snp*snp);
	
	for(double l=0.0; l<nobs; l++){
		Yhat = Y(l,_);
		Yghat = Y_g(l,_);
		
		double G = snp[l];
		double sumY = std::accumulate(Yhat.begin(), Yhat.end(), 0.0);
		double sumY2 = std::inner_product(Yhat.begin(), Yhat.end(), Yhat.begin(),0.0);
		
		double sumYg = std::accumulate(Yghat.begin(), Yghat.end(), 0.0);
		double sumY2g = std::inner_product(Yghat.begin(), Yghat.end(), Yghat.begin(),0.0);
		
		double sumProd = (sumY*sumY) - sumY2;
		double sumProd_g = (sumYg*sumYg) - sumY2g;
		
		double sbeta =  ((G-mG) * sumY)/(Eps + k * Tau);
		double dlde = (1/k) * ( (((k-1)*sumY2 - sumProd - (k*(k-1))*Eps) / (Eps*Eps)) + ((sumY*sumY)/((Eps+k*Tau)*(Eps+k*Tau))) - (k/(Eps+k*Tau)) );
		double sgamma = 0.5 * (((G*G)-mG2)*dlde); 
		
		double dlde_g = (1/k) * ( (((k-1)*sumY2g - sumProd_g - (k*(k-1))*Eps_g) / (Eps_g*Eps_g)) + ((sumYg*sumYg)/((Eps_g+k*Tau_g)*(Eps_g+k*Tau_g))) - (k/(Eps_g+k*Tau_g)) );
		double sgamma_g = 0.5 * (((G*G)-mG2)*dlde_g); 
		
		NumericVector psi(3);
		psi[0] = sbeta;
		psi[1] = sgamma;
		psi[2] = sgamma_g;
		
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
	NumericVector psi_g_gamma(nobs);
	
	psi_bta = psiMAT(_,0);
	psi_gamma = psiMAT(_,1);
	psi_g_gamma = psiMAT(_,2);
	
	double sum_psi_bta = std::accumulate(psi_bta.begin(), psi_bta.end(),0.0);
	double sum_psi_gamma = std::accumulate(psi_gamma.begin(), psi_gamma.end(),0.0);
	double var_psi_bta = var(psi_bta);
	double var_psi_gamma = var(psi_gamma);
	double U_bta = (1/ (nobs*var_psi_bta))*(sum_psi_bta*sum_psi_bta);
	double U_gamma = (1/ (nobs*var_psi_gamma)) * (sum_psi_gamma*sum_psi_gamma);
	double U_psi_ALL = ( -2*cvVAL*sum_psi_bta*sum_psi_gamma + (sum_psi_gamma*sum_psi_gamma*var_psi_bta) + (sum_psi_bta*sum_psi_bta*var_psi_gamma) ) / ( (var_psi_bta*var_psi_gamma)-(cvVAL*cvVAL) );
	double U_psi = U_psi_ALL/nobs;
	
	double sum_psi_g_gamma = std::accumulate(psi_g_gamma.begin(), psi_g_gamma.end(),0.0);
	double var_psi_g_gamma = var(psi_g_gamma);
	double U_g_gamma = (1/ (nobs*var_psi_g_gamma)) * (sum_psi_g_gamma*sum_psi_g_gamma);
	
	NumericVector pvals(4);
	pvals[0]=1-R::pchisq(U_bta,1.0,1,0); pvals[1]=1-R::pchisq(U_gamma,1.0,1,0); pvals[2]=1-R::pchisq(U_g_gamma,1.0,1,0);pvals[3]=1-R::pchisq(U_psi,2.0,1,0);
	pvals.names()=CharacterVector::create("U0_bta","U0_gamma","U_gamma","U_psi");
	return wrap(pvals);
}
