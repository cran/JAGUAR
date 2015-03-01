/*
Author: Chaitanya R. Acharya
Updated: March 1, 2015

This file has all the Rcpp functions that are part of the JAGUAR package
*/
#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include<vector>
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
#include <stdint.h>

// [[Rcpp::export]]

NumericVector RowSums(NumericMatrix x) {
  
	int nrow = x.nrow(), ncol = x.ncol();
	NumericVector out(nrow);
	
	for (int i = 0; i < nrow; i++) {
		double total = 0;
		for (int j = 0; j < ncol; j++) {
			total += x(i, j);
		}
		out[i] = total;
	}
	return out;
}

//[[Rcpp::export]]
NumericVector RowMin(NumericMatrix x){
  NumericVector out(x.nrow());
	for(int i=0; i<x.nrow(); i++){
		out[i] = min(x(i,_));
	}
	return wrap(out);
}

// [[Rcpp::export]]

double scoreTest(double Eps, double Tau, double k, NumericMatrix Y, NumericVector snp){	
	
	double nobs = Y.ncol();
	NumericVector Yhat(k);
	double V1 = (Eps+(k-1)*Tau) / ( (Eps*Eps)+k*Eps*Tau);
	double V2 = - (Tau)/(Eps*Eps+k*Eps*Tau);
	NumericVector Ug(k); 
	NumericVector Ubta(nobs);
	NumericMatrix Ugam(k,nobs);
	double mG = mean(snp);

	for(int i=0; i<nobs; i++){
		Yhat = Y(_,i);
		double G = snp[i]-mG;
		double Ysum = sum(Yhat);
		double sbeta =  (G * Ysum)/(Eps + k * Tau);
		for(int t=0; t<k; t++){
			Ug[t] = ( (V1*Yhat[t])+(V2*(Ysum-Yhat[t])) ) * G;
		}
		Ugam(_,i)=Ug;
		Ubta[i]=sbeta;
	}

	double Ugamma = 0.5 * sum(RowSums(Ugam)*RowSums(Ugam));
	double U2beta = sum(Ubta)*sum(Ubta);
	double snp2 = sum((snp-mG)*(snp-mG));
	double Db = snp2 * (V1+(k-1)*V2) * k;
	double Dg = 0.5 * sum((snp-mG)*(snp-mG)*V1*k);
	double Dg2 = 0.25 * snp2 * snp2 * (V1*V1+(k-1)*V2*V2) * k;
	double varB = 2 * (Db*Db);
	double varG = 2 * Dg2;
	double cov = (snp2 * (V1+(k-1)*V2))*(snp2 * (V1+(k-1)*V2))*k;	
	double abeta = (varB - cov)/(varB+varG-2*cov); 
	double agam = (varG - cov)/(varB+varG-2*cov);	
	double Upsi = abeta*U2beta + agam*Ugamma;
	double pval;
	Db = abeta * Db;
	Dg = agam * Dg;
	double c1 = Db+Dg; double c2 = Db*Db + (agam*agam * 0.5 * varG);
	if(c1 == 0 || c2==0){
		pval = 1;
	}else{
		double b = (c1*c1)/c2;
		double a = c2/c1;
		pval = 1 - R::pchisq(Upsi/a,b,1,0);
	}
	return pval;
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
		Rcpp::Rcout<<"Processing gene "<<i+1<<std::endl;
		Rcpp::NumericMatrix Ymat = Y[i];
		double est_eps = Eps[i];
		double est_tau = Tau[i];
		GENEout(i,_) = snpOUT(geno, est_eps, est_tau, k, Ymat);
	}
	return GENEout;
}

// [[Rcpp::export]]

double gamma_test(double Eps, double Tau, double k, NumericMatrix Y, NumericVector snp){  
	
  double nobs = Y.ncol();
	NumericVector Yhat(k);
	double V1 = (Eps+(k-1)*Tau) / ( (Eps*Eps)+k*Eps*Tau);
	double V2 = - (Tau)/(Eps*Eps+k*Eps*Tau);
	NumericVector Ug(k); NumericMatrix Ugam(k,nobs);
	double mG = mean(snp);

	for(int i=0; i<nobs; i++){
		Yhat = Y(_,i);
		double Ysum = sum(Yhat);
		double G = snp[i]-mG;
		for(int t=0; t<k; t++){
			Ug[t] = ( (V1*Yhat[t])+(V2*(Ysum-Yhat[t])) ) * G;
		}
		Ugam(_,i)=Ug;
	}
	
	double Ugamma = 0.5 * sum(RowSums(Ugam)*RowSums(Ugam));
	double snp2 = sum((snp-mG)*(snp-mG));
	double c1 = 0.5 * sum((snp-mG)*(snp-mG)*V1*k);
	double c2 = 0.25 * snp2 * snp2 * (V1*V1+(k-1)*V2*V2) * k;
	double pval;
	if(c1 == 0 || c2==0){
		pval = 1;
	}else{
		double b = (c1*c1)/c2;
		double a = c2/c1;
		pval = 1 - R::pchisq(Ugamma/a,b,1,0);
	}
	return pval;

}

// [[Rcpp::export]]

List jag_param(double Eps, double Tau, double k, NumericMatrix Y, NumericVector snp){
  double nobs = Y.ncol();
	NumericVector Yhat(k);
	double V1 = (Eps+(k-1)*Tau) / ( (Eps*Eps)+k*Eps*Tau);
	double V2 = - (Tau)/(Eps*Eps+k*Eps*Tau);
	NumericVector Ug(k); 
	NumericVector Ubta(nobs);
	NumericMatrix Ugam(k,nobs);
	double mG = mean(snp);
	
	for(int i=0; i<nobs; i++){
		Yhat = Y(_,i);
		double G = snp[i]-mG;
		double Ysum = sum(Yhat);
		double sbeta =  (G * Ysum)/(Eps + k * Tau);
		for(int t=0; t<k; t++){
			Ug[t] = ( (V1*Yhat[t])+(V2*(Ysum-Yhat[t])) ) * G;
		}
		Ugam(_,i)=Ug;
		Ubta[i]=sbeta;
	}
	double snp2 = sum((snp-mG)*(snp-mG));
	double Db = snp2 * (V1+(k-1)*V2) * k;
	double Dg = 0.5 * sum((snp-mG)*(snp-mG)*V1*k);
	double Dg2 = 0.25 * snp2 * snp2 * (V1*V1+(k-1)*V2*V2) * k;
	double varB = 2 * (Db*Db);
	double varG = 2 * Dg2;
	double cov = (snp2 * (V1+(k-1)*V2))*(snp2 * (V1+(k-1)*V2))*k;	
	double abeta = (varB - cov)/(varB+varG-2*cov); 
	double agam = (varG - cov)/(varB+varG-2*cov);	
	Db = abeta * Db;
	Dg = agam * Dg;
	double c1 = Db+Dg; double c2 = Db*Db + (agam*agam * 0.5 * varG);
	double b = (c1*c1)/c2;
	double a = c2/c1;
	return List::create(Ubta,Ugam,abeta,agam,a,b);
}

// [[Rcpp::export]]
NumericVector getMinP(NumericMatrix mcMAT, List Y, NumericMatrix geno, NumericVector Eps, NumericVector Tau, double k,bool v){
  bool verbose = v;
	int ngenes = Y.size();
	int nsnps = geno.nrow();
	int B = mcMAT.nrow();
	int nobs = geno.ncol();
	Rcpp::Rcout<<"Num of genes: "<<ngenes<<" Num of snps: "<<nsnps<<std::endl;
	NumericMatrix snpTEMP(B,nsnps);
	NumericMatrix geneTEMP(B,ngenes);
	NumericVector Ubta_star(nobs);
	NumericVector Ugam_star(k);
	NumericVector Gvec(nobs);
	NumericVector Upval_star(B);
	
	for(int j=0; j<ngenes; j++){
		if(verbose){
      Rcpp::Rcout<<"Processing Gene "<<j+1<<std::endl;
    }
		NumericMatrix Ymat = Y[j];
		double est_eps = Eps[j];
		double est_tau = Tau[j];

		for(int i=0; i<nsnps; i++){
			NumericVector snp=geno(i,_);
			List p = jag_param(est_eps,est_tau,k,Ymat,snp);
			NumericVector Ubta = p[0]; NumericMatrix UgamNEW = p[1];
			double abeta = as<double>(p[2]); double agam = as<double>(p[3]); double a1 = as<double>(p[4]); double a2 = as<double>(p[5]);
			for(int b=0; b<B; b++){
				Gvec=mcMAT(b,_);
				Ubta_star = Gvec*Ubta;
				for(int t=0; t<k; t++){
				 Ugam_star[t] = sum(Gvec * UgamNEW(t,_));
				}
				double Ugamma_star = 0.5 * sum(Ugam_star*Ugam_star);
				double U2beta_star = sum(Ubta_star)*sum(Ubta_star);
				double Upsi_star = abeta*U2beta_star + agam*Ugamma_star;
				Upval_star[b] = 1 - R::pchisq(Upsi_star/a1,a2,1,0);
			}
			snpTEMP(_,i) = Upval_star;
		}
		geneTEMP(_,j)=RowMin(snpTEMP);
	}
	return wrap(RowMin(geneTEMP));
 }
