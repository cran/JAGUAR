# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

jag_fun <- function(Eps, Tau, k, Y, snp, R) {
    .Call('JAGUAR_jag_fun', PACKAGE = 'JAGUAR', Eps, Tau, k, Y, snp, R)
}

GENEapply <- function(geno, Y, Eps, Tau, k, R, display_progress = TRUE) {
    .Call('JAGUAR_GENEapply', PACKAGE = 'JAGUAR', geno, Y, Eps, Tau, k, R, display_progress)
}

rowsumscpp <- function(x) {
    .Call('JAGUAR_rowsumscpp', PACKAGE = 'JAGUAR', x)
}

jagSIM <- function(Eps, Tau, k, Y, snp) {
    .Call('JAGUAR_jagSIM', PACKAGE = 'JAGUAR', Eps, Tau, k, Y, snp)
}

vcSIM <- function(Eps, Tau, k, Y, snp) {
    .Call('JAGUAR_vcSIM', PACKAGE = 'JAGUAR', Eps, Tau, k, Y, snp)
}

cis_eqtl <- function(SNP, Y, perm_mat, Eps, Tau, k, R, v, display_progress = TRUE) {
    .Call('JAGUAR_cis_eqtl', PACKAGE = 'JAGUAR', SNP, Y, perm_mat, Eps, Tau, k, R, v, display_progress)
}

