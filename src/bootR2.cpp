// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-

#include "bootR2.h"


double ave(VectorXd y) {
    const double n(y.size());    
    const double ySum(y.sum());
    const double yAve(ySum/n);
    return(yAve);
}


VectorXd dev(VectorXd y) {
    const int n(y.size());    
    const double yAve(ave(y));
    VectorXd devs(y);
    for(int i = 0; i < n; ++i) {
	devs[i] = y[i] - yAve;
    }
    return devs;
}

LLT<MatrixXd> CholX(MatrixXd X) {
    const int p(X.cols());
    const MatrixXd XtX(MatrixXd(p, p).setZero().selfadjointView<Lower>().rankUpdate(X.adjoint()));
    LLT<MatrixXd> out(XtX.selfadjointView<Lower>());
    return out;
}

// [[Rcpp::export]]
VectorXd betaHat(MatrixXd X, VectorXd y) {
    LLT<MatrixXd> Ch(CholX(X));
    const VectorXd betaOut = Ch.solve(X.adjoint() * y);
    return betaOut;
}

// [[Rcpp::export]]
// MatrixXd betaHatMat(MatrixXd X, MatrixXd Y) {
//     Rcpp::Rcout << "here we are" << std::endl;
//     LLT<MatrixXd> Ch(CholX(X));
//     const MatrixXd betaOut = Ch.solve(X.adjoint() * Y);
//     return betaOut;
// }


// [[Rcpp::export]]
double R2(MatrixXd X, VectorXd y) {
    // const double n(X.rows());  // FIXME:  not needed?
    const VectorXd coefHat(betaHat(X, y));
    const VectorXd fitted(X * coefHat);
    const VectorXd resid(y - fitted);
    const VectorXd residNull(dev(y));
    const double SSerr(pow(resid.norm(), 2));
    const double SStot(pow(residNull.norm(), 2));
    // Rcpp::Rcout << SSerr << "\n\n" << SStot << std::endl;
    const double out(1 - (SSerr/SStot));
    return out;
}

// // [[Rcpp::export]]
// double R2Mat(MatrixXd X, MatrixXd Y) {
//     // const double n(X.rows());  // FIXME:  not needed?
//     const MatrixXd coefHat(betaHatMat(X, Y));
//     const MatrixXd fitted(X * coefHat);
//     const MatrixXd resid(Y - fitted);
//     const VectorXd residNull(dev(y));
//     const double SSerr(pow(resid.norm(), 2));
//     const double SStot(pow(residNull.norm(), 2));
//     // Rcpp::Rcout << SSerr << "\n\n" << SStot << std::endl;
//     const double out(1 - (SSerr/SStot));
//     return out;
// }


// [[Rcpp::export]]
double R2pred(MatrixXd Xt, VectorXd yt, MatrixXd Xv, VectorXd yv) {
    const VectorXd coefHat(betaHat(Xt, yt));
    const VectorXd fitted(Xv * coefHat);
    const VectorXd resid(yv - fitted);
    const VectorXd residNull(dev(yv));
    const double SSerr(pow(resid.norm(), 2));
    const double SStot(pow(residNull.norm(), 2));
    // Rcpp::Rcout << SSerr << "\n\n" << SStot << std::endl;
    const double out(1 - (SSerr/SStot));
    return out;
}


// [[Rcpp::export]]
IntegerVector bootPerm(const int n) {
    RNGScope scope;
    NumericVector unRound(runif(n, 0, n));
    NumericVector rounded(floor(unRound));
    IntegerVector out = Rcpp::as< IntegerVector >(rounded);
    return out;
}

// [[Rcpp::export]]
MatrixXd shuffleMatrix(const MatrixXd X, const IntegerVector prm) {
    // TODO: test that prm is zero-based
    const int n(X.rows());
    const int m(X.cols());
    MatrixXd Xout(n,m);
    for(int i = 0; i < n; ++i) {
	for(int j = 0; j < m; ++j) {
	    Xout(i,j) = X(prm[i],j);
	}
    }
    return Xout;
}

// [[Rcpp::export]]
VectorXd shuffleVector(const VectorXd y, const IntegerVector prm) {
    // TODO: test that prm is zero-based
    const int n = y.size();
    VectorXd yOut(n);
    for(int i = 0; i < n; ++i) {
	yOut[i] = y[prm[i]];
    }
    return yOut;
}


// [[Rcpp::export]]
NumericVector rUnif(const int n) {
    RNGScope scope;
    NumericVector out(runif(1, 0, n));
    return out;
}


// [[Rcpp::export]]
MatrixXd bootCoef(const MatrixXd X, const VectorXd y, int nBoot){
    RNGScope scope;
    const int n(X.rows());
    const int p(X.cols());
    MatrixXd bootBeta(nBoot, p);
    MatrixXd Xi(X);
    VectorXd yi(y);
    IntegerVector prm(n);
    VectorXd betaHati(p);
    for(int i = 0; i < nBoot; ++i) {
	prm = bootPerm(n);
	Xi = shuffleMatrix(X, prm);
	yi = shuffleVector(y, prm);
	betaHati = betaHat(Xi, yi);
	for(int j = 0; j < p; ++j) {
	    bootBeta(i, j) = betaHati[j];
	}
    }
    return bootBeta;
}


// [[Rcpp::export]]
VectorXd bootR2(const MatrixXd X, const VectorXd y, int nBoot){
    RNGScope scope;
    const int n(X.rows());
    // const int p(X.cols());
    VectorXd R2s(nBoot);
    MatrixXd Xi(X);
    VectorXd yi(y);
    IntegerVector prm(n);
    // double R2i(R2(Xi, yi));
    for(int i = 0; i < nBoot; ++i) {
	prm = bootPerm(n);
	Xi = shuffleMatrix(X, prm);
	yi = shuffleVector(y, prm);
	R2s(i) = R2(Xi, yi);
    }
    return R2s;
}



// [[Rcpp::export]]
VectorXd bootR2pred(const MatrixXd X, const VectorXd y, int nBoot){
    RNGScope scope;
    const int n(X.rows());
    // const int p(X.cols());
    VectorXd R2s(nBoot);
    MatrixXd Xti(X);
    VectorXd yti(y);
    MatrixXd Xvi(X);
    VectorXd yvi(y);
    IntegerVector prmt(n);
    IntegerVector prmv(n);
    // double R2i(R2pred(Xti, yti, Xvi, yvi));
    for(int i = 0; i < nBoot; ++i) {
	prmt = bootPerm(n);
	prmv = bootPerm(n);
	Xti = shuffleMatrix(X, prmt);
	yti = shuffleVector(y, prmt);
	Xvi = shuffleMatrix(X, prmv);
	yvi = shuffleVector(y, prmv);
	R2s(i) = R2pred(Xti, yti, Xvi, yvi);
    }
    return R2s;
}

