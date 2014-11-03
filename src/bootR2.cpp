// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-

#include "bootR2.h"


double ave(MatrixXd y) {
    const int n(y.rows());
    const int m(y.cols());
    const double ySum(y.sum());
    const double yAve(ySum/(n*m));
    return(yAve);
}


MatrixXd dev(MatrixXd y) {
    const int n(y.rows());
    const int m(y.cols());
    const double yAve(ave(y));
    MatrixXd devs(y);
    for(int i = 0; i < n; ++i) {
	for(int j = 0; j < m; ++j) {
	    devs(i, j) = y(i, j) - yAve;
	}
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
MatrixXd betaHat(MatrixXd X, MatrixXd y) {
    LLT<MatrixXd> Ch(CholX(X));
    const MatrixXd betaOut = Ch.solve(X.adjoint() * y);
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
double R2(MatrixXd X, MatrixXd y) {
    // const double n(X.rows());  // FIXME:  not needed?
    const MatrixXd coefHat(betaHat(X, y));
    const MatrixXd fitted(X * coefHat);
    const MatrixXd resid(y - fitted);
    const MatrixXd residNull(dev(y));
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
double R2pred(MatrixXd Xt, MatrixXd yt, MatrixXd Xv, MatrixXd yv) {
    const MatrixXd coefHat(betaHat(Xt, yt));
    const MatrixXd fitted(Xv * coefHat);
    const MatrixXd resid(yv - fitted);
    const MatrixXd residNull(dev(yv));
    const double SSerr(pow(resid.norm(), 2));
    const double SStot(pow(residNull.norm(), 2));
    // Rcpp::Rcout << SSerr << "\n\n" << SStot << std::endl;
    const double out(1 - (SSerr/SStot));
    return out;
}


// [[Rcpp::export]]
IntegerVector bootPerm(const int n, const int N) {
    RNGScope scope;
    NumericVector unRound(runif(n, 0, N));
    NumericVector rounded(floor(unRound));
    IntegerVector out = Rcpp::as< IntegerVector >(rounded);
    return out;
}

// [[Rcpp::export]]
MatrixXd shuffleMatrix(const MatrixXd X, const IntegerVector prm) {
    // TODO: test that prm is zero-based
    // const int n(X.rows());
    const int n(prm.size());
    const int m(X.cols());
    MatrixXd Xout(n,m);
    for(int i = 0; i < n; ++i) {
	for(int j = 0; j < m; ++j) {
	    // Rcpp::Rcout << "\n Here! " << m << " " << i << " " << j << std::endl;
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


// // [[Rcpp::export]]
// MatrixXd bootCoef(const MatrixXd X, const MatrixXd y, int nBoot){
//     RNGScope scope;
//     const int n(X.rows());
//     const int p(X.cols());
//     const int m(y.cols());
//     MatrixXd bootBeta(nBoot, p);
//     MatrixXd Xi(X);
//     MatrixXd yi(y);
//     IntegerVector prm(n);
//     MatrixXd betaHati(p, m);
//     for(int i = 0; i < nBoot; ++i) {
// 	prm = bootPerm(n, n);
// 	Xi = shuffleMatrix(X, prm);
// 	yi = shuffleVector(y, prm);
// 	betaHati = betaHat(Xi, yi);
// 	for(int j = 0; j < p; ++j) {
// 	    bootBeta(i, j) = betaHati[j];
// 	}
//     }
//     return bootBeta;
// }


// [[Rcpp::export]]
VectorXd bootR2(const MatrixXd X, const MatrixXd y, int nBoot){
    RNGScope scope;
    const int n(X.rows());
    // const int p(X.cols());
    VectorXd R2s(nBoot);
    MatrixXd Xi(X);
    MatrixXd yi(y);
    IntegerVector prm(n);
    // double R2i(R2(Xi, yi));
    for(int i = 0; i < nBoot; ++i) {
	prm = bootPerm(n, n);
	Xi = shuffleMatrix(X, prm);
	yi = shuffleMatrix(y, prm);
	R2s(i) = R2(Xi, yi);
    }
    return R2s;
}



// [[Rcpp::export]]
VectorXd bootR2pred(const MatrixXd X, const MatrixXd y, int nBoot){
    RNGScope scope;
    const int n(X.rows());
    // const int p(X.cols());
    VectorXd R2s(nBoot);
    MatrixXd Xti(X);
    MatrixXd yti(y);
    MatrixXd Xvi(X);
    MatrixXd yvi(y);
    IntegerVector prmt(n);
    IntegerVector prmv(n);
    // double R2i(R2pred(Xti, yti, Xvi, yvi));
    for(int i = 0; i < nBoot; ++i) {
	prmt = bootPerm(n, n);
	prmv = bootPerm(n, n);
	Xti = shuffleMatrix(X, prmt);
	yti = shuffleMatrix(y, prmt);
	Xvi = shuffleMatrix(X, prmv);
	yvi = shuffleMatrix(y, prmv);
	R2s(i) = R2pred(Xti, yti, Xvi, yvi);
    }
    return R2s;
}


// [[Rcpp::export]]
MatrixXd simExperiment(const MatrixXd Xsamp, const MatrixXd Xpop,
		       const MatrixXd Ysamp, const MatrixXd Ypop,
		       const int pNoi, const int pSig) {
    const int n(Xsamp.rows());
    const int N(Xpop.rows());
    MatrixXd R2s(pNoi+1, 2);
    for(int j = 0; j < (pNoi + 1); ++j) {
	R2s(j, 0) = R2pred(Xsamp.block(0, 0, n, pSig + j - 1),
			   Ysamp,
			   Xpop.block(0, 0, N, pSig + j - 1),
			   Ypop);
	R2s(j, 1) = R2(Xsamp.block(0, 0, n, pSig + j - 1),
		       Ysamp);
    }
    return R2s;
}


// [[Rcpp::export]]
MatrixXd hellinger(MatrixXd X) {
    MatrixXd Xsum = X.rowwise().sum();
    int isPos = Xsum(0, 0) > 0;
    MatrixXd Xhell(X);
    for(int i = 0; i < X.rows(); ++i) {
	isPos = Xsum(i, 0) > 0;
	for(int j = 0; j < X.cols(); ++j) {
	    Xhell(i, j) = isPos ? X(i, j)/Xsum(i, 0) : 0.0;
	}
    }
    return Xhell.cwiseSqrt();
}


// [[Rcpp::export]]
MatrixXd iterSimExperiment(MatrixXd X, MatrixXd Y,
			   const int pNoi, const int pSig,
			   const int nIter, const int n) {
    MatrixXd Yhell = hellinger(Y);
    MatrixXd R2s(nIter, pNoi + 1);
    IntegerVector prm(n);
    const int N(X.rows());
    const int m(Y.cols());
    MatrixXd Xsamp(n, pNoi + pSig);
    MatrixXd Ysamp(n, m);
    for(int i = 0; i < nIter; ++i){
	prm = bootPerm(n, N);
	// Rcpp::Rcout << "\n" <<  << std::endl;
	Xsamp = shuffleMatrix(X, prm);
	Rcpp::Rcout << "\n" << Xsamp << std::endl;
	Ysamp = hellinger(shuffleMatrix(Y, prm));
	// Rcpp::Rcout << "\n" << Ysamp.array() << "\n" << Xsamp.array() << "\n" << std::endl;
	// Rcpp::Rcout << "\n" << Yhell.array() << "\n" << X.array() << "\n" << std::endl;
	for(int j = 0; j < (pNoi + 1); ++j) {
	    // R2s(i, j) = R2pred(Xsamp.block(0, 0, n, pSig + j - 1),
	    // 		       Ysamp,
	    // 		       X.block(0, 0, N, pSig + j - 1),
	    // 		       Yhell);
	    // Rcpp::Rcout << "\n" << Xsamp.block(0, 0, n, pSig + j - 1) << std::endl;
	    R2s(i, j) = R2(Xsamp.block(0, 0, n, pSig + j - 1), Ysamp);
	    // 	R2(Xsamp.block(0, 0, n, pSig + j - 1),
	    // 	   Ysamp);
	}
	//simExpRes.block(i*pNoi, 0, i * (pNoi + 1) + pNoi, 1) =
	//    simExperiment(Xsamp, X, Ysamp, Yhell, pNoi, pSig);
	//Rcpp::Rcout << simExperiment(Xsamp, X, Ysamp, Yhell, pNoi, pSig) << std::endl;
    }
    return Xsamp;
}
