// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "bootR2.h"

// mult
NumericVector mult(NumericVector a, NumericVector b);
RcppExport SEXP bootR2_mult(SEXP aSEXP, SEXP bSEXP) {
    BEGIN_RCPP
	SEXP __sexp_result;
    {
	Rcpp::RNGScope __rngScope;
	Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
	Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
	NumericVector __result = mult(a, b);
	PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
    END_RCPP
}


// trans
MatrixXd trans(MatrixXd X);
RcppExport SEXP bootR2_trans(SEXP XX) {
    BEGIN_RCPP
	SEXP __sexp_result;
    { 
	const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
	MatrixXd __result = trans(X);
	PROTECT(__sexp_result = wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
    END_RCPP
}


// crossProd
MatrixXd crossProd(MatrixXd X);
RcppExport SEXP bootR2_crossProd(SEXP XX) {
    BEGIN_RCPP
	SEXP __sexp_result;
    {
	const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
	MatrixXd __result = crossProd(X);
	PROTECT(__sexp_result = wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
    END_RCPP
}

// betahat
VectorXd betaHat(MatrixXd X, VectorXd y);
RcppExport SEXP bootR2_betaHat(SEXP XX, SEXP yy) {
    BEGIN_RCPP
	SEXP __sexp_result;
    {
	const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
	const Map<VectorXd> y(as<Map<VectorXd> >(yy));
	VectorXd __result = betaHat(X, y);
	PROTECT(__sexp_result = wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
    END_RCPP
}


// R2
double R2(MatrixXd X, VectorXd y);
RcppExport SEXP bootR2_R2(SEXP XX, SEXP yy) {
    BEGIN_RCPP
	SEXP __sexp_result;
    {
	const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
	const Map<VectorXd> y(as<Map<VectorXd> >(yy));
	double __result = R2(X, y);
	PROTECT(__sexp_result = wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
    END_RCPP
}


// R2pred
double R2pred(MatrixXd Xt, VectorXd yt, MatrixXd Xv, VectorXd yv);
RcppExport SEXP bootR2_R2pred(SEXP XXt, SEXP yyt, SEXP XXv, SEXP yyv) {
    BEGIN_RCPP
	SEXP __sexp_result;
    {
	const Map<MatrixXd> Xt(as<Map<MatrixXd> >(XXt));
	const Map<VectorXd> yt(as<Map<VectorXd> >(yyt));
	const Map<MatrixXd> Xv(as<Map<MatrixXd> >(XXv));
	const Map<VectorXd> yv(as<Map<VectorXd> >(yyv));
	double __result = R2pred(Xt, yt, Xv, yv);
	PROTECT(__sexp_result = wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
    END_RCPP
}



// rUnif
NumericVector rUnif(int n);
RcppExport SEXP bootR2_rUnif(SEXP nn) {
    BEGIN_RCPP
	SEXP __sexp_result;
    {
	Rcpp::traits::input_parameter< int >::type n(nn);
	NumericVector __result = rUnif(n);
	PROTECT(__sexp_result = wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
    END_RCPP
}


// bootPerm
IntegerVector bootPerm(int n);
RcppExport SEXP bootR2_bootPerm(SEXP nn) {
    BEGIN_RCPP
	SEXP __sexp_result;
    {
	Rcpp::traits::input_parameter< int >::type n(nn);
	IntegerVector __result = bootPerm(n);
	PROTECT(__sexp_result = wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
    END_RCPP
}

IntegerVector order_(NumericVector x);
RcppExport SEXP bootR2_order_(SEXP xx) {
    BEGIN_RCPP
	SEXP __sexp_result;
    {
	Rcpp::traits::input_parameter< NumericVector >::type x(xx);
	IntegerVector __result = order_(x);
	PROTECT(__sexp_result = wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
    END_RCPP
}



// shuffleMatrix
MatrixXd shuffleMatrix(const MatrixXd X, const IntegerVector prm);
RcppExport SEXP bootR2_shuffleMatrix(SEXP XX, SEXP prm_) {
    BEGIN_RCPP
	SEXP __sexp_result;
    {
	const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
	Rcpp::traits::input_parameter< IntegerVector >::type prm(prm_);
	MatrixXd __result = shuffleMatrix(X, prm);
	PROTECT(__sexp_result = wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
    END_RCPP
}

// shuffleVector
VectorXd shuffleVector(const VectorXd y, const IntegerVector prm);
RcppExport SEXP bootR2_shuffleVector(SEXP yy, SEXP prm_) {
    BEGIN_RCPP
	SEXP __sexp_result;
    {
	// Rcpp::traits::input_parameter< NumericVector >::type y(yy);
	const Map<VectorXd> y(as<Map<VectorXd> >(yy));
	Rcpp::traits::input_parameter< IntegerVector >::type prm(prm_);
	VectorXd __result = shuffleVector(y, prm);
	PROTECT(__sexp_result = wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
    END_RCPP
}


// bootCoef
MatrixXd bootCoef(const MatrixXd X, const VectorXd y, const int nBoot);
RcppExport SEXP bootR2_bootCoef(SEXP XX, SEXP yy, SEXP nnBoot) {
    BEGIN_RCPP
	SEXP __sexp_result;
    {
	const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
	const Map<VectorXd> y(as<Map<VectorXd> >(yy));
	Rcpp::traits::input_parameter< int >::type nBoot(nnBoot);
	MatrixXd __result = bootCoef(X, y, nBoot);
	PROTECT(__sexp_result = wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
    END_RCPP
}



// bootR2
VectorXd bootR2(const MatrixXd X, const VectorXd y, const int nBoot);
RcppExport SEXP bootR2_bootR2(SEXP XX, SEXP yy, SEXP nnBoot) {
    BEGIN_RCPP
	SEXP __sexp_result;
    {
	const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
	const Map<VectorXd> y(as<Map<VectorXd> >(yy));
	Rcpp::traits::input_parameter< int >::type nBoot(nnBoot);
	VectorXd __result = bootR2(X, y, nBoot);
	PROTECT(__sexp_result = wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
    END_RCPP
}




// bootR2pred
VectorXd bootR2pred(const MatrixXd X, const VectorXd y, const int nBoot);
RcppExport SEXP bootR2_bootR2pred(SEXP XX, SEXP yy, SEXP nnBoot) {
    BEGIN_RCPP
	SEXP __sexp_result;
    {
	const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
	const Map<VectorXd> y(as<Map<VectorXd> >(yy));
	Rcpp::traits::input_parameter< int >::type nBoot(nnBoot);
	VectorXd __result = bootR2pred(X, y, nBoot);
	PROTECT(__sexp_result = wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
    END_RCPP
}



