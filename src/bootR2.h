// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace RcppEigen;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::LLT;
using Eigen::Lower;

typedef Eigen::Map<MatrixXd> MMap;
