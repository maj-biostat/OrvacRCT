#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <RcppDist.h>
// [[Rcpp::depends(RcppDist)]]

#include <cmath>
#include <algorithm>


// column indices
#define COL_ID            0
#define COL_TRT           1
#define COL_ACCRT         2
#define COL_AGE           3
#define COL_SEROT2        4
#define COL_SEROT3        5
#define COL_PROBT3        6
#define COL_EVTT          7
#define COL_CEN           8
#define COL_OBST          9
#define COL_REASON        10
#define COL_IMPUTE        11
#define COL_REFTIME       12
#define NCOL              13

#define COL_THETA0        0
#define COL_THETA1        1
#define COL_DELTA         2

#define COL_LAMB0         0
#define COL_LAMB1         1
#define COL_RATIO         2

// function prototypes

arma::mat get_trial_dat(const Rcpp::List& cfg);

// Rcpp::List rcpp_clin(arma::mat& d, const Rcpp::List& cfg,
//                      const int look, const int idxsim);
// Rcpp::List rcpp_clin_set_state(arma::mat& d, const int look,
//                                const double fu,
//                                const Rcpp::List& cfg, const int idxsim);
// 
// 
// Rcpp::List rcpp_immu(const arma::mat& d, const Rcpp::List& cfg,
//                      const int look);
// int rcpp_n_obs(const arma::mat& d,
//                const int look,
//                const Rcpp::NumericVector looks,
//                const Rcpp::NumericVector months,
//                const double info_delay);
// Rcpp::List rcpp_lnsero(const arma::mat& d,
//                        const int nobs);
// void rcpp_immu_interim_post(const arma::mat& d,
//                             arma::mat& m,
//                             const int nobs,
//                             const int post_draw,
//                             const Rcpp::List& lnsero);
// Rcpp::List rcpp_immu_interim_ppos(const arma::mat& d,
//                                   const arma::mat& m,
//                                   const int look,
//                                   const int nobs,
//                                   const int nimpute,
//                                   const int post_draw,
//                                   const Rcpp::List& lnsero,
//                                   const Rcpp::List& cfg);
// Rcpp::List rcpp_immu_ppos_test(const arma::mat& d,
//                                const arma::mat& m,
//                                const int look,
//                                const int nobs,
//                                const int nimpute,
//                                const int post_draw,
//                                const Rcpp::List& lnsero,
//                                const Rcpp::List& cfg);
// 
// Rcpp::List rcpp_logrank(const arma::mat& d,
//                         const int look,
//                         const Rcpp::List& cfg);
// void rcpp_outer(const arma::vec& z,
//                 const arma::vec& t,
//                 arma::mat& out);
// arma::vec rcpp_gamma(const int n, const double a, const double b);
// void rcpp_test_1(arma::mat& d);
// void rcpp_test_sub_1(arma::mat& d);
// arma::mat rcpp_test_2(const arma::mat& d) ;
// arma::mat rcpp_test_sub_2(arma::mat& d);
// Rcpp::List rcpp_dotrial(const int idxsim, const Rcpp::List& cfg,
//                         const bool rtn_trial_dat);

// end function prototypes


#endif