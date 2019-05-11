#include <RcppDist.h>
// [[Rcpp::depends(RcppDist)]]

#include <cmath>
#include <algorithm>

#include "trial.h"
#include "simulator.h"


#define _DEBUG 0

#if _DEBUG
#define DBG( os, msg )                                \
(os) << "DBG: " << __FILE__ << "(" << __LINE__ << ") "\
     << msg << std::endl
#else
#define DBG( os, msg )
#endif
   
#define _INFO  1
   
#if _INFO
#define INFO( os, i, msg )                                \
   (os) << "INFO: " << __FILE__ << "(" << __LINE__ << ") "\
        << " sim = " << i << " " << msg << std::endl
#else
#define INFO( os, i, msg )
#endif



// [[Rcpp::export]]
arma::mat get_trial_dat(const Rcpp::List& cfg) {

  int n = (int)cfg["n_stop"];
  arma::mat d = arma::zeros(n, NCOL);

  for(int i = 0; i < n; i++){
    
    d(i, COL_ID) = i;
    d(i, COL_TRT) = ((i-1)%2 == 0) ? 0 : 1;
    // poisson proc
    if(i == 0) {
      d(i, COL_ACCRT) = R::rexp(1/(double)cfg["accrual"])  ;  
    } else {
      d(i, COL_ACCRT) = d(i-1, COL_ACCRT) + R::rexp(1/(double)cfg["accrual"])  ;  
    }
    
    d(i, COL_AGE) = R::runif((double)cfg["age_months_lwr"], (double)cfg["age_months_upr"]);
    
    d(i, COL_SEROT2) = R::rbinom(1, cfg["baseline_prob_sero"]);
    d(i, COL_SEROT3) = d(i, COL_SEROT2);
    d(i, COL_PROBT3) = d(i, COL_TRT) * (double)cfg["delta_sero_t3"];
    
    if(d(i, COL_SEROT2) == 0 && d(i, COL_TRT) == 1){
      d(i, COL_SEROT3) = R::rbinom(1, d(i, COL_PROBT3));
    }
    
    // tte - the paramaterisation of rexp uses SCALE NOTE RATE!!!!!!!!!!!
    // event time is the time from randomisation (not birth) at which first
    // medical presentation occurs
    if(d(i, COL_TRT) == 0){
      d(i, COL_EVTT) = R::rexp(1/(double)cfg["b0_tte"])  ;
    } else {
      double beta = (double)cfg["b0_tte"] + (double)cfg["b1_tte"];
      d(i, COL_EVTT) = R::rexp(1/beta)  ;
    }
    
  }
  
  d.col(COL_CEN) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_OBST) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_REASON) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_IMPUTE) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_REFTIME) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  
  return d;
}


// [[Rcpp::export]]
arma::mat get_interims(const arma::mat& d, const Rcpp::List& cfg){
  
  // interims occur at index INT_I_END at time INT_T_END
  // an interim uses kids from 0 to INT_I_END with 
  // kids at index INT_I_START to INT_I_END being introduced
  // for this interim.
  
  int n = (int)cfg["n_stop"];
  int n_start = (int)cfg["n_start"] - 1;
  
  // I resize later
  arma::mat intrms = arma::zeros(300, 4);
  
  intrms(0, INT_I_START) = 0;
  intrms(0, INT_I_END) = n_start;                // 69
  intrms(0, INT_T_START) = 0;
  intrms(0, INT_T_END) = d(n_start, COL_ACCRT);  // accrual for 69

  intrms(1, INT_I_START) = intrms(0, INT_I_END) + 1;    // 70
  intrms(1, INT_T_START) = intrms(0, INT_T_END);        // accrual for 69
  
  int k = 2;
  int i = 0;
  
  for(i = intrms(1, INT_I_START) + 1; i < n; i++){
    
    // every three months or next 50, whichever sooner
    // start the clock from the accrual of the previous kid
    if(d(i, COL_ACCRT) - intrms(k-1, INT_T_START) >= (double)cfg["int_mnth"]  |  
      d(i, COL_ID) - intrms(k-1, INT_I_START)  >= (int)cfg["int_n"]-1) {

      intrms(k-1, INT_I_END) = i;
      // start the clock from the accrual of the previous kid
      intrms(k-1, INT_T_END) = d(i, COL_ACCRT);
      
      intrms(k, INT_I_START) = intrms(k-1, INT_I_END) + 1;
      // start the clock from the accrual of the previous kid
      intrms(k, INT_T_START) = intrms(k-1, INT_T_END);

      k++;
    }
  }
  
  // Rcpp::Rcout  << " k-2 " << k - 2 << " intrms row k-2 " << intrms.row(k-2) << std::endl;
  // Rcpp::Rcout  << " k-1 " << k - 1 << " intrms row k-1 " << intrms.row(k-1) << std::endl;
  // Rcpp::Rcout  << " k   " << k  << " intrms row k   " << intrms.row(k) << std::endl;
  
  int z = 0;
  int ilast = (int)intrms(k, INT_I_END);

  while(ilast == 0 | ilast == n-1){
    
    // Rcpp::Rcout << "entry at intrms(" << k << ", INT_I_END) = " 
    //   << " is " << intrms(k, INT_I_END) << std::endl;
    
    k = k - 1;
    ilast = (int)intrms(k, INT_I_END);
    
    z++;
    if(z > 10) break;
  }
  
  intrms.resize(k+1, 4);
  
  return intrms;
  
}

// [[Rcpp::export]]
void nothing_much() {
  
  // Trial t;
  
  // int stop_ven_samp = 0;
  // int stop_immu_fut = 0;
  // int stop_clin_fut = 0;
  // int stop_clin_es = 0;
  // int inconclu = 0;
  // int nmaxsero = 250;
  // int nstartclin = 200;
  // int immu_ss = 0;
  // int clin_ss = 0;
  // 
  // bool i_final_win = 0;
  // bool c_final_win = 0;
  
  // Rcpp::Rcout << "test 1 do_immu(10) " << do_immu(10) == true << std::endl;
  
  // bool do_immu(10);
  // bool do_clin(int n_current);
  // 
  // int get_nmaxsero();
  // int get_nstartclin();
  // int get_immu_ss();
  // int get_clin_ss();
  // int get_immu_final();
  // int get_clin_final();
  // int get_startclin_n();
  // 
  // int is_v_samp_stopped();
  // int is_immu_fut();
  // int is_clin_fut();
  // int is_clin_es();
  // int is_inconclusive();
  // 
  // 
  // void set_immu_stopv();
  // void set_immu_fut();
  // void set_clin_fut();
  // void set_clin_es();
  // void set_inconclusive();
  // void set_immu_ss(int n);
  // void set_clin_ss(int n);
  // void set_immu_final_win(bool won);
  // void set_clin_final_win(bool won);
  // 
  // int is_v_samp_stopped();
  // int is_immu_fut();
  // int is_clin_fut();
  // int is_clin_es();
  // int is_inconclusive();
  // 
  // int get_nmaxsero();
  // int get_nstartclin();
  // int get_immu_ss();
  // int get_clin_ss();
  // int get_immu_final();
  // int get_clin_final();
  // int get_startclin_n();
  // 
  // 
  // 
  // 
  // Rcpp::Rcout << "test 1 " << t.get_clin_ss() << std::endl;
  // 
  // t.set_clin_ss(10);
  // 
  // Rcpp::Rcout << "test 2 " << t.get_clin_ss() << std::endl;
  
  // int a = 5;
  
  // return x * a;
  
  return;
}



