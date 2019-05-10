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

  int n = cfg["nstop"];
  arma::mat d = arma::zeros(n, NCOL);
  double tpp = (double)cfg["months_per_person"];

  for(int i = 0; i < n; i++){
    
    d(i, COL_ID) = i+1;
    d(i, COL_TRT) = ((i-1)%2 == 0) ? 0 : 1;
    // poisson proc
    if(i == 0) {
      d(i, COL_ACCRT) = R::rexp(1/(double)cfg["people_per_interim_period"])  ;  
    } else {
      d(i, COL_ACCRT) = d(i-1, COL_ACCRT) + R::rexp(1/(double)cfg["people_per_interim_period"])  ;  
    }
    
    d(i, COL_AGE) = R::runif((double)cfg["age_months_lwr"], (double)cfg["age_months_upr"]);
    
    d(i, COL_SEROT2) = R::rbinom(1, cfg["baselineprobsero"]);
    d(i, COL_SEROT3) = d(i, COL_SEROT2);
    d(i, COL_PROBT3) = d(i, COL_TRT) * (double)cfg["deltaserot3"];
    
    if(d(i, COL_SEROT2) == 0 && d(i, COL_TRT) == 1){
      d(i, COL_SEROT3) = R::rbinom(1, d(i, COL_PROBT3));
    }
    
    // tte - the paramaterisation of rexp uses SCALE NOTE RATE!!!!!!!!!!!
    // event time is the time from randomisation (not birth) at which first
    // medical presentation occurs
    if(d(i, COL_TRT) == 0){
      d(i, COL_EVTT) = R::rexp(1/(double)cfg["b0tte"])  ;
    } else {
      double beta = (double)cfg["b0tte"] + (double)cfg["b1tte"];
      d(i, COL_EVTT) = R::rexp(1/beta)  ;
    }
    
  }
  
  d.col(COL_CEN) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_OBST) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_REASON) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_IMPUTE) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_REFTIME) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  
  return d;
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



