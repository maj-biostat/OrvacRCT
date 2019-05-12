#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include <algorithm>
#include "simulator.h"
#include "trial.h"


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
Rcpp::List simulate_trial(int idxsim, Rcpp::List& cfg, bool rtn_trial_dat) {


  Rcpp::Rcout << " STARTING TRIAL ID " << idxsim << std::endl;
  
  Trial t(cfg, idxsim);
  
  // Rcpp::Rcout << " get_thresh_pp_fut " << t.get_thresh_pp_fut() << std::endl;
  // Rcpp::Rcout << " get_thresh_pp_es  " << t.get_thresh_pp_es() << std::endl;
  // Rcpp::Rcout << " get_thresh_p_sup  " << t.get_thresh_p_sup() << std::endl;
  t.print_cfg();

  t.run_interims();
  
  // Rcpp::List res_fin = run_fin(d, cfg, t, idxsim);
  // Rcpp::Rcout << " back in simulate trial t.get_c_ppos_n() " 
  //                 << t.get_c_ppos_n() << std::endl;
  
  Rcpp::List ret;
  ret["trial"] = 1;
  // ret["trial"] = t.as_list();
  
  return ret;
}


// the following are just hacks to support access to Trial class

// // [[Rcpp::export]]
// void trial_print(const Rcpp::List& cfg){
//   Trial t((int)cfg["n_max_sero"], (int)cfg["n_start_clin"]);
//   t.print_state(1);
//   return;
// }
// 
// // [[Rcpp::export]]
// Rcpp::List trial_as_list(const Rcpp::List& cfg){
//   Trial t((int)cfg["n_max_sero"], (int)cfg["n_start_clin"]);
//   Rcpp::List ret = t.as_list();
//   return ret;
// }

