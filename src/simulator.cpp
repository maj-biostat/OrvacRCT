#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include <algorithm>
#include "simulator.h"
#include "trial.h"



// [[Rcpp::export]]
arma::mat get_interims(Rcpp::List& cfg){
  
  Trial t(cfg, 1);
  
  Rcpp::Rcout << " STARTING TRIAL ID " << 1 << std::endl;
  
  return t.get_interims();
}
   
// [[Rcpp::export]]
Rcpp::List test_set_state(Rcpp::List& cfg,
  int n_target, double ref_time, bool dofu, int cur_intrm){
  
  Trial t(cfg, 1);
  Rcpp::List l ;
  if(dofu == 1){
    Rcpp::Rcout << " doing with fu " << std::endl;
    l = t.clin_set_state(n_target, ref_time);
  } else {
    Rcpp::Rcout << " doing without fu " << std::endl;
    
    t.set_curr_intrm_idx(cur_intrm);
    l = t.clin_set_state();
  }
  
  return l;
}


   
 
// [[Rcpp::export]]
Rcpp::List simulate_trial(int idxsim, Rcpp::List& cfg, bool rtn_trial_dat) {


  Rcpp::Rcout << " STARTING TRIAL ID " << idxsim << std::endl;
  
  Trial t(cfg, idxsim);

  t.print_cfg();
  t.run_interims();
  t.run_final();
  t.print_state();
  
  
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

