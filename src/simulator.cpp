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
Rcpp::List test_set_state(Rcpp::List& cfg, bool dofu, int cur_intrm_idx, 
                          int n_target, double ref_time){
  
  Trial t(cfg, 1);

  Rcpp::List l ;
  arma::mat dcop;
  arma::mat icop;
  
  if(dofu == 0){
    
    icop = t.get_interims();
    int n_targ = icop(cur_intrm_idx, INT_N);
    double ref_t = icop(cur_intrm_idx, INT_T_END);
    
    Rcpp::Rcout << " doing without fu n_target " << n_targ << " ref_time " << ref_t << " intrm_idx " << cur_intrm_idx << std::endl;
    
    t.set_curr_intrm_idx(cur_intrm_idx);
    l = t.clin_censoring();
    dcop = t.get_data();

    //arma::uvec uimpute = arma::find(d.col(COL_IMPUTE) == 1);
    
  } else {
    
    icop = t.get_interims();

    Rcpp::Rcout << " doing with fu n_target " << n_target << " ref_time " << ref_time << " intrm_idx " << cur_intrm_idx << std::endl;

    t.set_curr_intrm_idx(cur_intrm_idx);
    l = t.clin_censoring(n_target, ref_time);
    dcop = t.get_data();
  }
  
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("state") = l,
                     Rcpp::Named("dat") = dcop,
                     Rcpp::Named("intrm") = icop);
  
  return ret;
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

