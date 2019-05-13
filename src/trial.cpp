#include <RcppArmadillo.h>


#include <cmath>
#include <algorithm>
#include "simulator.h"
#include "trial.h"


Trial::Trial(Rcpp::List& trial_cfg, int idxsim){
  
  sim_id = idxsim;
  cfg = trial_cfg;
  d = init_trial_dat();
  intrms = init_interims();
  
  print_interims();
}


arma::mat Trial::init_trial_dat(){

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
  
  d.col(COL_SEROIMPUTE) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_CEN) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_OBST) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_REASON) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_IMPUTE) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_REFTIME) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  
  return d;
}

arma::mat Trial::init_interims(){
  
  // get_interims
  // interims occur at index INT_I_END at time INT_T_END
  // an interim uses kids from 0 to INT_I_END with 
  // kids at index INT_I_START to INT_I_END being introduced
  // for this interim.
  
  int n = (int)cfg["n_stop"];
  int n_start = (int)cfg["n_start"] - 1;
  
  // I resize later
  arma::mat intrms = arma::zeros(300, 5);
  
  intrms(0, INT_I_START) = 0;
  intrms(0, INT_I_END) = n_start;                // 69
  intrms(0, INT_T_START) = 0;
  intrms(0, INT_T_END) = d(n_start, COL_ACCRT);  // accrual for 69
  intrms(0, INT_N) = intrms(0, INT_I_END) +1;   

  intrms(1, INT_I_START) = intrms(0, INT_I_END) + 1;    // 70
  intrms(1, INT_T_START) = intrms(0, INT_T_END);        // accrual for 69
  
  int k = 2;
  int i = 0;
  
  for(i = intrms(1, INT_I_START) + 1; i < n; i++){
    
    // every three months or next 50 enrollees, whichever sooner
    // start the clock from the accrual of the previous kid
    
    // Rcpp::Rcout << " comparing d(" << i << ", COL_ACCRT) "
    //   << d(i, COL_ACCRT) << " with intrms(" 
    //   << k - 1 << ", INT_T_START) " 
    //   << intrms(k-1, INT_T_START) << " diff is "
    //   << d(i, COL_ACCRT) - intrms(k-1, INT_T_START) << " n "
    //   << d(i, COL_ID) - intrms(k-1, INT_I_START)
    //   << std::endl;
    
    // first we figure out where the previous interim should end.
    // once we have that number we know where the current one should
    // start.
    
    if(d(i, COL_ACCRT) - intrms(k-1, INT_T_START) >= (double)cfg["int_mnth"]  |  
      d(i, COL_ID) - intrms(k-1, INT_I_START)  >= (int)cfg["int_n"]-1) {

      intrms(k-1, INT_I_END) = i;
      // start the clock from the accrual of the previous kid
      intrms(k-1, INT_T_END) = d(i, COL_ACCRT);
      intrms(k-1, INT_N) = intrms(k-1, INT_I_END) + 1;   
      
      intrms(k, INT_I_START) = intrms(k-1, INT_I_END) + 1;
      // start the clock from the accrual of the previous kid
      intrms(k, INT_T_START) = intrms(k-1, INT_T_END);

      k++;
    }
  }
  
  
  // now we trim down the interims so that we have a complete set
  // that does not include the final.
  int z = 0;
  int ilast = (int)intrms(k, INT_I_END);

  while(ilast == 0 | ilast == n-1){
    
    k = k - 1;
    ilast = (int)intrms(k, INT_I_END);
    
    z++;
    if(z > 10) break;
  }
  
  intrms.resize(k+1, 5);
  
  return intrms;
  
}

arma::mat Trial::get_interims(){
  
  arma::mat intrm_copy = arma::mat(intrms);
  return intrm_copy;
}

void Trial::run_final(){
  
  // we are now doing the final analysis
  is_final = 1;
  
  immu_fin();
  clin_fin();
  
}

void Trial::clin_fin(){
  
  double a = (double)cfg["prior_gamma_a"];
  double b = (double)cfg["prior_gamma_b"];
  int n_target = 0;
  
  // if no early stopping then look at max size
  if(is_clin_fut() == 0 && is_clin_es() == 0 && is_immu_fut() == 0){
    n_target = (int)cfg["n_stop"];
    set_enrolled_ss(n_target);
    INFO(Rcpp::Rcout, get_sim_id(), "final clin no early stop n = " << n_target);
  } else {
    n_target = get_enrolled_ss();
    INFO(Rcpp::Rcout, get_sim_id(), "final clin early stop n = " << n_target);
  }
  
  set_clin_ss(n_target);
  c_suf_fin = clin_censoring(n_target, 0);
  
  int tot_imp = 0;
  for(int i = 0; i < n_enrolled; i++){
    tot_imp = tot_imp + d(i, COL_IMPUTE);
    
  }
  
  INFO(Rcpp::Rcout, get_sim_id(), "final clin imputing " << tot_imp << " (should be 0)");
  
  arma::mat m = arma::zeros((int)cfg["post_draw"] , 3);  
  
  for(int i = 0; i < (int)cfg["post_draw"]; i++){
    
    m(i, COL_LAMB0) = R::rgamma(a + (double)c_suf_fin["n_evnt_0"], 
      1/(b + (double)c_suf_fin["tot_obst_0"]));
    m(i, COL_LAMB1) = R::rgamma(a + (double)c_suf_fin["n_evnt_1"], 
      1/(b + (double)c_suf_fin["tot_obst_1"]));
    m(i, COL_RATIO) = m(i, COL_LAMB1) / m(i, COL_LAMB0);
    
  }

  arma::uvec ugt1 = arma::find(m.col(COL_RATIO) < 1);
  c_p_fin =  (double)ugt1.n_elem / (double)cfg["post_draw"];
  
  if(c_p_fin > (double)cfg["thresh_p_sup"]){
    set_clin_fin_decision(1);
  } else if(c_p_fin < (double)cfg["thresh_p_fut"]){
    set_clin_fin_decision(-1);
  } else {
    set_clin_fin_decision(0);
  }
  
  INFO(Rcpp::Rcout, get_sim_id(), "final clin c_p_fin " << c_p_fin
          << " vs sup thresh " << (double)cfg["thresh_p_sup"]
          << " decision = " << get_clin_fin_decision() 
          << " (-1 is fut, 0 is nd, 1 is won)");
  
}








void Trial::immu_fin(){
  
  double a = (double)cfg["prior_beta_a"];
  double b = (double)cfg["prior_beta_b"];
  int n_target = 0;
  
  // if no early stopping then look at max size
  if(is_immu_fut() == 0 && is_v_samp_stopped() == 0){
    n_target = (int)cfg["n_max_sero"];
    INFO(Rcpp::Rcout, get_sim_id(), "final immu no early stop n = " << n_target);
  } else {
    n_target = get_immu_ss();
    INFO(Rcpp::Rcout, get_sim_id(), "final immu early stop n = " << n_target);
  }
  
  set_immu_ss(n_target);
  i_suf_fin = immu_observed(n_target);
  
  int tot_imp = 0;
  for(int i = 0; i < n_target; i++){
    tot_imp = tot_imp + d(i, COL_SEROIMPUTE);
    
  }
  
  INFO(Rcpp::Rcout, get_sim_id(), "final immu imputing " << tot_imp << " (should be 0)");
  
  arma::mat m = arma::zeros((int)cfg["post_draw"] , 3);  
  
  for(int i = 0; i < (int)cfg["post_draw"]; i++){
    
    m(i, COL_THETA0) = R::rbeta(a + (double)i_suf_fin["n_evnt_0"], 
      b + (double)i_suf_fin["n_tot_0"] - (double)i_suf_fin["n_evnt_0"]);
    m(i, COL_THETA1) = R::rbeta(a + (double)i_suf_fin["n_evnt_1"], 
      b + (double)i_suf_fin["n_tot_1"] - (double)i_suf_fin["n_evnt_1"]);
    m(i, COL_DELTA) = m(i, COL_THETA1) - m(i, COL_THETA0);
    
  }
  
  arma::uvec ugt0 = arma::find(m.col(COL_DELTA) > 0);
  i_p_fin =  (double)ugt0.n_elem / (double)cfg["post_draw"];
  
  if(i_p_fin > (double)cfg["thresh_p_sup"]){
    set_immu_fin_decision(1);
  } else if(i_p_fin < (double)cfg["thresh_p_fut"]){
    set_immu_fin_decision(-1);
  } else {
    set_immu_fin_decision(0);
  }
  
  INFO(Rcpp::Rcout, get_sim_id(), "final immu i_p_fin " << i_p_fin
    << " vs. sup thresh " << (double)cfg["thresh_p_sup"]
    << " decision = " << get_immu_fin_decision() 
    << " (-1 is fut, 0 is nd, 1 is won)");
  
}













void Trial::run_interims(){

  int i = 0;
  for(i = 0; i < intrms.n_rows ;  i++){
    
    // keep track of how many currently enrolled
    this->idx_cur_intrm = i;
    this->n_enrolled = intrms(i, INT_N);
    
    // set clin_ss here but not immu_ss as we do not want 
    // to continue updating immu_ss after 250.
    this->clin_ss = intrms(i, INT_N);
    
    INFO(Rcpp::Rcout, get_sim_id(), "STARTING: intrm " << idx_cur_intrm
      << " n_enrolled " << n_enrolled
      << " do immu " << do_immu()
      << " do clin " << do_clin()
      << " thresh_pp_fut " << get_thresh_pp_fut()
      << " thresh_pp_es " << get_thresh_pp_es()
      << " thresh_p_sup " << get_thresh_p_sup() 
      << " post_draw " << (int)cfg["post_draw"]);
    
    
    if(do_immu()){
      
      // setting immu_ss here ensures that we do not overshoot the
      // max ss of 250.
      this->immu_ss = intrms(i, INT_N);
      // set idx_cur_intrm prior to use
      immu_interim();
      
      INFO(Rcpp::Rcout, get_sim_id(), "intrm " << idx_cur_intrm 
        << " immu i_p_n = " << i_p_n 
        << ", i_ppos_n = " << i_ppos_n 
        << " vs. pp_es_thresh = " << (double)cfg["thresh_pp_es"]
        << ", i_ppos_max = " << i_ppos_max
        << " vs. pp_fut_thresh = " << (double)cfg["thresh_pp_fut"]);
      
      if(i_ppos_max < (double)cfg["thresh_pp_fut"]){
        INFO(Rcpp::Rcout, get_sim_id(), "intrm " << idx_cur_intrm 
          << " immu hit futile - stop now, ppmax " 
          << i_ppos_max
          << " fut_thresh " << (double)cfg["thresh_pp_fut"]);
        set_immu_fut();
        break;
      }
      
      if (i_ppos_n > (double)cfg["thresh_pp_es"] && !is_immu_fut()){
        INFO(Rcpp::Rcout, get_sim_id(), "intrm " << idx_cur_intrm 
          << " immu es - stop v samp now, ppn " 
          << i_ppos_n
          << " es thresh " << (double)cfg["thresh_pp_es"]);
        set_immu_stopv();
        // do not break yet
      }
    }
    
    if(do_clin()){
      
      // set idx_cur_intrm prior to use
      clin_interim();
  
      INFO(Rcpp::Rcout, get_sim_id(), "intrm " << idx_cur_intrm 
        << " clin c_p_n " << c_p_n 
        << "  c_ppos_n = " << c_ppos_n 
        << " vs. pp_es_thresh = " << (double)cfg["thresh_pp_es"]
        << ", c_ppos_max = " << c_ppos_max
        << " vs. pp_fut_thresh = " << (double)cfg["thresh_pp_fut"]);
      
      if(c_ppos_max < (double)cfg["thresh_pp_fut"]){
        INFO(Rcpp::Rcout, get_sim_id(), "intrm " << idx_cur_intrm 
          << " clin hit futile - stop now, ppmax " 
          << c_ppos_max
          << " fut_thresh " << (double)cfg["thresh_pp_fut"]);
        set_clin_fut();
        break;
      }
      
      if (c_ppos_n > (double)cfg["thresh_pp_es"] && !is_clin_fut()){
        INFO(Rcpp::Rcout, get_sim_id(), "intrm " << idx_cur_intrm 
          << " clin hit es - stop now, ppn " 
          << c_ppos_n
          << " es_thresh " << (double)cfg["thresh_pp_es"]);
        set_clin_es();
        break;
      }
    }
  }
  
  if(is_v_samp_stopped() == 0 && is_immu_fut() == 0 &&
     is_clin_fut() == 0 && is_clin_es() == 0){
    
    INFO(Rcpp::Rcout, get_sim_id(), "intrms are inconclusive - no stopping rules hit ");
    set_inconclusive();
  }

}

void Trial::immu_interim(){
  

  // stores samples from posterior, 
  arma::mat m = arma::zeros((int)cfg["post_draw"] , 3);
  arma::mat m_pp_int = arma::zeros((int)cfg["post_draw"] , 3);
  arma::mat m_pp_max = arma::zeros((int)cfg["post_draw"] , 3);
  
  arma::uvec ugt0;
  int int_win = 0;
  int max_win = 0;
  
  double a = (double)cfg["prior_beta_a"];
  double b = (double)cfg["prior_beta_b"];

  i_suf_intrm = immu_observed();
  
  // subjs that require imputation
  arma::uvec uimpute = arma::find(d.col(COL_SEROIMPUTE) == 1);
  
  INFO(Rcpp::Rcout, get_sim_id(), "intrm " << idx_cur_intrm 
    << " immu using up to n = " << (int)intrms(idx_cur_intrm, INT_N) 
    << " n_evnt_0 " << (int)i_suf_intrm["n_evnt_0"]
    << " n_evnt_1 " << (int)i_suf_intrm["n_evnt_1"]
    << " n_tot_0 " << (int)i_suf_intrm["n_tot_0"]
    << " n_tot_0 " << (int)i_suf_intrm["n_tot_1"]
    << " unobserved " << (int)uimpute.n_elem
  );

  // keep a copy of the original state
  arma::mat d_orig = arma::zeros((int)cfg["n_stop"], 1);
  d_orig.col(0) = arma::vec(d.col(COL_SEROT3));

  for(int i = 0; i < (int)cfg["post_draw"]; i++){
    
    // draw from the posterior at current time.
    m(i, COL_THETA0) = R::rbeta(a + (double)i_suf_intrm["n_evnt_0"], 
      b + (double)i_suf_intrm["n_tot_0"] - (double)i_suf_intrm["n_evnt_0"]);
    m(i, COL_THETA1) = R::rbeta(a + (double)i_suf_intrm["n_evnt_1"], 
      b + (double)i_suf_intrm["n_tot_1"] - (double)i_suf_intrm["n_evnt_1"]);
    m(i, COL_DELTA) = m(i, COL_THETA1) - m(i, COL_THETA0);
    
    // Rcpp::Rcout << " post1 m(i, COL_THETA0) " << m(i, COL_THETA0) 
    //             << " m(i, COL_THETA1) " << m(i, COL_THETA1) << std::endl;
    
    // impute the ones that we do not yet have events for up to the 
    // current time
    for(int j = 0; j < (int)uimpute.n_elem; j++){
      int sub_idx = uimpute(j);
      if(d(sub_idx, COL_TRT) == 0){
        d(sub_idx, COL_SEROT3) = R::rbinom(1, m(i, COL_THETA0));
      } else {
        d(sub_idx, COL_SEROT3) = R::rbinom(1, m(i, COL_THETA1));
      }
    }
    
    // compute pred prob at the time of interim
    i_suf_intrm_fu = immu_observed((int)intrms(idx_cur_intrm, INT_N));
    
    // Rcpp::Rcout << " with impute to n n_evnt_0 " << (int)i_suf_intrm["n_evnt_0"]
    //             << " n_evnt_1 " << (int)i_suf_intrm["n_evnt_1"]
    //             << " n_tot_0 " << (int)i_suf_intrm["n_tot_0"]
    //             << " n_tot_0 " << (int)i_suf_intrm["n_tot_1"]
    //             << std::endl;
    
    for(int j = 0; j < (int)cfg["post_draw"]; j++){
      m_pp_int(j, COL_THETA0) = R::rbeta(a + (double)i_suf_intrm_fu["n_evnt_0"], 
        b + (double)i_suf_intrm_fu["n_tot_0"] - (double)i_suf_intrm_fu["n_evnt_0"]);
      m_pp_int(j, COL_THETA1) = R::rbeta(a + (double)i_suf_intrm_fu["n_evnt_1"], 
        b + (double)i_suf_intrm_fu["n_tot_1"] - (double)i_suf_intrm_fu["n_evnt_1"]);
      m_pp_int(j, COL_DELTA) = m_pp_int(j, COL_THETA1) - m_pp_int(j, COL_THETA0);
    }
    
    // Rcpp::Rcout << " post2 m_pp_int(999, COL_THETA0) " << m_pp_int(999, COL_THETA0) 
    //             << " m_pp_int(999, COL_THETA1) " << m_pp_int(999, COL_THETA1) << std::endl;

    // empirical posterior probability that ratio_lamb > 1
    ugt0 = arma::find(m_pp_int.col(COL_DELTA) > 0);
    if((double)ugt0.n_elem / (double)cfg["post_draw"] > (double)cfg["thresh_p_sup"]){
      int_win++;
    }
    
    // impute the ones that we do not yet have events for up to the 
    // max samp size of 250
    for(int k = (int)intrms(idx_cur_intrm, INT_I_START); k < (int)cfg["n_max_sero"]; k++){
      if(d(k, COL_TRT) == 0){
        d(k, COL_SEROT3) = R::rbinom(1, m(i, COL_THETA0));
      } else {
        d(k, COL_SEROT3) = R::rbinom(1, m(i, COL_THETA1));
      }
    }
    
    // compute pred prob assuming we go up to nmax
    // what does the posterior at max sample size say?
    i_suf_max_fu = immu_observed((int)cfg["n_max_sero"]);
    
    
    // Rcpp::Rcout << " at max n_evnt_0 " << (int)i_suf_max_fu["n_evnt_0"]
    //             << " n_evnt_1 " << (int)i_suf_max_fu["n_evnt_1"]
    //             << " n_tot_0 " << (int)i_suf_max_fu["n_tot_0"]
    //             << " n_tot_0 " << (int)i_suf_max_fu["n_tot_1"]
    //             << std::endl;
    
    
    for(int j = 0; j < (int)cfg["post_draw"]; j++){
      m_pp_max(j, COL_THETA0) = R::rbeta(a + (double)i_suf_max_fu["n_evnt_0"], 
               b + (double)i_suf_max_fu["n_tot_0"] - (double)i_suf_max_fu["n_evnt_0"]);
      m_pp_max(j, COL_THETA1) = R::rbeta(a + (double)i_suf_max_fu["n_evnt_1"], 
               b + (double)i_suf_max_fu["n_tot_1"] - (double)i_suf_max_fu["n_evnt_1"]);
      m_pp_max(j, COL_DELTA) = m_pp_max(j, COL_THETA1) - m_pp_max(j, COL_THETA0);
    }
    
    // Rcpp::Rcout << " post3 m_pp_max(999, COL_THETA0) " << m_pp_max(999, COL_THETA0) 
    //             << " m_pp_max(999, COL_THETA1) " << m_pp_max(999, COL_THETA1) << std::endl;    

    
    // empirical posterior probability that ratio_lamb > 1
    ugt0 = arma::find(m_pp_max.col(COL_DELTA) > 0);
    if((double)ugt0.n_elem / (double)cfg["post_draw"] > (double)cfg["thresh_p_sup"]){
      max_win++;
    }
    
    // reset to original state ready for the next posterior draw
    d.col(COL_SEROT3) = arma::vec(d_orig.col(0));
    
  }
  
  // reset so have clean state for next interim
  d.col(COL_SEROIMPUTE) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  
  // Rcpp::Rcout << " int_win      " << max_win << std::endl;
  // Rcpp::Rcout << " max_win      " << max_win << std::endl;

  ugt0 = arma::find(m.col(COL_DELTA) > 0);
  this->i_p_n = (double)ugt0.n_elem/(double)cfg["post_draw"];
  this->i_ppos_n = (double)int_win/(double)cfg["post_draw"];
  this->i_ppos_max = (double)max_win/(double)cfg["post_draw"];

  INFO(Rcpp::Rcout, get_sim_id(), "intrm " << idx_cur_intrm 
    << " immu int_win " << int_win << " max_win " << max_win 
    << " out of " << (double)cfg["post_draw"]);
  
  return;
}

// determine:
// 1. how many observations do we have in both arms at the current time 
//    when we factor in that there is a delay in information for the 
//    seroconversion status
// 2. how many successes in each arm
Rcpp::List Trial::immu_observed(){
  
  int n_evnt_0 = 0;
  int n_evnt_1 = 0;
  int n_tot_0 = 0;
  int n_tot_1 = 0;
  
  int n_target = intrms(idx_cur_intrm, INT_N);
  double ref_time = intrms(idx_cur_intrm, INT_T_END);
  
  for(int sub_idx = 0; sub_idx < n_target; sub_idx++){
    
    if(d(sub_idx, COL_ACCRT) + (double)cfg["sero_info_delay"] > ref_time){
      
      d(sub_idx, COL_SEROIMPUTE) = 1;

      continue;
    }
    
    d(sub_idx, COL_SEROIMPUTE) = 0;
    
    if(d(sub_idx, COL_TRT) == 0){
      n_tot_0++;
      if(d(sub_idx, COL_SEROT3) == 1){
        n_evnt_0++;
      }
      
    } else {
      n_tot_1++;
      if(d(sub_idx, COL_SEROT3) == 1){
        n_evnt_1++;
      }
    }
    
  }
  
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("n_evnt_0") = n_evnt_0,
                                      Rcpp::Named("n_tot_0") = n_tot_0,
                                      Rcpp::Named("n_evnt_1") = n_evnt_1,
                                      Rcpp::Named("n_tot_1") = n_tot_1);
  
  return ret;
}
  
  
Rcpp::List Trial::immu_observed(int n_target){
  
  int n_evnt_0 = 0;
  int n_evnt_1 = 0;
  int n_tot_0 = 0;
  int n_tot_1 = 0;
  
  for(int sub_idx = 0; sub_idx < n_target; sub_idx++){
    
    d(sub_idx, COL_SEROIMPUTE) = 0;
    
    if(d(sub_idx, COL_TRT) == 0){
      n_tot_0++;
      if(d(sub_idx, COL_SEROT3) == 1){
        n_evnt_0++;
      }
      
    } else {
      n_tot_1++;
      if(d(sub_idx, COL_SEROT3) == 1){
        n_evnt_1++;
      }
    }
    
  }
  
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("n_evnt_0") = n_evnt_0,
                                      Rcpp::Named("n_tot_0") = n_tot_0,
                                      Rcpp::Named("n_evnt_1") = n_evnt_1,
                                      Rcpp::Named("n_tot_1") = n_tot_1);
  
  return ret;
}  
  
  
  
void Trial::clin_interim(){

  // stores samples from posterior, 
  arma::mat m = arma::zeros((int)cfg["post_draw"] , 3);
  arma::mat m_pp_int = arma::zeros((int)cfg["post_draw"] , 3);
  arma::mat m_pp_max = arma::zeros((int)cfg["post_draw"] , 3);

  arma::uvec ugt1;
  int int_win = 0;
  int max_win = 0;
  
  double a = (double)cfg["prior_gamma_a"];
  double b = (double)cfg["prior_gamma_b"];
  
  c_suf_intrm = clin_censoring();
  // subjs that require imputation
  arma::uvec uimpute = arma::find(d.col(COL_IMPUTE) == 1);
  
  INFO(Rcpp::Rcout, get_sim_id(), "intrm " << idx_cur_intrm 
    << " clin using up to n = " << (int)intrms(idx_cur_intrm, INT_N) 
    << " n_evnt_0 " << (int)c_suf_intrm["n_evnt_0"]
    << " n_evnt_1 " << (int)c_suf_intrm["n_evnt_1"]
    << " tot_obst_0 " << (double)c_suf_intrm["tot_obst_0"]
    << " tot_obst_1 " << (double)c_suf_intrm["tot_obst_1"]
    << " unobserved " << (int)uimpute.n_elem
  );
  
  
  // keep a copy of the original state
  arma::mat d_orig = arma::zeros((int)cfg["n_stop"], 6);
  d_orig.col(0) = arma::vec(d.col(COL_EVTT));
  d_orig.col(1) = arma::vec(d.col(COL_CEN));
  d_orig.col(2) = arma::vec(d.col(COL_OBST));
  d_orig.col(3) = arma::vec(d.col(COL_REASON));
  d_orig.col(4) = arma::vec(d.col(COL_IMPUTE));
  d_orig.col(5) = arma::vec(d.col(COL_REFTIME));
  
  
  // for i in postdraws do posterior predictive trials
  // 1. for the interim (if we are at less than 50 per qtr)
  // 2. for the max sample size
  for(int i = 0; i < (int)cfg["post_draw"]; i++){

    // compute the posterior based on the __observed__ data to the time of the interim
    // take single draw
    m(i, COL_LAMB0) = R::rgamma(a + (double)c_suf_intrm["n_evnt_0"], 
      1/(b + (double)c_suf_intrm["tot_obst_0"]));
    m(i, COL_LAMB1) = R::rgamma(a + (double)c_suf_intrm["n_evnt_1"], 
      1/(b + (double)c_suf_intrm["tot_obst_1"]));
    m(i, COL_RATIO) = m(i, COL_LAMB1) / m(i, COL_LAMB0);
    
    //Rcpp::Rcout << " m(i, COL_RATIO)      " << m(i, COL_RATIO) << std::endl;

    // use memoryless prop of exponential and impute enrolled kids that have not
    // yet had event.
    for(int j = 0; j < (int)uimpute.n_elem; j++){
      int sub_idx = uimpute(j);
      if(d(sub_idx, COL_TRT) == 0){

        d(sub_idx, COL_EVTT) = d(sub_idx, COL_OBST) + R::rexp(1/m(i, COL_LAMB0))  ;

      } else {

        d(sub_idx, COL_EVTT) = d(sub_idx, COL_OBST) + R::rexp(1/m(i, COL_LAMB1))  ;

      }
    }
  
    // update view of the sufficent stats using enrolled
    // kids that have all now been given an event time
    c_suf_intrm_fu = clin_censoring((int)intrms(idx_cur_intrm, INT_N), 
      (double)intrms(idx_cur_intrm, INT_T_END));
    
    for(int j = 0; j < (int)cfg["post_draw"]; j++){
      m_pp_int(j, COL_LAMB0) = R::rgamma(a + (double)c_suf_intrm_fu["n_evnt_0"], 
        1/(b + (double)c_suf_intrm_fu["tot_obst_0"]));
      m_pp_int(j, COL_LAMB1) = R::rgamma(a + (double)c_suf_intrm_fu["n_evnt_1"], 
        1/(b + (double)c_suf_intrm_fu["tot_obst_1"]));
      m_pp_int(j, COL_RATIO) = m_pp_int(j, COL_LAMB1) / m_pp_int(j, COL_LAMB0);
    }
    // empirical posterior probability that ratio_lamb < 1
    ugt1 = arma::find(m_pp_int.col(COL_RATIO) < 1);
    if((double)ugt1.n_elem / (double)cfg["post_draw"] > (double)cfg["thresh_p_sup"]){
      int_win++;
    }

    // impute the remaining kids up to max sample size
    for(int k = (int)intrms(idx_cur_intrm, INT_I_START); k < d.n_rows; k++){
      if(d(k, COL_TRT) == 0){
        d(k, COL_EVTT) = R::rexp(1/m(i, COL_LAMB0))  ;
      } else {
        d(k, COL_EVTT) = R::rexp(1/m(i, COL_LAMB1))  ;
      }
    }

    // set the state up to the max sample size at time of the final analysis
    // what we put in for ref_time is irrelevant as it will be updated
    // inside clin_censoring
    c_suf_max_fu = clin_censoring((int)d.n_rows, 0);

    // what does the posterior at max sample size say?
    for(int j = 0; j < (int)cfg["post_draw"]; j++){
      m_pp_max(j, COL_LAMB0) = R::rgamma(a + (double)c_suf_max_fu["n_evnt_0"], 
        1/(b + (double)c_suf_max_fu["tot_obst_0"]));
      m_pp_max(j, COL_LAMB1) = R::rgamma(a + (double)c_suf_max_fu["n_evnt_1"], 
        1/(b + (double)c_suf_max_fu["tot_obst_1"]));
      m_pp_max(j, COL_RATIO) = m_pp_max(j, COL_LAMB1) / m_pp_max(j, COL_LAMB0);
    }
    // empirical posterior probability that ratio_lamb < 1
    ugt1 = arma::find(m_pp_max.col(COL_RATIO) < 1);
    if((double)ugt1.n_elem / (double)cfg["post_draw"] > (double)cfg["thresh_p_sup"]){
      max_win++;
    }

    // reset to original state ready for the next posterior draw
    d.col(COL_EVTT) = arma::vec(d_orig.col(0));
    d.col(COL_CEN) = arma::vec(d_orig.col(1));
    d.col(COL_OBST) = arma::vec(d_orig.col(2));
    d.col(COL_REASON) = arma::vec(d_orig.col(3));
    d.col(COL_IMPUTE) = arma::vec(d_orig.col(4));
    d.col(COL_REFTIME) = arma::vec(d_orig.col(5));

    
  }
  
  
  // reset so have clean state for next interim
  d.col(COL_CEN) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_OBST) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_REASON) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_IMPUTE) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_REFTIME) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  
  
  INFO(Rcpp::Rcout, get_sim_id(), "intrm " << idx_cur_intrm 
    << " clin int_win " << int_win << " max_win " << max_win 
    << " out of " << (double)cfg["post_draw"]);

  ugt1 = arma::find(m.col(COL_RATIO) < 1);
  
  this->c_p_n = (double)ugt1.n_elem/(double)cfg["post_draw"];
  this->c_ppos_n = (double)int_win/(double)cfg["post_draw"];
  this->c_ppos_max = (double)max_win/(double)cfg["post_draw"];

  return;
}




  
// set censoring status up to participant.
// to be used purely for reference time being the current interim
Rcpp::List Trial::clin_censoring(){

  int n_evnt_0 = 0;
  int n_evnt_1 = 0;
  double tot_obst_0 = 0;
  double tot_obst_1 = 0;
  
  int n_target = intrms(idx_cur_intrm, INT_N);
  double ref_time = intrms(idx_cur_intrm, INT_T_END);
  
  // set censoring and event times up to the current interim
  // in this instance we do no imputation - this is purely for the
  // posterior at the current interim (no predictive probs.)
  for(int sub_idx = 0; sub_idx < n_target; sub_idx++){
    
    DBG(Rcpp::Rcout, "STARTING censoring state for sub_idx " << sub_idx);
    
    d(sub_idx, COL_REFTIME) = ref_time;
    
    DBG(Rcpp::Rcout, " ref_time     "  << d(sub_idx, COL_REFTIME)); 
    DBG(Rcpp::Rcout, " acc + evt    "  << d(sub_idx, COL_ACCRT) + d(sub_idx, COL_EVTT)); 
    DBG(Rcpp::Rcout, " age + evt    "  << d(sub_idx, COL_AGE) + d(sub_idx, COL_EVTT));
    DBG(Rcpp::Rcout, " ref - acc    "  << d(sub_idx, COL_REFTIME) - d(sub_idx, COL_ACCRT));
    DBG(Rcpp::Rcout, " agemax - age "  << (double)cfg["max_age_fu"] - d(sub_idx, COL_AGE));
    
    clin_censoring_subj(sub_idx);
    
    if(d(sub_idx, COL_CEN) == 0 && d(sub_idx, COL_TRT) == 0) {
      n_evnt_0 += 1;
    } 
    if(d(sub_idx, COL_CEN) == 0 && d(sub_idx, COL_TRT) == 1) {
      n_evnt_1 += 1;
    } 
    
    if(d(sub_idx, COL_TRT) == 0) {
      tot_obst_0 = tot_obst_0 + d(sub_idx, COL_OBST);
    } 
    if(d(sub_idx, COL_TRT) == 1) {
      tot_obst_1 = tot_obst_1 + d(sub_idx, COL_OBST);
    } 

    DBG(Rcpp::Rcout, " trt          "  << d(sub_idx, COL_TRT)); 
    DBG(Rcpp::Rcout, " obst         "  << d(sub_idx, COL_OBST)); 
    DBG(Rcpp::Rcout, " cen          "  << d(sub_idx, COL_CEN));
    DBG(Rcpp::Rcout, " reas         "  << d(sub_idx, COL_REASON));
    DBG(Rcpp::Rcout, " impute       "  << d(sub_idx, COL_IMPUTE));
    

  }
  
  // Rcpp::Rcout << d << std::endl;

  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("n_evnt_0") = n_evnt_0,
                           Rcpp::Named("tot_obst_0") = tot_obst_0,
                           Rcpp::Named("n_evnt_1") = n_evnt_1,
                           Rcpp::Named("tot_obst_1") = tot_obst_1);

  return ret;
}
      

// for final:
// set n_target to n_enrolled (= nmax if no early stopping) and ref_time to zero 
//  => the max fu for each subj will be computed
// 
// for interim:
// set n_target to interim n and ref_time to time of interim
Rcpp::List Trial::clin_censoring(int n_target, double ref_time){

  int n_evnt_0 = 0;
  int n_evnt_1 = 0;
  double tot_obst_0 = 0;
  double tot_obst_1 = 0;
  double fudge = 0.00001;

  // set censoring and event times up to the current interim
  // in this instance we do no imputation - this is purely for the
  // posterior at the current interim (no predictive probs.)
  for(int sub_idx = 0; sub_idx < n_target; sub_idx++){
    
    DBG(Rcpp::Rcout, "STARTING set censoring state for sub_idx " << sub_idx);
    
    if(ref_time < (double)cfg["max_age_fu"] - d(sub_idx, COL_AGE) + d(sub_idx, COL_ACCRT)){
      
      DBG(Rcpp::Rcout, " max fu is after reftime for sub " << sub_idx); 
      d(sub_idx, COL_REFTIME) = (double)cfg["max_age_fu"] - 
        d(sub_idx, COL_AGE) + d(sub_idx, COL_ACCRT) + fudge;
      
      
    } else{
      
      d(sub_idx, COL_REFTIME) = ref_time + fudge;
    }
    
    DBG(Rcpp::Rcout, " ref_time     "  << d(sub_idx, COL_REFTIME)); 
    DBG(Rcpp::Rcout, " acc + evt    "  << d(sub_idx, COL_ACCRT) + d(sub_idx, COL_EVTT)); 
    DBG(Rcpp::Rcout, " age + evt    "  << d(sub_idx, COL_AGE) + d(sub_idx, COL_EVTT));
    DBG(Rcpp::Rcout, " ref - acc    "  << d(sub_idx, COL_REFTIME) - d(sub_idx, COL_ACCRT));
    DBG(Rcpp::Rcout, " agemax - age "  << (double)cfg["max_age_fu"] - d(sub_idx, COL_AGE));
    
    clin_censoring_subj(sub_idx);
    
    if(d(sub_idx, COL_CEN) == 0 && d(sub_idx, COL_TRT) == 0) {
      n_evnt_0 += 1;
    } 
    if(d(sub_idx, COL_CEN) == 0 && d(sub_idx, COL_TRT) == 1) {
      n_evnt_1 += 1;
    } 
      
    if(d(sub_idx, COL_TRT) == 0) {
      tot_obst_0 = tot_obst_0 + d(sub_idx, COL_OBST);
    } 
    if(d(sub_idx, COL_TRT) == 1) {
      tot_obst_1 = tot_obst_1 + d(sub_idx, COL_OBST);
    } 
        
    DBG(Rcpp::Rcout, " trt          "  << d(sub_idx, COL_TRT)); 
    DBG(Rcpp::Rcout, " obst         "  << d(sub_idx, COL_OBST)); 
    DBG(Rcpp::Rcout, " cen          "  << d(sub_idx, COL_CEN));
    DBG(Rcpp::Rcout, " reas         "  << d(sub_idx, COL_REASON));
    DBG(Rcpp::Rcout, " impute       "  << d(sub_idx, COL_IMPUTE));
  }
  
  // Rcpp::Rcout << d << std::endl;

  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("n_evnt_0") = n_evnt_0,
                           Rcpp::Named("tot_obst_0") = tot_obst_0,
                           Rcpp::Named("n_evnt_1") = n_evnt_1,
                           Rcpp::Named("tot_obst_1") = tot_obst_1);
  
  return ret;
}



void Trial::clin_censoring_subj(int sub_idx){
  
  if(d(sub_idx, COL_ACCRT) + d(sub_idx, COL_EVTT) <= d(sub_idx, COL_REFTIME)){
    
    DBG(Rcpp::Rcout, " evtt before reftime " );
    
    if(d(sub_idx, COL_AGE) + d(sub_idx, COL_EVTT) <= (double)cfg["max_age_fu"]){
      // observed event
      DBG(Rcpp::Rcout, " evtt occurred before too old" );
      
      d(sub_idx, COL_CEN) = 0;
      d(sub_idx, COL_OBST) = d(sub_idx, COL_EVTT);
      d(sub_idx, COL_REASON) = 10;
      d(sub_idx, COL_IMPUTE) = 0;
      
    } else {
      
      DBG(Rcpp::Rcout, " evtt before reftime but kid too old" );
      
      if(d(sub_idx, COL_REFTIME) <
        d(sub_idx, COL_ACCRT) - d(sub_idx, COL_AGE) + (double)cfg["max_age_fu"] ){
        
        // censor at mnth - accrual
        DBG(Rcpp::Rcout, " censored - but at ref_time not yet older than max age" );
        
        d(sub_idx, COL_CEN) = 1;
        d(sub_idx, COL_OBST) = d(sub_idx, COL_REFTIME) - d(sub_idx, COL_ACCRT) ;
        d(sub_idx, COL_REASON) = 11;
        d(sub_idx, COL_IMPUTE) = 1;
        
      } else {
        
        DBG(Rcpp::Rcout, " censored - and kid is already older than max age at ref_time" );
        
        d(sub_idx, COL_CEN) = 1;
        d(sub_idx, COL_OBST) = (double)cfg["max_age_fu"] - d(sub_idx, COL_AGE) ;
        d(sub_idx, COL_REASON) = 12;
        d(sub_idx, COL_IMPUTE) = 0;
      }
      
    }
    
  } else {
    
    DBG(Rcpp::Rcout, " evtt after reftime " );
    
    if(d(sub_idx, COL_REFTIME) <
      d(sub_idx, COL_ACCRT) - d(sub_idx, COL_AGE) + (double)cfg["max_age_fu"] ){
      
      // censor at mnth - accrual
      DBG(Rcpp::Rcout, " censored - but at ref_time not yet older than max age" );
      
      d(sub_idx, COL_CEN) = 1;
      d(sub_idx, COL_OBST) = d(sub_idx, COL_REFTIME) - d(sub_idx, COL_ACCRT) ;
      d(sub_idx, COL_REASON) = 21;
      d(sub_idx, COL_IMPUTE) = 1;
      
    } else {
      
      DBG(Rcpp::Rcout, " censored - and at ref_time older than max age" );
      
      d(sub_idx, COL_CEN) = 1;
      d(sub_idx, COL_OBST) = (double)cfg["max_age_fu"] - d(sub_idx, COL_AGE) ;
      d(sub_idx, COL_REASON) = 22;
      d(sub_idx, COL_IMPUTE) = 0;
    }
    
  }
  
  return;
}

bool Trial::do_immu(){
  if(stop_ven_samp == 1){
    return false;
  }
  if(n_enrolled > nmaxsero){
    return false;
  }
  if(stop_immu_fut == 1){
    return false;
  }
  if(stop_clin_fut == 1){
    return false;
  }
  if(stop_clin_es == 1){
    return false;
  }
  return true;
}
bool Trial::do_clin(){
  
  if(n_enrolled < nstartclin){
    return false;
  }
  if(stop_clin_fut == 1){
    return false;
  }
  if(stop_clin_es == 1){
    return false;
  }
  if(stop_immu_fut == 1){
    return false;
  }
  return true;
}

void Trial::set_nmaxsero(int nmaxs){nmaxsero = nmaxs;}
void Trial::set_nstartclin(int nstartc){nstartclin = nstartc;}
void Trial::set_immu_ss(int n){immu_ss = n;}
void Trial::set_clin_ss(int n){clin_ss = n;}
void Trial::set_enrolled_ss(int n){n_enrolled = n;}
void Trial::set_immu_fin_decision(int decision){i_final_decision = decision;}
void Trial::set_clin_fin_decision(int decision){c_final_decision = decision;}

void Trial::set_curr_intrm_idx(int cur_intrm){idx_cur_intrm = cur_intrm;}

void Trial::set_i_ppos_n(double ppos_n){i_ppos_n = ppos_n;}
void Trial::set_i_ppos_max(double ppos_max){i_ppos_max = ppos_max;}
void Trial::set_c_ppos_n(double ppos_n){c_ppos_n = ppos_n;}
void Trial::set_c_ppos_max(double ppos_max){c_ppos_max = ppos_max;}

void Trial::set_immu_stopv(){stop_ven_samp = 1;}
void Trial::set_immu_fut(){stop_immu_fut = 1;}
void Trial::set_clin_fut(){stop_clin_fut = 1;}
void Trial::set_clin_es(){stop_clin_es = 1;}
void Trial::set_inconclusive(){inconclu = 1;}

// getters

arma::mat Trial::get_data(){arma::mat dcopy = arma::mat(d); return dcopy;}

int Trial::get_sim_id(){return sim_id;}
int Trial::get_nmaxsero(){return nmaxsero;}
int Trial::get_nstartclin(){return nstartclin;}
int Trial::get_immu_ss(){return immu_ss;}
int Trial::get_clin_ss(){return clin_ss;}
int Trial::get_enrolled_ss(){return n_enrolled;}
int Trial::get_immu_fin_decision(){return i_final_decision;}
int Trial::get_clin_fin_decision(){return c_final_decision;}

double Trial::get_i_ppos_n(){return i_ppos_n;}
double Trial::get_i_ppos_max(){return i_ppos_max;} 
double Trial::get_c_ppos_n(){return c_ppos_n;}
double Trial::get_c_ppos_max(){return c_ppos_max;} 

double Trial::get_i_p_n(){return i_p_n;}
double Trial::get_c_p_n(){return c_p_n;} 

double Trial::get_i_p_fin(){return i_p_fin;}
double Trial::get_c_p_fin(){return c_p_fin;} 

double Trial::get_thresh_pp_fut(){return (double)cfg["thresh_pp_fut"];} 
double Trial::get_thresh_pp_es(){return (double)cfg["thresh_pp_es"];} 
double Trial::get_thresh_p_sup(){return (double)cfg["thresh_p_sup"];} 

int Trial::get_num_interims(){return intrms.n_rows;}
int Trial::get_interim_n(int i){return intrms(i, INT_N);}
int Trial::get_interim_time(int i){return intrms(i, INT_T_END);}
int Trial::get_curr_intrm_idx(){return idx_cur_intrm;}

int Trial::is_v_samp_stopped(){return stop_ven_samp;}
int Trial::is_immu_fut(){return stop_immu_fut;}
int Trial::is_clin_fut(){return stop_clin_fut;}
int Trial::is_clin_es(){return stop_clin_es;}
int Trial::is_inconclusive(){return inconclu;}

Rcpp::List Trial::as_list(){
  Rcpp::List ret;
  ret["n_max_sero"] = get_nmaxsero();
  ret["n_start_clin"] = get_nstartclin();
  ret["immu_ss"] = get_immu_ss();
  ret["clin_ss"] = get_clin_ss();
  ret["n_enrolled"] = get_enrolled_ss();
  ret["immu_final"] = get_immu_fin_decision();
  ret["clin_final"] = get_clin_fin_decision();
  ret["v_samp_stopped"] = is_v_samp_stopped();
  ret["immu_fut"] = is_immu_fut();
  ret["clin_fut"] = is_clin_fut();
  ret["clin_es"] = is_clin_es();
  ret["inconclusive"] = is_inconclusive();
  
  return ret;
}

void Trial::print_state(){
  
  INFO(Rcpp::Rcout, get_sim_id(), 
    " SUMMARY - Interims: immu stop_v_samp " << is_v_samp_stopped()
    << " immu_fut " << is_immu_fut()
    << " clin_es " << is_clin_es()
    << " clin_fut " << is_clin_fut()
    << " inconclusive " << is_inconclusive() );
    
  INFO(Rcpp::Rcout, get_sim_id(), 
    " SUMMARY - Interims: immu ppos at intrm = " << get_i_ppos_n()
    << " immu ppos at max = " << get_i_ppos_max()
    << " clin ppos at intrm = " << get_c_ppos_n()
    << " clin ppos at max = " << get_c_ppos_max()
  );
    
  INFO(Rcpp::Rcout, get_sim_id(), 
    " SUMMARY - Final analysis: immu_decision " 
    << get_immu_fin_decision()
    << " post prob that delta in seroconv > 0 = " << i_p_fin
    << " clin_decision " 
    << get_clin_fin_decision() 
    << " post prob of hazard ratio < 1 = " << c_p_fin);
  
}

void Trial::print_cfg(){
  
  INFO(Rcpp::Rcout, get_sim_id(), "Trial configuration: " );
  
  int n = cfg.size();
  for(int i = 0; i < n; i++){
    Rcpp::CharacterVector le = cfg[i];
    INFO(Rcpp::Rcout, get_sim_id(), le );
  }

  INFO(Rcpp::Rcout, get_sim_id(), "End of configuration." << std::endl );
}



void Trial::print_interims(){
  
  INFO(Rcpp::Rcout, get_sim_id(), "Interims for this trial: " );
  
  Rcpp::Rcout <<  intrms << std::endl;

  
}



//