#include <RcppArmadillo.h>

#include <string>
#include <sstream>
#include "trial.h"


Trial::Trial(Rcpp::List& trial_cfg, int idxsim){
  
  sim_id = idxsim;
  cfg = trial_cfg;
  d = init_trial_dat();
  intrms = get_interims();
  
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
  
  d.col(COL_CEN) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_OBST) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_REASON) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_IMPUTE) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  d.col(COL_REFTIME) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["n_stop"]));
  
  return d;
}

arma::mat Trial::get_interims(){
  
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
    
    // every three months or next 50, whichever sooner
    // start the clock from the accrual of the previous kid
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



void Trial::run_interims(){

  int i = 0;
  for(i = 0; i < intrms.n_rows ;  i++){
    
    // keep track of how many currently enrolled
    this->idx_cur_intrm = i;
    this->n_enrolled = intrms(i, INT_N);
    this->immu_ss = intrms(i, INT_N);
    this->clin_ss = intrms(i, INT_N);
    
    INFO(Rcpp::Rcout, get_sim_id(), "intrm  " << idx_cur_intrm
      << " n_enrolled " << n_enrolled
      << " do immu " << do_immu()
      << " do clin " << do_clin());
    
    
    if(do_clin()){
      
      // set idx_cur_intrm prior to use
      clin_interim();
      
      
      INFO(Rcpp::Rcout, get_sim_id(), "clin c_ppos_n " << c_ppos_n 
        << " clin c_ppos_max " << c_ppos_max);
      
      if(c_ppos_max < (double)cfg["thresh_pp_fut"]){
        INFO(Rcpp::Rcout, get_sim_id(), "clin futile - stopping now, ppmax " 
          << c_ppos_max
          << " fut thresh " << (double)cfg["thresh_pp_fut"]);
        set_clin_fut();
        break;
      }
      
      if (c_ppos_n > (double)cfg["thresh_pp_es"] && !is_clin_fut()){
        INFO(Rcpp::Rcout, get_sim_id(), "clin es - stopping now, ppn " 
          << c_ppos_n
          << " es thresh " << (double)cfg["thresh_pp_es"]);
        set_clin_es();
        break;
      }
    }
  }

}

void Trial::clin_interim(){
  

  INFO(Rcpp::Rcout, get_sim_id(), "doing clin intrm, to n = "
        << (int)intrms(idx_cur_intrm, INT_N) 
        << ", enrld " << get_enrolled_ss()
        << ", thresh_pp_fut " << get_thresh_pp_fut()
        << ", thresh_pp_es " << get_thresh_pp_es()
        << ", thresh_p_sup " << get_thresh_p_sup() 
        << " post_draw " << (int)cfg["post_draw"]);
  
      
  // stores samples from posterior, 
  arma::mat m = arma::zeros((int)cfg["post_draw"] , 3);
  arma::mat m_pp_int = arma::zeros((int)cfg["post_draw"] , 3);
  arma::mat m_pp_max = arma::zeros((int)cfg["post_draw"] , 3);

  arma::uvec ugt1;
  int int_win = 0;
  int max_win = 0;
  
  double a = (double)cfg["prior_gamma_a"];
  double b = (double)cfg["prior_gamma_a"];
  
  c_suf_intrm = clin_set_state((int)intrms(idx_cur_intrm, INT_N),
    (double)intrms(idx_cur_intrm, INT_T_END), 0);
  
  // Rcpp::Rcout << " interim ss " << (int)intrms(idx_cur_intrm, INT_N)
  //             << " fu " << 0
  //             << " n_event_0 " << (int)c_suf_intrm["n_evnt_0"]
  //             << " n_event_1 " << (int)c_suf_intrm["n_evnt_1"]
  //             << " tot_obst_0 " << (double)c_suf_intrm["tot_obst_0"]
  //             << " tot_obst_1 " << (double)c_suf_intrm["tot_obst_1"]
  //             << std::endl;
  
  // keep a copy of the original state
  arma::mat d_orig = arma::zeros((int)cfg["n_stop"], 6);
  d_orig.col(0) = arma::vec(d.col(COL_EVTT));
  d_orig.col(1) = arma::vec(d.col(COL_CEN));
  d_orig.col(2) = arma::vec(d.col(COL_OBST));
  d_orig.col(3) = arma::vec(d.col(COL_REASON));
  d_orig.col(4) = arma::vec(d.col(COL_IMPUTE));
  d_orig.col(5) = arma::vec(d.col(COL_REFTIME));
  
  // subjs that require imputation
  arma::uvec uimpute = arma::find(d.col(COL_IMPUTE) == 1);

 
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
    m(i, COL_RATIO) = m(i, COL_LAMB0) / m(i, COL_LAMB1);
    

    // Rcpp::Rcout << " m(i, COL_RATIO)      " << m(i, COL_RATIO) << std::endl;
    //Rcpp::Rcout << " (int)uimpute.n_elem  " << (int)uimpute.n_elem << std::endl;
    
    // use memoryless prop of exponential and impute enrolled kids that have not
    // yet had event.
    for(int j = 0; j < (int)uimpute.n_elem; j++){
      int sub_idx = uimpute(j);
      if(d(sub_idx, COL_TRT) == 0){
        
        // Rcpp::Rcout << "ctl d(sub_idx, COL_EVTT) was " 
        //         << d(sub_idx, COL_EVTT)
        //         << std::endl;
        
        d(sub_idx, COL_EVTT) = d(sub_idx, COL_OBST) + R::rexp(1/m(i, COL_LAMB0))  ;
        
        // Rcpp::Rcout << "ctl d(sub_idx, COL_EVTT) is " 
        //         << d(sub_idx, COL_EVTT)
        //         << std::endl;
      } else {
        // Rcpp::Rcout << "trt d(sub_idx, COL_EVTT) was " 
        //         << d(sub_idx, COL_EVTT)
        //         << std::endl;
        d(sub_idx, COL_EVTT) = d(sub_idx, COL_OBST) + R::rexp(1/m(i, COL_LAMB1))  ;
        
        // Rcpp::Rcout << "trt d(sub_idx, COL_EVTT) is " 
        //         << d(sub_idx, COL_EVTT)
        //         << std::endl;
      }
    }

    
    
    // update view of the sufficent stats using enrolled
    // kids that have all now been given an event time
    c_suf_intrm_fu = clin_set_state((int)intrms(idx_cur_intrm, INT_N), 
      (double)intrms(idx_cur_intrm, INT_T_END),
      (double)cfg["max_age_fu"]);
    
    // Rcpp::Rcout << "c_suf_intrm_fu " 
    //             << " n_event_0 " << (int)c_suf_intrm_fu["n_evnt_0"]
    //             << " n_event_1 " << (int)c_suf_intrm_fu["n_evnt_1"]
    //             << " tot_obst_0 " << (double)c_suf_intrm_fu["tot_obst_0"]
    //             << " tot_obst_1 " << (double)c_suf_intrm_fu["tot_obst_1"]
    //             << std::endl;

    for(int j = 0; j < (int)cfg["post_draw"]; j++){
      m_pp_int(j, COL_LAMB0) = R::rgamma(a + (double)c_suf_intrm_fu["n_evnt_0"], 
        1/(b + (double)c_suf_intrm_fu["tot_obst_0"]));
      m_pp_int(j, COL_LAMB1) = R::rgamma(a + (double)c_suf_intrm_fu["n_evnt_1"], 
        1/(b + (double)c_suf_intrm_fu["tot_obst_1"]));
      m_pp_int(j, COL_RATIO) = m_pp_int(j, COL_LAMB0) / m_pp_int(j, COL_LAMB1);
    }
    // empirical posterior probability that ratio_lamb > 1
    ugt1 = arma::find(m_pp_int.col(COL_RATIO) > 1);
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
    // inside clin_set_state
    c_suf_max_fu = clin_set_state((int)d.n_rows, 0, (double)cfg["max_age_fu"]);

    // what does the posterior at max sample size say?
    for(int j = 0; j < (int)cfg["post_draw"]; j++){
      m_pp_max(j, COL_LAMB0) = R::rgamma(a + (double)c_suf_max_fu["n_evnt_0"], 
        1/(b + (double)c_suf_max_fu["tot_obst_0"]));
      m_pp_max(j, COL_LAMB1) = R::rgamma(a + (double)c_suf_max_fu["n_evnt_1"], 
        1/(b + (double)c_suf_max_fu["tot_obst_1"]));
      m_pp_max(j, COL_RATIO) = m_pp_max(j, COL_LAMB0) / m_pp_max(j, COL_LAMB1);
    }
    // empirical posterior probability that ratio_lamb > 1
    ugt1 = arma::find(m_pp_max.col(COL_RATIO) > 1);
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
  
  INFO(Rcpp::Rcout, get_sim_id(), "clin int_win " << int_win 
    << " max_win " << max_win);

  this->c_ppos_n = (double)int_win/(double)cfg["post_draw"];
  this->c_ppos_max = (double)max_win/(double)cfg["post_draw"];

  return;
}

// set censoring status up to participant n_target (1 to 1000)
// with or without fu.
Rcpp::List Trial::clin_set_state(int n_target, double ref_time, double fu){

  int n_evnt_0 = 0;
  int n_evnt_1 = 0;
  double tot_obst_0 = 0;
  double tot_obst_1 = 0;

  // set censoring and event times up to current enrolled
  // these kids were all enrolled prior to the current look
  for(int sub_idx = 0; sub_idx < n_target; sub_idx++){

    // if we are looking at the max sample size it means
    // we are either imputing to the max or doing a final 
    // analysis. either way the correct ref_time is individual 
    // specific and unknowable until we look at the age and 
    // accrual time for that individual
    if(n_target == d.n_rows){
      //Rcpp::Rcout <<  " updating ref time " << std::endl;
      ref_time = (double)cfg["max_age_fu"] - d(sub_idx, COL_AGE) +
        d(sub_idx, COL_ACCRT);
    }
    
    if(fu == 0){
      d(sub_idx, COL_REFTIME) = ref_time;
      
    } else {

      if(fu - d(sub_idx, COL_AGE) + d(sub_idx, COL_ACCRT) > ref_time){
        d(sub_idx, COL_REFTIME) = fu - d(sub_idx, COL_AGE) + d(sub_idx, COL_ACCRT);
        
      } else {
        
        d(sub_idx, COL_REFTIME) = ref_time;
        
      }

    }

    if(d(sub_idx, COL_ACCRT) + d(sub_idx, COL_EVTT) <= d(sub_idx, COL_REFTIME) &&
      d(sub_idx, COL_AGE) + d(sub_idx, COL_EVTT) <= (double)cfg["max_age_fu"]){
      // observed event
      // dont impute
      d(sub_idx, COL_CEN) = 0;
      d(sub_idx, COL_OBST) = d(sub_idx, COL_EVTT);
      d(sub_idx, COL_REASON) = 1;
      d(sub_idx, COL_IMPUTE) = 0;

      if(d(sub_idx, COL_TRT) == 0) {
        n_evnt_0 += 1;
      } else {
        n_evnt_1 += 1;
      }

    } else if (d(sub_idx, COL_ACCRT) + d(sub_idx, COL_EVTT) <= d(sub_idx, COL_REFTIME) &&
      d(sub_idx, COL_AGE) + d(sub_idx, COL_EVTT) > (double)cfg["max_age_fu"]){
      // censor at max age
      // dont impute
      d(sub_idx, COL_CEN) = 1;
      d(sub_idx, COL_OBST) = (double)cfg["max_age_fu"] - d(sub_idx, COL_AGE);
      d(sub_idx, COL_REASON) = 2;
      d(sub_idx, COL_IMPUTE) = 0;

    } else if (d(sub_idx, COL_ACCRT) + d(sub_idx, COL_EVTT) > d(sub_idx, COL_REFTIME) &&
      d(sub_idx, COL_REFTIME) - d(sub_idx, COL_ACCRT) <=
                  (double)cfg["max_age_fu"] - d(sub_idx, COL_AGE)){
      // censor at mnth - accrual (t2)
      // impute
      d(sub_idx, COL_CEN) = 1;
      d(sub_idx, COL_OBST) = d(sub_idx, COL_REFTIME) - d(sub_idx, COL_ACCRT) ;
      d(sub_idx, COL_REASON) = 3;
      d(sub_idx, COL_IMPUTE) = 1;

    } else { // mnth - accrual >= max age - age at accrual
      // censor at max age
      // dont impute
      d(sub_idx, COL_CEN) = 1;
      d(sub_idx, COL_OBST) = (double)cfg["max_age_fu"] - d(sub_idx, COL_AGE) ;
      d(sub_idx, COL_REASON) = 4;
      d(sub_idx, COL_IMPUTE) = 0;

    }

    if(d(sub_idx, COL_TRT) == 0) {
      tot_obst_0 = tot_obst_0 + d(sub_idx, COL_OBST);
    } else {
      tot_obst_1 = tot_obst_1 + d(sub_idx, COL_OBST);
    }
  }
  
  // Rcpp::Rcout << d << std::endl;

  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("n_evnt_0") = n_evnt_0,
                           Rcpp::Named("tot_obst_0") = tot_obst_0,
                           Rcpp::Named("n_evnt_1") = n_evnt_1,
                           Rcpp::Named("tot_obst_1") = tot_obst_1,
                           Rcpp::Named("fu") = fu);
  
  return ret;
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
void Trial::set_immu_final(bool won){i_final_win = won;}
void Trial::set_clin_final(bool won){c_final_win = won;}

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

int Trial::get_sim_id(){return sim_id;}
int Trial::get_nmaxsero(){return nmaxsero;}
int Trial::get_nstartclin(){return nstartclin;}
int Trial::get_immu_ss(){return immu_ss;}
int Trial::get_clin_ss(){return clin_ss;}
int Trial::get_enrolled_ss(){return n_enrolled;}
int Trial::get_immu_final(){return i_final_win;}
int Trial::get_clin_final(){return c_final_win;}

double Trial::get_i_ppos_n(){return i_ppos_n;}
double Trial::get_i_ppos_max(){return i_ppos_max;} 
double Trial::get_c_ppos_n(){return c_ppos_n;}
double Trial::get_c_ppos_max(){return c_ppos_max;} 

double Trial::get_thresh_pp_fut(){return (double)cfg["thresh_pp_fut"];} 
double Trial::get_thresh_pp_es(){return (double)cfg["thresh_pp_es"];} 
double Trial::get_thresh_p_sup(){return (double)cfg["thresh_p_sup"];} 

int Trial::get_num_interims(){return intrms.n_rows;}
int Trial::get_interim_n(int i){return intrms(i, INT_N);}
int Trial::get_interim_time(int i){return intrms(i, INT_T_END);}

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
  ret["immu_final"] = get_immu_final();
  ret["clin_final"] = get_clin_final();
  ret["v_samp_stopped"] = is_v_samp_stopped();
  ret["immu_fut"] = is_immu_fut();
  ret["clin_fut"] = is_clin_fut();
  ret["clin_es"] = is_clin_es();
  ret["inconclusive"] = is_inconclusive();
  
  return ret;
}

void Trial::print_state(){
  Rcpp::Rcout << "INFO: " << __FILE__ << "(" << __LINE__ << ") "
    << " sim = " << get_sim_id() << " "
    << "immu ep state: intrm stop v samp " << is_v_samp_stopped()
    << " futile " << is_immu_fut()
    << " inconclusive " << is_inconclusive()
    << " fin analy win " << get_clin_final()
    << std::endl;
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