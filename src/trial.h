#ifndef TRIAL_H
#define TRIAL_H

#include <RcppArmadillo.h>


// column indices
#define COL_ID            0
#define COL_TRT           1
#define COL_ACCRT         2
#define COL_AGE           3
#define COL_SEROT2        4
#define COL_SEROT3        5
#define COL_PROBT3        6
#define COL_IMPUTESERO    7
#define COL_EVTT          8
#define COL_CEN           9
#define COL_OBST          10
#define COL_REASON        11
#define COL_IMPUTEEVTT    12
#define COL_REFTIME       13
#define NCOL              14

#define COL_THETA0        0
#define COL_THETA1        1
#define COL_DELTA         2

#define COL_LAMB0         0
#define COL_LAMB1         1
#define COL_RATIO         2

#define INT_I_START       0
#define INT_I_END         1
#define INT_T_START       2
#define INT_T_END         3
#define INT_N             4

class Trial {
private:

  int stop_ven_samp = 0;
  int stop_immu_fut = 0;
  int stop_clin_fut = 0;
  int stop_clin_es = 0;
  int inconclu = 0;
  
  int nmaxsero = 250;
  int nstartclin = 200;
  
  int i_final_decision = 0;
  int c_final_decision = 0;

protected:
  
  int sim_id = 0;
  Rcpp::List cfg;
  arma::mat d;
  arma::mat intrms;
  
  // dynamic record sample size 
  int immu_ss = 0;
  int clin_ss = 0;
  int n_enrolled = 0;
  int idx_cur_intrm = 0;
  bool is_final = 0;
  
  // sufficient stats for clinical
  Rcpp::List c_suf_fin;
  Rcpp::List c_suf_intrm;
  Rcpp::List c_suf_intrm_fu;
  Rcpp::List c_suf_max_fu;
  
  // sufficient stats for immu
  Rcpp::List i_suf_fin;
  Rcpp::List i_suf_intrm;
  Rcpp::List i_suf_intrm_fu;
  Rcpp::List i_suf_max_fu;
  
  // current posterior probability of trt effect
  double i_p_n = 0;
  double c_p_n = 0;
  
  // predictive prob at time of current interim
  double i_ppos_n = 0;
  double i_ppos_max = 0; 
  double c_ppos_n = 0;
  double c_ppos_max = 0;
  
  // posterior prob at final analysis
  double i_p_fin = 0; 
  double c_p_fin = 0;
  

public:
  
  Trial(Rcpp::List& trial_cfg, int idxsim);
  
  arma::mat init_trial_dat();
  arma::mat init_interims();
  arma::mat get_interims();
  
  void clin_interim();
  Rcpp::List clin_censoring();
  Rcpp::List clin_censoring(int n_target, double ref_time);
  void clin_censoring_subj(int sub_idx);
  void clin_fin();
  
  void immu_interim();
  Rcpp::List immu_observed();
  Rcpp::List immu_observed(int n_target);
  void immu_fin();
  
  void run_interims();
  void run_final();
  
  bool do_immu();
  bool do_clin();
  
  void set_nmaxsero(int nmaxs);
  void set_nstartclin(int nstartc);
  void set_immu_ss(int n);
  void set_clin_ss(int n);
  void set_enrolled_ss(int n);
  void set_immu_fin_decision(int decision);
  void set_clin_fin_decision(int decision); 
  
  void set_curr_intrm_idx(int cur_intrm);

  void set_i_ppos_n(double ppos_n);
  void set_i_ppos_max(double ppos_max);
  void set_c_ppos_n(double ppos_n);
  void set_c_ppos_max(double ppos_max);  
  
  void set_immu_stopv();
  void set_immu_fut();
  void set_clin_fut();
  void set_clin_es();
  void set_inconclusive();

  
  arma::mat get_data();
  int get_sim_id();
  int get_nmaxsero();
  int get_nstartclin();
  int get_immu_ss();
  int get_clin_ss();
  int get_enrolled_ss();
  int get_immu_fin_decision();
  int get_clin_fin_decision();

  // getters for predictive probs
  double get_i_ppos_n();
  double get_i_ppos_max(); 
  double get_c_ppos_n();
  double get_c_ppos_max(); 
  // posterior prob at current interim
  double get_i_p_n();
  double get_c_p_n();
  // final posterior prob
  double get_i_p_fin();
  double get_c_p_fin();
  
  double get_thresh_pp_fut();
  double get_thresh_pp_es();
  double get_thresh_p_sup();

  int get_num_interims();
  int get_interim_n(int i);
  int get_interim_time(int i);
  int get_curr_intrm_idx();
  
  int is_v_samp_stopped();
  int is_immu_fut();
  int is_clin_fut();
  int is_clin_es();
  int is_inconclusive();
  
  Rcpp::List as_list();
  
  void print_state();
  void print_cfg();
  void print_interims();

};

#endif