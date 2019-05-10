#ifndef TRIAL_H
#define TRIAL_H

#include <Rcpp.h>


class Trial {
private:
  int stop_ven_samp = 0;
  int stop_immu_fut = 0;
  int stop_clin_fut = 0;
  int stop_clin_es = 0;
  int inconclu = 0;
  int nmaxsero = 250;
  int nstartclin = 200;
  int immu_ss = 0;
  int clin_ss = 0;
  
  bool i_final_win = 0;
  bool c_final_win = 0;
  
public:
  Trial();
  Trial(int nmaxs, int nstartc);
  Trial(int vstop, int ifut, int cfut, int ces, int inc, 
        int nmaxs, int nstartc);
  bool do_immu(int n_current);
  bool do_clin(int n_current);
  
  void set_nmaxsero(int nmaxs);
  void set_nstartclin(int nstartc);
  void set_immu_ss(int n);
  void set_clin_ss(int n);
  void set_immu_final(bool won);
  void set_clin_final(bool won); 
  
  void set_immu_stopv();
  void set_immu_fut();
  void set_clin_fut();
  void set_clin_es();
  void set_inconclusive();
  
  int get_nmaxsero();
  int get_nstartclin();
  int get_immu_ss();
  int get_clin_ss();
  int get_immu_final();
  int get_clin_final();
  
  int is_v_samp_stopped();
  int is_immu_fut();
  int is_clin_fut();
  int is_clin_es();
  int is_inconclusive();
};


#endif