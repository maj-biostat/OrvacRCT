#include <string>
#include <sstream>
#include <Rcpp.h>

#include "trial.h"



Trial::Trial()
{
  nmaxsero = 250;
  nstartclin = 200;
}
Trial::Trial(int nmaxs, int nstartc)
{
  nmaxsero = nmaxs;
  nstartclin = nstartc;
}
Trial::Trial(int vstop, int ifut, int cfut, int ces, int inc, 
             int nmaxs, int nstartc)
{
  nmaxsero = nmaxs;
  nstartclin = nstartc;
  stop_ven_samp = vstop;
  stop_immu_fut = ifut;
  stop_clin_fut = cfut;
  stop_clin_es = ces;
  inconclu = inc;
}

bool Trial::do_immu(int n_current){
  if(stop_ven_samp == 1){
    return false;
  }
  if(n_current > nmaxsero){
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
bool Trial::do_clin(int n_current){
  
  if(n_current < nstartclin){
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
void Trial::set_immu_final(bool won){i_final_win = won;}
void Trial::set_clin_final(bool won){c_final_win = won;}

void Trial::set_immu_stopv(){stop_ven_samp = 1;}
void Trial::set_immu_fut(){stop_immu_fut = 1;}
void Trial::set_clin_fut(){stop_clin_fut = 1;}
void Trial::set_clin_es(){stop_clin_es = 1;}
void Trial::set_inconclusive(){inconclu = 1;}

int Trial::get_nmaxsero(){return nmaxsero;}
int Trial::get_nstartclin(){return nstartclin;}
int Trial::get_immu_ss(){return immu_ss;}
int Trial::get_clin_ss(){return clin_ss;}
int Trial::get_immu_final(){return i_final_win;}
int Trial::get_clin_final(){return c_final_win;}

int Trial::is_v_samp_stopped(){return stop_ven_samp;}
int Trial::is_immu_fut(){return stop_immu_fut;}
int Trial::is_clin_fut(){return stop_clin_fut;}
int Trial::is_clin_es(){return stop_clin_es;}
int Trial::is_inconclusive(){return inconclu;}

// Basic test cases/sanity checking.

// [[Rcpp::export]]
void test_trial1() {
  
  Trial t;
  bool test_result = false;
  
  assert(t.do_immu(10) == 1);
  assert(t.do_immu(1000) == 1);
  
  assert(t.do_clin(10) == 1);
  assert(t.do_clin(300) == 1);
  
  assert(t.get_nmaxsero() == 250);
  t.set_nmaxsero(10);
  assert(t.get_nmaxsero() == 10);
  
  assert(t.get_nstartclin() == 200);
  t.set_nstartclin(11);
  assert(t.get_nstartclin() == 11);

  assert(t.get_immu_ss() == 0);
  t.set_immu_ss(13);
  assert(t.get_immu_ss() == 13);
  
  assert(t.get_clin_ss() == 0);
  t.set_clin_ss(14);
  assert(t.get_clin_ss() == 14);
  
  assert(t.get_immu_final() == 0);
  t.set_immu_final(1);
  assert(t.get_immu_final() == 1);
  t.set_immu_final(0);
  assert(t.get_immu_final() == 0);
  
  assert(t.get_clin_final() == 0);
  t.set_clin_final(1);
  assert(t.get_clin_final() == 1);
  t.set_clin_final(0);
  assert(t.get_clin_final() == 0);

  Rcpp::Rcout << "all tests passed"  << std::endl;
  // 
  // t.set_clin_ss(10);
  // 
  // Rcpp::Rcout << "test 2 " << t.get_clin_ss() << std::endl;
  
  return ;
}
// [[Rcpp::export]]
void test_trial2() {
  
  Trial t;
  bool test_result = false;
  
  assert(t.is_v_samp_stopped() == 0);
  t.set_immu_stopv();
  assert(t.is_v_samp_stopped() == 1);
  
  assert(t.is_immu_fut() == 0);
  t.set_immu_fut();
  assert(t.is_immu_fut() == 1);
  
  assert(t.is_clin_fut() == 0);
  t.set_clin_fut();
  assert(t.is_clin_fut() == 1);
  
  assert(t.is_clin_es() == 0);
  t.set_clin_es();
  assert(t.is_clin_es() == 1);  
  
  assert(t.is_inconclusive() == 0);
  t.set_inconclusive();
  assert(t.is_inconclusive() == 1);  
  
  Rcpp::Rcout << "all tests passed"  <<  std::endl;
  
  // 
  // Rcpp::Rcout << "test 1 " << t.get_clin_ss() << std::endl;
  // 
  // t.set_clin_ss(10);
  // 
  // Rcpp::Rcout << "test 2 " << t.get_clin_ss() << std::endl;
  
  return ;
}
// [[Rcpp::export]]
void test_trial3() {
  
  Trial t1;
  assert(t1.do_immu(10) == 1);
  t1.set_immu_fut();
  assert(t1.do_immu(10) == 0);
  
  Trial t2;
  assert(t2.do_immu(10) == 1);
  t2.set_clin_fut();
  assert(t2.do_immu(10) == 0);
  
  Trial t3;
  assert(t3.do_immu(10) == 1);
  t3.set_clin_es();
  assert(t3.do_immu(10) == 0);
  
  Rcpp::Rcout << "all tests passed"  <<  std::endl;
  
  return ;
}
// [[Rcpp::export]]
void test_trial4() {
  
  Trial t1;
  assert(t1.do_clin(300) == 1);
  t1.set_immu_fut();
  assert(t1.do_immu(300) == 0);
  
  Trial t2;
  assert(t2.do_clin(300) == 1);
  t2.set_clin_es();
  assert(t2.do_immu(300) == 0);
  
  Trial t3;
  assert(t3.do_immu(300) == 1);
  t3.set_immu_fut();
  assert(t3.do_immu(300) == 0);
  
  Rcpp::Rcout << "all tests passed"  <<  std::endl;
  
  return ;
}






//