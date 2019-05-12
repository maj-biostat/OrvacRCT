// #include <Rcpp.h>
// 
// #include "trial.h"
// using namespace Rcpp;


// namespace Rcpp {
// 
//   template <typename T> SEXP wrap_(const Trial& t) {
//     
//     
//     
//     std::string klass = "dgCMatrix";
//     S4 s(klass);
//    
//     return s;
//   }
// 
// }
// 



// // see http://dirk.eddelbuettel.com/code/rcpp/Rcpp-modules.pdf
// RCPP_MODULE(trial) {
//   Rcpp::class_<Trial>( "Trial" )
//   .constructor<int,int>()
//   .method( "is_v_samp_stopped", &Trial::is_v_samp_stopped )
//   .method( "is_immu_fut", &Trial::is_immu_fut )
//   .method( "is_clin_fut", &Trial::is_clin_fut )
//   .method( "is_clin_es", &Trial::is_clin_es )
//   .method( "is_inconclusive", &Trial::is_inconclusive )
//   .method( "get_immu_final", &Trial::get_immu_final )
//   .method( "get_clin_final", &Trial::get_clin_final )
//   ;
// }
