// misc.h
//#include "RcppArmadillo.h"

// Soft-threshold operator
arma::vec soft_thresh(const arma::vec & z, double l){
  arma::vec S = arma::vec(z.n_elem);
  for(int i=0; i<z.n_elem; i++){
    if(std::abs(z[i]) <= l){
      S(i) = 0;
    }else{
      if(z[i]<=0){
        S[i] = z[i] + l;
      }else{
        S[i] = z[i] - l;
      }
    }
  }
  return S;
}

void one_rank_svd(arma::vec & v, arma::vec & u,  arma::mat & Eta_n){
   // We will to compute the 1-rank svd... we are going to iterate using a "power" method
  v = arma::randu(Eta_n.n_cols);
  u = arma::randu(Eta_n.n_rows);
  u = u/arma::norm(u, 2);

  arma::vec v_old;
  for (int iter = 0; iter < 100; iter++)
  {
      v_old = v;
      v = Eta_n.t()*u;
      u = Eta_n*v;
      if (arma::norm(u, 2) > 0)
      {
        u = u/arma::norm(u, 2);
      } 
      if(arma::approx_equal(v, v_old, "reldiff", 0.1)){
        break;
      }
  }
}

