// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


#include "losses.h"
#include "misc.h"

#define TYPE_PREDICT_DEFAULT 0
#define TYPE_PREDICT_RESPONSE 1
#define TYPE_PREDICT_PROB 2
#define TYPE_PREDICT_CLASS 3

#define FISTA_MAX_ITER_INNER 500
#define FISTA_TOL 1e-4
#define FISTA_T0 2
#define FISTA_STEP 0.5

#define GLASP_MAX_ITER 200
#define GLASP_SGL_MAX_ITER 100

#define MAX_ITER_OPT 500

// glasp is the main class.
class glasp{
private:
    // Data matrix
    arma::mat X;
    // Labels
    arma::mat y;
    // Latent groups
    arma::mat T;
    // Cluster information matrix
    arma::mat W;
    // variable order
    arma::uvec ord;
    // groups
    arma::vec grp_start;
    arma::vec grp_end;
    arma::vec grp_len;

    // Eta
    arma::mat Eta;

    //scale params
    arma::rowvec X_col_means;
    arma::rowvec X_col_stds;
    arma::uvec rm_cols;

    // hyperparameters
    double lambda1;
    double lambda2;
    double lambda3;
    // Number of true groups
    int K;
    // Number of desired groups (max.)
    int K_max;
    // number of features
    int p;
    // number of observations
    int N;
    // Loss function
    Loss R;
    // coefficients
    arma::vec beta;
    double intercept;

    double mean_y;

    // type of loss
    int loss;

    // FISTA config
    bool use_warmstart;
    int Arg_FISTA_MAX_ITER_INNER;
    double Arg_FISTA_TOL;
    double Arg_FISTA_T0;
    double Arg_FISTA_STEP;
    bool use_default;

    // Null model's objective
    double obj_null;

    void update_svds(const arma::vec & eta_times_u, const arma::vec & beta, arma::vec & v, double gamma, double C_b, double C_v);
    void update_grp_order();
    void reorder_grp();
    void optimization_internal(int group_k);
    void optimization_external();
    void update_Eta_k(arma::vec & beta_k, int k);
    void update_Eta();

    void set_beta(arma::vec & beta0){
      beta = beta0;
      update_Eta();
    }
    void set_beta_k(arma::vec & beta0, int k){
      beta.subvec(grp_start(k), grp_end(k)) = beta0;
      update_Eta_k(beta0, k);
    }

    double L(arma::vec & beta_k, arma::mat & X_k, arma::mat & eta_minus_k, arma::mat & y){
      double glasp_penalty = 0;
      int p_k = beta_k.n_elem;
      for (int j = 0; j < p_k; j++)
      {
        glasp_penalty +=  arma::sum(arma::square(X_k.col(j)*beta_k(j) - T*W.row(j).t()));
      }
      return R.R2(beta_k, X_k, eta_minus_k, y) + lambda3/2*glasp_penalty;
    }

    arma::vec gradL(arma::vec & beta_k, arma::mat & X_k, arma::mat & eta_minus_k, arma::mat & y){
      arma::vec gradR = R.gradR2(beta_k, X_k, eta_minus_k, y);
      int p_k = beta_k.n_elem;
      arma::vec gradP(p_k);
      for (int j = 0; j < p_k; j++)
      {
        gradP(j) = arma::as_scalar(lambda3*X_k.col(j).t()*(X_k.col(j)*beta_k(j) - T*W.row(j).t()));
      }
      return gradR + gradP;
    }

    arma::vec Update(const arma::vec & beta0, const arma::vec & gradL_beta0, double t){
      int p_k = beta0.n_elem;
      arma::vec S1 = soft_thresh(beta0 - t*gradL_beta0, t*lambda1);
      return std::max(1 - lambda2*std::sqrt(p_k)*t/arma::norm(S1, 2), 0.0)*S1;
    }

public:
// Constructor
    glasp(const Rcpp::List gData, int loss);
// glasp objective function
    double objective();
// Main optimization
    void glasp_optimization(const arma::vec & params);

// Compute the groups
    void glasp_subproblem2();

// predictions
    arma::vec predict_response(const arma::mat & newX);
    arma::vec predict_class(const arma::mat & newX);
    arma::vec predict_probability(const arma::mat & newX);
    arma::vec predict(const arma::mat & newX, int type);

    // get/set
    arma::vec get_beta(){
      //re-scale beta
      arma::vec coefs = beta/X_col_stds.t();

      // add-removed variables
      if(rm_cols.n_elem > 0){
        for (unsigned int i = 0; i < rm_cols.n_elem; i++)
        {
            coefs.insert_rows(rm_cols(i), 1);
        }
      }
      return coefs;
    }
    double get_intercept(){
      //re-scale intercept
      double c0 = intercept - arma::as_scalar((X_col_means/X_col_stds) * beta);
      return c0;}
    arma::mat get_X(){return X;}
    arma::mat get_y(){return y;}
    arma::mat get_T(){return T;}
    arma::mat get_W(){return W;}
    arma::vec get_clusters(){
      int K = W.n_cols;
      // Find the groups
      arma::vec grp_index = arma::zeros(p);
      for (int j = 0; j < p; j++)
      {
        for (int k = 0; k < K; k++)
        {
          if( W(j,k)!=0 ){
            grp_index(j) = k;
            break;
          }
        }
      }
      // add-removed variables
      if(rm_cols.n_elem > 0){
        for (unsigned int i = 0; i < rm_cols.n_elem; i++)
        {
            grp_index.insert_rows(rm_cols(i), arma::ones(1)*K);
        }
      }
      return grp_index;
    }
    arma::uvec get_ord(){return ord;}
};

/* Constructor */
glasp::glasp(const Rcpp::List glaspData, int loss)
{
    this->X = Rcpp::as<arma::mat>(glaspData["X"]);
    this->y = Rcpp::as<arma::mat>(glaspData["y"]);
    this->loss = loss;

    // X must be scaled!
    // 1 - remove constant columns
    X_col_stds = arma::stddev(X, 0);
    rm_cols = arma::find(X_col_stds < 1e-16);
    X.shed_cols(rm_cols);

    // 2- scale
    this->X_col_means = arma::mean(X, 0);
    this->X_col_stds = arma::stddev(X, 0);
    this->X = X.each_row() - X_col_means;
    this->X = X.each_row()/X_col_stds;


    p = X.n_cols;
    N = X.n_rows;
    this->R = Loss(loss);

    // Default values
    lambda1 = 0;
    lambda2 = 0;
    lambda3 = 0;
    K = 1;
    K_max = 3;
    ord = arma::regspace<arma::uvec>(0, p-1);

    use_warmstart = FALSE;
    use_default = TRUE;

    W = arma::mat(p, 1, arma::fill::ones);
    T = arma::mat(N, 1, arma::fill::zeros);

    beta = arma::zeros(p);
    Eta = arma::zeros(N, p);

    update_grp_order();

    double obj_null = objective();
    // Compute intercept... or not
    arma::uvec y_ord;
    switch (loss)
    {
    case LOSS_LINEAR:
        intercept = arma::mean(arma::mean(y));
        y = y - intercept;
        mean_y = 0;
        break;
    case LOSS_LOGIT:
        intercept = 0;
        mean_y = arma::mean(arma::mean(y));
        break;
    case LOSS_COX:
        intercept = 0;
        // We have to sort y, and X, accordingly
        y_ord = arma::sort_index(y.col(0), "descend");
        this->y = y.rows(y_ord);
        this->X = X.rows(y_ord);
        break;
    default:
        break;
    }
}


// Optimization

void glasp::glasp_optimization(const arma::vec & params){
  // Update params
    lambda1 = params[0];
    lambda2 = params[1];
    lambda3 = params[2];
    K_max = params[3];

    // lambda1, lambda2, lambda3 are transformed

    // double R_null = R.R(arma::zeros(p), X, y);
    // lambda1 = lambda1/R_null;
    // lambda2 = lambda2/R_null;
    // lambda3 = lambda2*lambda3;

    // is warm start used?
    if(!use_warmstart){
        // initial beta
        beta = arma::zeros(p);
        update_Eta();
        W = arma::mat(p, 1, arma::fill::ones);
        T = arma::mat(N, 1, arma::fill::zeros);
    }
    obj_null = objective();

    double obj_new;
    double obj_old = objective();
    //printf("initial objective %f\n", obj_old);
    for (int iter = 0; iter < GLASP_MAX_ITER; iter++){
      // Compute groups
      if (lambda3 > 0)
      {
        glasp_subproblem2();
      }
      update_grp_order();

      // Compute beta
      optimization_external();

      obj_new = objective();
      //printf("objective %f\n", obj_new);
      //beta.print("beta");

      if(obj_old - obj_new < FISTA_TOL*obj_null){
        //printf("Exiting main loop because %f < %f\n", obj_old - obj_new, FISTA_TOL*obj_null);
        break;
      }
      obj_old = obj_new;
    }
    //ord.print("ord");
    reorder_grp();
}

void glasp::optimization_external(){
  // T,W,K are fixed
  K = W.n_cols;
  double obj_new;
  double obj_old = objective();
  for (int iter = 0; iter < GLASP_SGL_MAX_ITER; iter++)
  {
    // Iterate through the groups
    for (int k = 0; k < K; k++)
    {
      optimization_internal(k);
    }
    obj_new = objective();
    if(obj_old - obj_new < FISTA_TOL*obj_null){break;}
    obj_old = obj_new;
  }
}

void glasp::optimization_internal(int group_k){
  arma::mat X_k = X.cols(grp_start(group_k), grp_end(group_k));
  double p_k = X_k.n_cols;

  arma::vec theta_old = arma::zeros(p_k);
  arma::vec beta_0 = beta.subvec(grp_start(group_k), grp_end(group_k));
  arma::vec theta_new = beta_0;
  arma::mat eta_minus_k = X*beta - X_k*beta_0;
  // Compute the gradient at 0
  arma::vec g = gradL(theta_old, X_k, eta_minus_k, y);

  // Condition 1
  if(arma::norm(soft_thresh(g, lambda1), 2) <= lambda2*sqrt(p_k)){
    // Set beta 0
    set_beta_k(theta_old, group_k);
  }else{
    double t = FISTA_T0;
    double l_new = 1;
    double l_old = 1;
    double L_beta;
    double L_beta_new = L(beta_0, X_k, eta_minus_k, y);
    double obj_new;
    double obj_old = obj_null;

    for (int iter = 0; iter < FISTA_MAX_ITER_INNER; iter++){
        // This is U_{t_k-1}(\bbeta_{k-1})
        theta_old = theta_new;
        l_old = l_new;

        // Compute the gradient in beta_k
        g = gradL(beta_0, X_k, eta_minus_k, y);

        // Compute the risk in beta_k
        L_beta = L_beta_new;

        //Find t such that R_updated <= R + t(g)*(beta_updated - beta) + 1/2t||beta_updated - beta||_2^2
        theta_new = Update(beta_0, g, t);

        while (  L(theta_new, X_k, eta_minus_k, y) > L_beta + arma::as_scalar(g.t()*(theta_new - beta_0)) + 1/(2*t)*arma::sum(arma::square(theta_new - beta_0))){
            t = FISTA_STEP*t;
            theta_new = Update(beta_0, g, t);
        }

        // Compute the aceleration term
        l_new = (1 + sqrt(1 + 4*pow(l_old, 2)))/2;

        // Update beta
        beta_0 = theta_new + (l_old - 1)*(theta_new - theta_old)/l_new;
        set_beta_k(beta_0, group_k);

        // Compute the new risk
        L_beta_new = L(beta_0, X_k, eta_minus_k, y);

        // COmpute the new objective
        obj_new = objective();
        if(obj_old - obj_new < FISTA_TOL*obj_null){break;}
        obj_old = obj_new;
    }
  }
}

double glasp::objective(){
  double lasso_penalty = arma::norm(beta, 1);
  double grp_lasso_penalty = 0;
  for (unsigned int j = 0; j < grp_end.n_elem; j++)
  {
    grp_lasso_penalty += arma::norm(beta.subvec(grp_start(j), grp_end(j)), 2);
  }


  double glasp_penalty = std::pow(arma::norm(Eta - T*W.t(), "fro"), 2);

  return R.R(beta, X, y) +
         lambda1*lasso_penalty +
         lambda2*grp_lasso_penalty +
         lambda3/2*glasp_penalty;
}

void glasp::update_grp_order(){
   // K is fixed from W
  K = W.n_cols;

  // Find the groups
  arma::vec grp_index = arma::zeros(p);
  grp_len = arma::zeros(K);


  for (int j = 0; j < p; j++)
  {
    for (int k = 0; k < K; k++)
    {
      if( W(j,k)!=0 ){
        grp_index(j) = k;
        grp_len(k) = grp_len(k) + 1;
        break;
      }
    }
  }
  // Order X, W according to groups
  arma::uvec new_order = arma::sort_index(grp_index);
  ord = ord(new_order);
  X = X.cols(new_order);
  beta = beta(new_order);
  W = W.rows(new_order);
  Eta.cols(new_order);

  // Create the group structures
  grp_start = grp_len;
  grp_end = arma::cumsum(grp_len) - 1;
  grp_start(0) = 0;

  if (K > 1){
    grp_start.subvec(1, K-1) = arma::cumsum(grp_len.subvec(0,K-2));
  }

}

void glasp::reorder_grp(){
  arma::uvec reord = arma::sort_index(ord);

  X = X.cols(reord);
  beta = beta(reord);
  W = W.rows(reord);
  ord = ord(reord);
}

void glasp::glasp_subproblem2(){

  double gamma = 2*lambda2/lambda3;

  // 1- Remove null columns in Eta, creating Eta_n with the columns removed
  // beta_eq_zero stores the indices of zero elements of beta
  arma::vec beta_eq_zero;
  // beta_n stores the elements of beta different from zero.
  arma::vec beta_n;
  arma::mat Eta_n;
  int p_Eta_n = 0;
  int p_beta_eq_zero = 0;

  // Find sizes
  for (int i = 0; i < p; i++)
  {
    if( std::abs(beta(i)) < 1e-16 ){
      p_beta_eq_zero++;
    }else{
      p_Eta_n++;
    }
  }

  beta_eq_zero = arma::vec(p_beta_eq_zero);
  Eta_n = arma::mat(N, p_Eta_n);
  beta_n = arma::vec(p_Eta_n);
  p_Eta_n = 0;
  p_beta_eq_zero = 0;

  for (int i = 0; i < p; i++)
  {
    if( std::abs(beta(i)) < 1e-16 ){
      beta_eq_zero(p_beta_eq_zero++) = i;
    }else{
      Eta_n.col(p_Eta_n) = Eta.col(i);
      beta_n(p_Eta_n++) = beta(i);
    }
  }
  // There is an exception (is it?) is the number of desired groups exceeds the number of columns of Eta_n

  // 2 - Define the left and right singular vector matrices, as well as W, T

  arma::mat V, U;
  arma::vec v, u, v_old;
  arma::mat W, T, eta_times_u;

  // 3 - Compute the SparseSVD
  for (int k = 0; k < K_max; k++)
  {
    // Initial solution
    // We have to compute the 1-rank svd... we are going to iterate using a "power" method
    one_rank_svd(v, u, Eta_n);

    double C_b = 0.0;
    double C_v = 0.0;
    for (int j = 0; j < p_Eta_n; j++)
    {
      if(v(j)!= 0){
        C_b += pow(beta_n(j),2);
        C_v ++;
      }
    }
    // loop
    int  iter = 0;
    do
    {
      v_old = v;
      // update v
      eta_times_u = Eta_n.t()*u;
      update_svds(eta_times_u, beta_n, v, gamma, C_b, C_v);

      u = Eta_n*v;
      if (arma::norm(u, 2) > 0)
      {
        u = u/arma::norm(u, 2);
      }
    } while (arma::approx_equal(v_old, v, "absdiff", FISTA_TOL) && ++iter < MAX_ITER_OPT);

    double norm = arma::norm(v, 2);
    if(norm > 0){
      v = v/norm;
      u = u*norm;
    }

    W.insert_cols(k, v);
    T.insert_cols(k, u);
    Eta_n = Eta_n - u*v.t();
  }

  T.insert_cols(K_max, arma::zeros(N));


  // 4 - insert the zero rows in W
  for (unsigned int i = 0; i < beta_eq_zero.n_elem; i++)
  {
    W.insert_rows(beta_eq_zero(i), 1, true);
  }

  arma::mat W_sparse = arma::zeros<arma::mat>(p, K_max + 1);

  for (int j = 0; j < p; j++){
    double m = arma::abs(W.row(j)).max();
    if(m < 1e-8){
      W_sparse(j, K_max) = 1.0;
    }else{
      for (int k = 0; k < K_max; k++)
      {
        if ( std::abs(W(j, k)) >= m)
        {
          W_sparse(j, k) = W(j, k);
          break;
        }
      }
    }
  }
  // Remove null columns in W (just in case)

  int k = 0;
  K = K_max;
  while(k <= K){
    double d = arma::norm(W_sparse.col(k), "fro");
    if(d > FISTA_TOL){
      k++;
    }else{
      W_sparse.shed_col(k);
      T.shed_col(k);
      K--;
    }
  }
  this->W = W_sparse;
  this->T = T;

}

void glasp::update_svds(const arma::vec & eta_times_u, const arma::vec & beta, arma::vec & v, double gamma, double C_b, double C_v){
  // Maybe we don't have to loop at all...
  bool loop = FALSE;
  for (int l = 0; l < v.n_elem; l++)
  {
    if(pow(eta_times_u(l), 2) > gamma*std::abs(beta(l))){
      loop = TRUE;
      break;
    }
  }
  if(!loop){
    v = arma::zeros(v.n_elem);
    C_b = 0;
    C_v = 0;
  }

  // if we have to loop...
  int num_iter = 0;
  while(loop){
    loop = FALSE;
    for (int l = 0; l < v.n_elem; l++)
    {
      // Compute the updated C_v, C_b
      if(v(l) == 0){
        if ( pow(eta_times_u(l), 2) > gamma*( sqrt(C_b + pow(beta(l), 2))*sqrt(C_v + 1) - sqrt(C_b*C_v)))
        {
          // update v to v != 0
          v(l) = eta_times_u(l);
          // update C_v, C_b
          C_b = C_b + pow(beta(l), 2);
          C_v = C_v + 1;
          loop = TRUE;
        }else{
          // remains the same
          v(l) = 0;
        }
      }else{
        if ( pow(eta_times_u(l), 2) <= gamma*( sqrt(C_b)*sqrt(C_v) - sqrt(C_b - pow(beta(l), 2))*sqrt(C_v - 1) ) )
        {
          // update v to v = 0
          v(l) = 0;
          // update C_v, C_b
          C_b = C_b -  pow(beta(l), 2);
          C_v = C_v - 1;
          loop = TRUE;
        }else{
          //remains the same
          v(l) = eta_times_u(l);
        }
      }
    }
    if(num_iter++ > MAX_ITER_OPT){
      break;
    }
  }


}

void glasp::update_Eta_k(arma::vec & beta_k, int k){
  int p_k = grp_len(k);
  for (int j = 0; j < p_k; j++)
  {
    Eta.col( grp_start(k) + j) = X.col( grp_start(k) + j)*beta_k(j);
  }
}

void glasp::update_Eta(){
  for (int j = 0; j < p; j++)
  {
    Eta.col(j) = X.col(j)*beta(j);
  }

}

arma::vec glasp::predict_response(const arma::mat & newX){
  // Asume that newX in the same space as xL, xU
  arma::vec eta = newX*beta + intercept;
  return eta;
}

arma::vec glasp::predict_probability(const arma::mat & newX){
  arma::vec eta = predict_response(newX);
  arma::vec prob = 1/(1 + arma::exp(-eta));
  return prob;
}

arma::vec glasp::predict_class(const arma::mat & newX){
  arma::vec prob = predict_probability(newX);
  prob.for_each([](double & p_i){
    if(p_i > 0.5){ // or not?
      p_i = 1;
    }else{
      p_i = 0;
    }});
  return prob;
}

arma::vec glasp::predict(const arma::mat & newX, int type){
  switch (type)
  {
  case TYPE_PREDICT_RESPONSE:
    return predict_response(newX);
    break;
  case TYPE_PREDICT_PROB:
    return predict_probability(newX);
    break;
  case TYPE_PREDICT_CLASS:
    return predict_class(newX);
    break;
  default:
    switch (loss)
    {
    case LOSS_LOGIT:
      return predict_probability(newX);
      break;
    default:
      return predict_response(newX);
    break;
    }
  break;
  }
}

// Expose class glasp
RCPP_MODULE(Rcpp_glasp_export){
    Rcpp::class_<glasp>("glasp")

    .constructor<const Rcpp::List , int>()

    .method("fit", &glasp::glasp_optimization, "Fit model")
    .property("beta", &glasp::get_beta)
    .property("intercept", &glasp::get_intercept)
    .property("clusters", &glasp::get_clusters)
    ;
}

