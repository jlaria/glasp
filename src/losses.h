// losses.h

//#include "RcppArmadillo.h"
#include<string> // for string class

#define LOSS_LINEAR 0
#define LOSS_LOGIT 1
#define LOSS_COX 2



double R_linear(const arma::vec & beta,
                const arma::mat & X,
                const arma::mat & y){
    arma::vec r = y - X*beta;
    return arma::sum(arma::square(r))/(2*y.n_elem);
}
double R_linear2(const arma::vec & beta,
                const arma::mat & X,
                const arma::mat & eta_minus_k,
                const arma::mat & y){
    arma::vec r = y - X*beta - eta_minus_k;
    return arma::sum(arma::square(r))/(2*X.n_rows);
}

arma::vec grad_R_linear(const arma::vec & beta,
                        const arma::mat & X,
                        const arma::mat & y){
    return -X.t()*(y - X*beta)/y.n_elem;
}

arma::vec grad_R_linear2(const arma::vec & beta,
                        const arma::mat & X,
                        const arma::mat & eta_minus_k,
                        const arma::mat & y){
    return -X.t()*(y - X*beta - eta_minus_k)/X.n_elem;
}

// Logit
double log1pexp(double eta){
    if(eta <= 18 && eta >= -37){
        return std::log(1 + std::exp(eta));
    }
    if(eta > 18 && eta <= 33.3){
        return eta + std::exp(-eta);
    }
    if(eta > 33.3){
        return eta;
    }
    return exp(eta);
}

double R_logit(const arma::vec & beta,
               const arma::mat & X,
               const arma::mat & y){
    arma::vec eta = X*beta;
    // return arma::mean(arma::log(1 + arma::exp(eta)) - y % eta);

    // There is a more efficient and stable computation
    double rlogit = 0;
    for (int i = 0; i < eta.n_elem; i++)
    {
        rlogit += eta[i] = log1pexp(eta[i]) - y[i] * eta[i];
    }
    return rlogit/eta.n_elem;
}
double R_logit2(const arma::vec & beta,
                const arma::mat & X,
                const arma::mat & eta_minus_k,
                const arma::mat & y){
    arma::vec eta = X*beta + eta_minus_k;
    // return arma::mean(arma::log(1 + arma::exp(eta)) - y % eta);

    // There is a more efficient and stable computation
    double rlogit = 0;
    for (int i = 0; i < eta.n_elem; i++)
    {
        rlogit += eta[i] = log1pexp(eta[i]) - y[i] * eta[i];
    }
    return rlogit/eta.n_elem;
}

double inv1pexpm(double eta){
    if(eta >= -30){
        return 1/(1 + std::exp(-eta));
    }else{
        return std::exp(eta);
    }
}
arma::vec grad_R_logit(const arma::vec & beta,
                       const arma::mat & X,
                       const arma::mat & y){
    // arma::vec eta = X*beta;
    // arma::vec d = 1/(1 + arma::exp(-eta)) - y;
    // return X.t() * d/y.n_elem;

    // More stable
    arma::vec tmp = X*beta;
    for (int i = 0; i < tmp.n_elem; i++)
    {
        tmp[i] = inv1pexpm(tmp[i]) - y[i];
    }
    return X.t()*tmp/y.n_elem;
}
arma::vec grad_R_logit2(const arma::vec & beta,
                        const arma::mat & X,
                        const arma::mat & eta_minus_k,
                        const arma::mat & y){
    // More stable
    arma::vec tmp = X*beta + eta_minus_k;
    for (int i = 0; i < tmp.n_elem; i++)
    {
        tmp[i] = inv1pexpm(tmp[i]) - y[i];
    }
    return X.t()*tmp/y.n_elem;
}

// ph model
arma::vec logCumSumExp(arma::vec & eta){
    // Find the maximum
    double a = eta.max();

    arma::vec eta_updated = eta - a;
    eta_updated = arma::exp(eta_updated);
    eta_updated = arma::cumsum(eta_updated);
    return a + arma::log(eta_updated);
}


double R_phcox(const arma::vec & beta, const arma::mat & X, const arma::mat & y){
    // We assume y(,0) is in decresing order
    arma::vec eta = X*beta;
    arma::vec log_S = logCumSumExp(eta);

    return arma::as_scalar(y.col(1).t()*(-eta + log_S))/X.n_rows;
}
double R_phcox2(const arma::vec & beta, const arma::mat & X, const arma::mat & eta_minus_k, const arma::mat & y){
    // We assume y(,0) is in decresing order
    arma::vec eta = X*beta + eta_minus_k;

    arma::vec log_S = logCumSumExp(eta);

    return arma::as_scalar(y.col(1).t()*(-eta + log_S))/X.n_rows;
}
arma::vec grad_R_phcox(const arma::vec & beta, const arma::mat & X, const arma::mat & y){

    arma::vec eta = X*beta;
    double a = eta.max();

    arma::vec exp_eta_updated = arma::exp(eta - a);
    double sum_exp_eta_updated = 0.0;
    arma::rowvec sum_exp_eta_updated_times_x = arma::zeros(1, X.n_cols);
    arma::rowvec grad2 = arma::zeros(1, X.n_cols);

    for (size_t j = 0; j < X.n_rows; j++){
        /* code */
        sum_exp_eta_updated += exp_eta_updated(j);
        sum_exp_eta_updated_times_x += exp_eta_updated(j)*X.row(j);

        if(y(j,1) != 0){
            grad2 += sum_exp_eta_updated_times_x/sum_exp_eta_updated - X.row(j);
        }
    }
    return grad2.t()/X.n_rows;
}

arma::vec grad_R_phcox2(const arma::vec & beta, const arma::mat & X, const arma::mat & eta_minus_k, const arma::mat & y){

    arma::vec eta = X*beta + eta_minus_k;
    double a = eta.max();

    arma::vec exp_eta_updated = arma::exp(eta - a);
    double sum_exp_eta_updated = 0.0;
    arma::rowvec sum_exp_eta_updated_times_x = arma::zeros(1, X.n_cols);
    arma::rowvec grad2 = arma::zeros(1, X.n_cols);

    for (size_t j = 0; j < X.n_rows; j++){
        /* code */
        sum_exp_eta_updated += exp_eta_updated(j);
        sum_exp_eta_updated_times_x += exp_eta_updated(j)*X.row(j);

        if(y(j,1) != 0){
            grad2 += sum_exp_eta_updated_times_x/sum_exp_eta_updated - X.row(j);
        }
    }
    return grad2.t()/X.n_rows;
}


/* // COx

double R_cox(const arma::vec & beta,
               const arma::mat & X,
               const arma::mat & y){
    // IMPORTANT! y: is N x 2 with y(,0) = t and y(,1) = sigma
  // In addition, t is in decreasing order.

  arma::vec eta =  X*beta;
  int N = X.n_rows;

  // Compute S1
  arma::vec S1 = arma::cumsum(arma::exp(eta));

  for(int i = N-1; i>0; i--){
    if(y(i,0) == y(i-1, 0)){
      S1(i-1) = S1(i);
    }
  }

  return arma::as_scalar(y.col(1).t()*(-eta + arma::log(S1)))/X.n_rows;
}
double R_cox2(const arma::vec & beta,
                const arma::mat & X,
                const arma::mat & eta_minus_k,
                const arma::mat & y){
    // IMPORTANT! y: is N x 2 with y(,0) = t and y(,1) = sigma
  // In addition, t is in decreasing order.

  arma::vec eta =  X*beta + eta_minus_k;
  int N = X.n_rows;

  // Compute S1
  arma::vec S1 = arma::cumsum(arma::exp(eta));

//   for(int i = N-1; i>0; i--){
//     if(y(i,0) == y(i-1, 0)){
//       S1(i-1) = S1(i);
//     }
//   }

  return arma::as_scalar(y.col(1).t()*(-eta + arma::log(S1)))/X.n_rows;
}

arma::vec grad_R_cox(const arma::vec & beta,
                       const arma::mat & X,
                       const arma::mat & y){
  // IMPORTANT! y: is N x 2 with y(,0) = t and y(,1) = sigma
  // In addition, t is in decreasing order.

  arma::vec eta = X*beta;
  int N = X.n_rows;
  int p_k = X.n_cols;

  // Compute S1
  arma::vec S1 = arma::exp(eta);

  // Compute S2
  arma::mat S2 = X.t();
  for(int i = 0; i < N; i++){
    S2.col(i) *= S1(i);
  }
  S2 = arma::cumsum(S2, 1);

  // Update S1
  S1 = arma::cumsum(S1);

  // Compute S3
  arma::colvec S3 = arma::zeros<arma::colvec>(p_k);

  int i = N-1;
  while(i>=0){
    if(y(i,1)==1){S3 += S2.col(i)/S1(i);}
    i--;
  }
  arma::vec g = (X.t()*y.col(1) - S3);
  return -g/X.n_rows;
}
arma::vec grad_R_cox2(const arma::vec & beta,
                        const arma::mat & X,
                        const arma::mat & eta_minus_k,
                        const arma::mat & y){
  // IMPORTANT! y: is N x 2 with y(,0) = t and y(,1) = sigma
  // In addition, t is in decreasing order.

  arma::vec eta = X*beta + eta_minus_k;
  int N = X.n_rows;
  int p_k = X.n_cols;

  // Compute S1
  arma::vec S1 = arma::exp(eta);

  // Compute S2
  arma::mat S2 = X.t();
  for(int i = 0; i < N; i++){
    S2.col(i) *= S1(i);
  }
  S2 = arma::cumsum(S2, 1);

  // Update S1
  S1 = arma::cumsum(S1);

  // Compute S3
  arma::colvec S3 = arma::zeros<arma::colvec>(p_k);

  int i = N-1;
  while(i>=0){
    if(y(i,1)==1){S3 += S2.col(i)/S1(i);}
    i--;
  }
  arma::vec g = (X.t()*y.col(1) - S3);
  return -g/X.n_rows;
} */

//===============================
// Class: Loss

class Loss
{
public:
    double (*R)(const arma::vec & beta, const arma::mat & X, const arma::mat & y);
    arma::vec (*gradR)(const arma::vec & beta, const arma::mat & X, const arma::mat & y);
    double (*R2)(const arma::vec & beta, const arma::mat & X,const arma::mat & eta_minus_k, const arma::mat & y);
    arma::vec (*gradR2)(const arma::vec & beta, const arma::mat & X, const arma::mat & eta_minus_k, const arma::mat & y);

    std::string name;

    Loss(){
        name = "linear";
        R = R_linear;
        gradR = grad_R_linear;
        R2 = R_linear2;
        gradR2 = grad_R_linear2;

    }

    Loss(int type){
        switch (type)
        {
        case LOSS_LOGIT:
            this->name = "logit";
            this->R = R_logit;
            this->gradR = grad_R_logit;
            this->R2 = R_logit2;
            this->gradR2 = grad_R_logit2;
            break;
        case LOSS_COX:
            this->name = "cox";
            this->R = R_phcox;
            this->gradR = grad_R_phcox;
            this->R2 = R_phcox2;
            this->gradR2 = grad_R_phcox2;
            break;
        default:
            this->name = "linear";
            this->R = R_linear;
            this->gradR = grad_R_linear;
            this->R2 = R_linear2;
            this->gradR2 = grad_R_linear2;
            break;
            }
    }
};

// This class is not exposed...


//' Mean squared error
//'
//' Computes the mean squared error of two vectors
//'
//' @param ypred vector of predictions
//' @param ytrue vector of true labels
//' @export
// [[Rcpp::export]]
double mean_squared_error(arma::vec ypred, arma::vec ytrue){
    return arma::mean(arma::square(ypred - ytrue));
}

//' Binary cross-entropy loss
//'
//' Computes the binary cross-entropy loss (log-loss) of two vectors
//'
//' @param ypred vector of predicted probabilities
//' @param ytrue vector of true labels {0, 1}
//' @export
// [[Rcpp::export]]
double binary_cross_entropy(arma::vec ypred, arma::vec ytrue){
    unsigned int N = ypred.n_elem;
    double binXent = 0;
    for (unsigned int i = 0; i < N; i++)
    {
        if (ypred(i) < 1e-8){ypred(i) = 1e-8;}
        if (ypred(i) > 1.0 - 1e-8){ypred(i) = 1.0 - 1e-8;}

        binXent -= (ytrue(i)*std::log(ypred(i)) + (1 - ytrue(i))*std::log(1 - ypred(i)));
    }
    return binXent/N;
}
