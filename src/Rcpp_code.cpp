#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

#include <iostream>
using namespace Rcpp;
using namespace arma;

int lr_(const arma::mat& A, 
        const arma::mat& X, 
        const arma::colvec& y,
        double tau,
        arma::colvec& u,
        arma::colvec& v,
        arma::colvec& beta,
        arma::colvec& u_distance,
        double& d,
        double l, 
        double sigma,
        double tol,
        int max_iter,
        int verbose,
        bool uv_is_init);

int irlba (
    arma::mat& U, 
    arma::vec& s, 
    arma::mat& V,  
    const arma::mat& A,
    int k);

double prox(double tau,
            double xi,
            double alpha);

int update_z(const arma::colvec& y,
             const arma::mat& X,
             const arma::colvec& u,
             const arma::colvec& v,
             const arma::colvec& beta,
             const arma::colvec& theta,
             double sigma,
             double tau,
             arma::colvec& z);

int update_beta(const arma::colvec& y,
                const arma::mat& X,
                const arma::colvec& u,
                const arma::colvec& v,
                arma::colvec& beta,
                const arma::colvec& theta,
                double sigma,
                const arma::colvec& z);

double update_d(const arma::mat& A,
                const arma::colvec& u, 
                const arma::colvec& v) ;

int update_u(const arma::mat& A,
             const arma::colvec& y,
             const arma::mat& X,
             arma::colvec& u,
             const arma::colvec& v,
             const arma::colvec& beta,
             const arma::colvec& theta,
             double d,
             double sigma,
             double tau,
             double l,
             const arma::colvec& z);

int update_v(const arma::mat& A,
             const arma::colvec& y,
             const arma::mat& X,
             const arma::colvec& u,
             arma::colvec& v,
             const arma::colvec& beta,
             const arma::colvec& theta,
             double d,
             double sigma,
             double tau,
             double l,
             const arma::colvec& z);

int update_theta(arma::colvec& theta,
                 double sigma,
                 const arma::mat& X,
                 const arma::colvec& u,
                 const arma::colvec& v,
                 const arma::colvec& y,
                 const arma::colvec& z,
                 const arma::colvec& beta);

double spec_norm(const arma::mat& A) ;

bool stop_condition(const arma::colvec& u, 
                    const arma::colvec& u_old,
                    const arma::colvec& v, 
                    const arma::colvec& v_old,
                    double tol);

// [[Rcpp::depends(RcppArmadillo)]]
  
// [[Rcpp::export]]
Rcpp::List lr(const arma::mat& A,
              const arma::mat& X,
              const arma::colvec& y,
              double tau,
              double l,
              double sigma = 1,
              double tol = 1e-4,
              int max_iter = 1000,
              int verbose = 0,
              int scaled = 1){
  
  int n_obs = X.n_rows, n = A.n_rows, k = X.n_cols;
  
  // output
  

 
  Rcpp::List ret_list;
  arma::colvec u, v, beta(k + 2, fill::zeros), residuals(n), u_distance(max_iter);
  double d;
  int iter;
  Rcpp::StringVector method = "lr";
  
  iter = lr_(A, X, y, tau, u, v, beta, u_distance, d, l, sigma, tol, max_iter, verbose, false);
  residuals = y - join_rows(X, u.head(n_obs), v.head(n_obs)) * beta;
  
  ret_list = Rcpp::List::create(_("u") = u,
                                _("v") = v,
                                _("beta") = beta,
                                _("d") = d,
                                _("l") = l,
                                _("iter") = iter,
                                _("residuals") = residuals,
                                _("method") = method,
                                _("max_iter") = max_iter,
                                _("u_distance") = u_distance.head(iter));
  
  return ret_list;
}


// [[Rcpp::export]]
int lr_(const arma::mat& A, 
        const arma::mat& X, 
        const arma::colvec& y,
        double tau,
        arma::colvec& u,
        arma::colvec& v,
        arma::colvec& beta,
        arma::colvec& u_distance,
        double& d,
        double l, 
        double sigma,
        double tol,
        int max_iter,
        int verbose,
        bool uv_is_init){
  int n_obs = X.n_rows, n = A.n_rows;
  
  // aux
  bool cond = 1; 
  arma::colvec u_old, v_old, beta_old, z(n_obs), z_old, theta, theta_old;
  double d_old = 0;
  int iter = 0;
  
  // initiate 
  if(!uv_is_init) {
    arma::mat U, V;
    arma::vec s;
    irlba(U, s, V, A, 1);
    u = U.col(0)*sqrt(n); v = V.col(0)*sqrt(n); d = s(0)/n;
  }

  theta = arma::randn<arma::colvec>(n);
  arma::mat WW = arma::join_rows(X, u.head(n_obs));
  arma::mat W = arma::join_rows(WW, v.head(n_obs));
  beta = arma::solve(W, y);  
  update_z(y,X,u,v,beta,theta,sigma,tau,z);
  
  
  
  while(cond & (iter < max_iter)) {
    u_old = u; v_old = v; beta_old = beta; d_old = d; z_old = z; theta_old = theta;
    
    
    
    update_beta(y, X, u, v, beta, theta, sigma, z);
    update_z(y,X,u,v,beta,theta,sigma,tau,z);
    d = update_d(A, u, v);
    d /= pow(n, 2);
    // if(d < 0) {u = -u; d = -d;}
    update_u(A, y, X, u, v, beta, theta, d, sigma, tau, l, z);
    update_v(A, y, X, u, v, beta, theta, d, sigma, tau, l, z);
    update_theta(theta, sigma, X, u, v, y, z, beta);
    
    
    if(verbose) {
      Rcpp::Rcout << "\t\t" << iter << ": " << beta.t() << std::endl;
      if(verbose == 3) {
        Rcpp::Rcout << "\t\t" << "u: " << spec_norm(u_old * u_old.t() - u * u.t()) << "; v: " << spec_norm(v_old * v_old.t() - v * v.t()) << "; beta: " << norm(beta - beta_old) << "; d: " << sqrt(pow(d - d_old, 2)) << std::endl;
      }
    }
    
    arma::mat diff_u = (u_old * u_old.t() - u * u.t());
    u_distance(iter) = spec_norm(diff_u);
    
    cond = stop_condition(u, u_old, v, v_old, tol);
    
    iter++;
  }
  
  return iter;
}

int irlba (
    arma::mat& U, 
    arma::vec& s, 
    arma::mat& V,  
    const arma::mat& A,
    int k)
{
  
  // Obtain environment containing function
  Rcpp::Environment base("package:irlba"); 
  
  // Make function callable from C++
  Rcpp::Function svd_r = base["irlba"];    
  
  // Call the function and receive its list output
  Rcpp::List res = svd_r(Rcpp::_["A"] = A,
                         Rcpp::_["nu"]  = k,
                         Rcpp::_["nv"]  = k); 
  
  U = as<arma::mat>(res["u"]);
  V = as<arma::mat>(res["v"]);
  s = as<arma::vec>(res["d"]);
  
  return res["iter"];
}

double prox(double tau,
            double xi,
            double alpha){
  if(xi>tau/alpha){
    return xi-tau/alpha;
  }else if(xi<(tau-1)/alpha){
    return xi-(tau-1)/alpha;
  }else{
    return 0;
  }
}

int update_z(const arma::colvec& y,
             const arma::mat& X,
             const arma::colvec& u,
             const arma::colvec& v,
             const arma::colvec& beta,
             const arma::colvec& theta,
             double sigma,
             double tau,
             arma::colvec& z){
  int n_obs = X.n_rows;
  arma::mat W = join_rows(X, u.head(n_obs), v.head(n_obs));
  arma::colvec xi = y -  W * beta + theta/sigma;
  for(int i = 0; i < n_obs; i++){
    z(i) = prox(tau, xi(i), n_obs * sigma);
  }
  return 1;
};

int update_beta(const arma::colvec& y,
                const arma::mat& X,
                const arma::colvec& u,
                const arma::colvec& v,
                arma::colvec& beta,
                const arma::colvec& theta,
                double sigma,
                const arma::colvec& z){
  int n_obs = X.n_rows;
  arma::mat WW = arma::join_rows(X, u.head(n_obs));
  arma::mat W = arma::join_rows(WW, v.head(n_obs));
  beta = arma::solve(W, y - z + theta / sigma);
  return 1;
}

double update_d(const arma::mat& A,
                const arma::colvec& u, 
                const arma::colvec& v) 
{
  double d = as_scalar(u.t() * A * v);
  return d;
}

int update_u(const arma::mat& A,
             const arma::colvec& y,
             const arma::mat& X,
             arma::colvec& u,
             const arma::colvec& v,
             const arma::colvec& beta,
             const arma::colvec& theta,
             double d,
             double sigma,
             double tau,
             double l,
             const arma::colvec& z){
  int n_obs = X.n_rows, k = X.n_cols, n = A.n_rows, n_c = n - n_obs;
  arma::colvec rmat = beta(k)*(sigma*(y - X*beta.head(k) - beta(k+1)*v.head(n_obs) - z) + theta) + 2*l*d*A.head_rows(n_obs)*v/n;
  u.head(n_obs) = rmat / (pow(beta(k),2)*sigma + 2*l*pow(d,2));
  
  if(n_obs < n) {
    u.tail(n_c) = A.tail_rows(n_c)*v/d/n;
  }
  
  // normalize u
  u = arma::normalise(u) * sqrt(n);
  
  return 1;
}

int update_v(const arma::mat& A,
             const arma::colvec& y,
             const arma::mat& X,
             const arma::colvec& u,
             arma::colvec& v,
             const arma::colvec& beta,
             const arma::colvec& theta,
             double d,
             double sigma,
             double tau,
             double l,
             const arma::colvec& z){
  int n_obs = X.n_rows, k = X.n_cols, n = A.n_rows, n_c = n - n_obs;
  arma::colvec rmat = beta(k+1)*(sigma*(y - X*beta.head(k) - beta(k)*u.head(n_obs) - z) + theta) + 2*l*d*A.head_cols(n_obs).t()*u/n;
  v.head(n_obs) = rmat / (pow(beta(k+1),2)*sigma + 2*l*pow(d,2));
  
  if(n_obs < n) {
    v.tail(n_c) = A.tail_cols(n_c).t()*u/d/n;
  }
  
  // normalize u
  v = arma::normalise(v) * sqrt(n);
  
  return 1;
}

int update_theta(arma::colvec& theta,
                 double sigma,
                 const arma::mat& X,
                 const arma::colvec& u,
                 const arma::colvec& v,
                 const arma::colvec& y,
                 const arma::colvec& z,
                 const arma::colvec& beta){
  int n_obs = X.n_rows;
  arma::mat WW = arma::join_rows(X, u.head(n_obs));
  arma::mat W = arma::join_rows(WW, v.head(n_obs));
  theta -= sigma*(W*beta+z-y);
  return 1;
}

double spec_norm (
    const arma::mat& A
) 
{
  int n = A.n_rows;
  arma::mat U_diff, V_diff;
  arma::vec s_diff;
  
  irlba(U_diff, s_diff, V_diff, A, 1);
  
  return s_diff[0]/n;
}

bool stop_condition(const arma::colvec& u, 
                    const arma::colvec& u_old,
                    const arma::colvec& v, 
                    const arma::colvec& v_old,
                    double tol)
{
  int N = u.n_elem;
  
  arma::mat diff_u = (u_old * u_old.t() - u * u.t());
  arma::mat diff_v = (v_old * v_old.t() - v * v.t());
  
  return ((spec_norm(diff_u) > tol) | (spec_norm(diff_v) > tol));
}




