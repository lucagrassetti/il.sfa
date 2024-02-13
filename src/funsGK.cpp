// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
#include <Rcpp.h>
using namespace Rcpp;
#include <cmath> // std::pow


class IntegGN: public Numer::Func {
private:
  double mu;
  double sigma;
  double alpha;
  double lambda;
  double y;
public:
  IntegGN(double mu_, double sigma_, double alpha_, double lambda_, double y_) : mu(mu_), sigma(sigma_), alpha(alpha_), lambda(lambda_), y(y_) {}

  double operator()(const double& x) const
  {
    double som = 0;
    som += R::dnorm(y, mu - x,  sigma, 1) + R::dgamma(x, alpha, lambda, 1);
    return exp(som);
  }
};

// for gradient wrt mu
class IntegGN1: public Numer::Func {
private:
  double mu;
  double sigma;
  double alpha;
  double lambda;
  double y;
public:
  IntegGN1(double mu_, double sigma_, double alpha_, double lambda_, double y_) : mu(mu_), sigma(sigma_), alpha(alpha_), lambda(lambda_), y(y_) {}

  double operator()(const double& x) const
  {
    double som = R::dnorm(y, mu - x,  sigma, 1) + R::dgamma(x, alpha, lambda, 1);
    return exp(som) * (y - mu + x) / pow(sigma, 2);
  }
};



// for hessian wrt mu
class IntegGN2: public Numer::Func {
private:
  double mu;
  double sigma;
  double alpha;
  double lambda;
  double y;
public:
  IntegGN2(double mu_, double sigma_, double alpha_, double lambda_, double y_) : mu(mu_), sigma(sigma_), alpha(alpha_), lambda(lambda_), y(y_) {}

  double operator()(const double& x) const
  {
    double som = R::dnorm(y, mu - x,  sigma, 1) + R::dgamma(x, alpha, lambda, 1);
    return exp(som) * (-pow(sigma, -2) + pow(y - mu + x, 2) / pow(sigma, 4));
  }
};


// [[Rcpp::plugins(cpp11)]]
template<class T> class trans_func: public T {
public:
  using T::T;

  double operator()(const double& t) const {
    double x = (1-t)/t;
    return (T::operator()(x))/pow(t, 2);
  }
};



// [[Rcpp::export]]
double dnormgam(double y, double mu, double sigma, double alpha, double lambda)
{
  trans_func<IntegGN> f(mu, sigma, alpha, lambda, y);
  double err_est;
  int err_code;
  int subdiv = 100;
  double eps_abs = 0.00000001; //maybe even too small
  double eps_rel = 0.000001;   //same here
  const Numer::Integrator<double>::QuadratureRule rule = Numer::Integrator<double>::GaussKronrod61;
  const double result = Numer::integrate(f, 0, 1, err_est, err_code, subdiv, eps_abs, eps_rel,
                                         rule);
  return log(result);
}



// [[Rcpp::export]]
DoubleVector Ii12(double y, double mu, double sigma, double alpha, double lambda)
{
  trans_func<IntegGN1> g(mu, sigma, alpha, lambda, y);
  trans_func<IntegGN2> h(mu, sigma, alpha, lambda, y);
  double err_est;
  int err_code;
  int subdiv = 100;
  double eps_abs = 0.00000001; //maybe even too small
  double eps_rel = 0.000001;   //same here
  const Numer::Integrator<double>::QuadratureRule rule = Numer::Integrator<double>::GaussKronrod61;
  const double result1 = Numer::integrate(g, 0, 1, err_est, err_code, subdiv, eps_abs, eps_rel,
                                         rule);
  const double result2 = Numer::integrate(h, 0, 1, err_est, err_code, subdiv, eps_abs, eps_rel,
                                          rule);
  DoubleVector  out (2);
  out[0] = result1;
  out[1] = result2;
  return out;
}


// [[Rcpp::export]]
double logh_G(double u, int m, DoubleVector y, DoubleVector eta, double sigma, double alpha, double lambda)
{
  double f = 0.0;
  for (int j=0; j<m; j++){
    double mu = u + eta[j];
    double fmu = dnormgam(y[j], mu, sigma, alpha, lambda);
    f += fmu;
  }
  return(f);
}





// [[Rcpp::export]]
double logHess(double u, int m, DoubleVector y, DoubleVector eta, double sigma, double alpha, double lambda)
{
  double H = 0.0;
  for (int j=0; j<m; j++){
    double mu = u + eta[j];
    double fmu = dnormgam(y[j], mu, sigma, alpha, lambda);
    double Ij = exp(fmu);
    DoubleVector I12 = Ii12(y[j], mu, sigma, alpha, lambda);
    H -= (I12[1] * Ij - pow(I12[0], 2)) / pow(Ij, 2);
  }
  return(H);
}


// [[Rcpp::export]]
DoubleVector getEffectsG(double lnsigma, double lnalpha, double lnlambda,
                          Rcpp::List list_eta, Rcpp::List list_y, int niter,
                          DoubleVector uinit)
{
  int len = list_y.size();
  double alpha = exp(lnalpha);
  double lambda = exp(lnlambda);
  double sigma = exp(lnsigma);
  DoubleVector  out (len);
  // first locate the modes in u
  for(int i=0; i<len; i++){
    DoubleVector etai = list_eta[i];
    DoubleVector yi = list_y[i];
    int ni = yi.size();
    double ui = uinit[i];
    double H = 0.0;
    for (int k=0; k<niter; k++){
      double g = 0.0;
      double f = 0.0;
      H = 0.0;
      for (int j=0; j<ni; j++){
        double mu = ui + etai[j];
        double fmu = dnormgam(yi[j], mu, sigma, alpha, lambda);
        double Ij = exp(fmu);
        DoubleVector I12 = Ii12(yi[j], mu, sigma, alpha, lambda);
        f += fmu;
        g += I12[0] / Ij;
        H -= (I12[1] * Ij - pow(I12[0], 2)) / pow(Ij, 2);
      }
      double db = g / H;
      double uiatt = ui + db;
      double ln = logh_G(uiatt, ni, yi, etai, sigma, alpha, lambda);
      while (ln < f + 0.001 * db * g){
        db *= 0.5;
        uiatt = ui + db;
        ln = logh_G(uiatt, ni, yi, etai, sigma, alpha, lambda);}
       ui += db;}
    out[i] = ui;
  }
  return(out);
}


// [[Rcpp::export]]
  double minlogh(double u, int m, DoubleVector y, DoubleVector eta, double sigma, double alpha, double lambda)
  {
    double maxk = pow(10, 6);
    double f = 0.0;
    for (int j=0; j<m; j++){
      double mu = u + eta[j];
      double fmu = dnormgam(y[j], mu, sigma, alpha, lambda);
      f -= fmu;
    }
   if(f == R_NegInf) f = (-1.0) * maxk;
   if(f == R_PosInf) f = maxk;
    return(f);
  }


// [[Rcpp::export]]
DoubleVector getEffectsGS(double lnsigma, double lnalpha, double lnlambda,
                             Rcpp::List list_eta, Rcpp::List list_y,
                             DoubleVector ulow, DoubleVector uup, double eps){
      int len = list_y.size();
      double alpha = exp(lnalpha);
      double lambda = exp(lnlambda);
      double sigma = exp(lnsigma);
      DoubleVector  out (len);
      double gc = (3.0 - sqrt(5.0)) * 0.5;
      // the value of gc is taken from Gray
      for(int i=0; i<len; i++){
        DoubleVector etai = list_eta[i];
        DoubleVector yi = list_y[i];
        int ni = yi.size();
        double xl = ulow[i];
        double xu = uup[i];
        double tmp = gc * (xu - xl);
        double xml = xl + tmp;
        double xmu = xu - tmp;
        double fl = minlogh(xml, ni, yi, etai, sigma, alpha, lambda);
        double fu = minlogh(xmu, ni, yi, etai, sigma, alpha, lambda);
        while (abs(xu - xl)>(0.00001 + abs(xl)) * eps){
            if(fl < fu){
              xu = xmu;
              xmu = xml;
              fu = fl;
              xml = xl + gc * (xu - xl);
              fl = minlogh(xml, ni, yi, etai, sigma, alpha, lambda);
            } else {
              xl = xml;
              xml = xmu;
              fl = fu;
              xmu = xu - gc * (xu - xl);
              fu = minlogh(xmu, ni, yi, etai, sigma, alpha, lambda);
              }
            }
        if (fl < fu) {out[i] = xml;} else {out[i] = xmu;}
        }
       return(out);
}


// [[Rcpp::export]]
double likG(double lnsigma, double lnalpha, double lnlambda,
            Rcpp::List list_eta, Rcpp::List list_y,  DoubleVector ulow,
            DoubleVector uup, double eps, DoubleVector ws, DoubleVector z,
            int niter, DoubleVector uinit, int method)
{
  int len = list_y.size();
  double alpha = exp(lnalpha);
  double lambda = exp(lnlambda);
  double sigma = exp(lnsigma);
  double ll = 0.0;
  DoubleVector  u (len);
  if(method == 0) {u = getEffectsGS(lnsigma, lnalpha, lnlambda, list_eta, list_y,
                   ulow, uup, eps);}
  else {u = getEffectsG(lnsigma, lnalpha, lnlambda, list_eta, list_y, niter, uinit);}
  // now computes AGH
  int nq = ws.size();
  for(int i=0; i<len; i++){
    DoubleVector etai = list_eta[i];
    DoubleVector yi = list_y[i];
    int ni = yi.size();
    double Ii = 0.0;
    double H = logHess(u[i], ni, yi, etai, sigma, alpha, lambda);
    double se =  1.0 / sqrt(H);
    for (int k=0; k<nq; k++){
      double zk = sqrt(2) * se * z[k] + u[i];
      double f = 0.0;
      for (int j=0; j<ni; j++){
        double mu = zk + etai[j];
        double fmu = dnormgam(yi[j], mu, sigma, alpha, lambda);
        f += fmu;
      }
      Ii += ws[k] * exp(f);
    }
    Ii *= se * sqrt(2.0);
    ll += log(Ii);
  }
  return(ll);
}
