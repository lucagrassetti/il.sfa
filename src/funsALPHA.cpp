#include <Rcpp.h>
using namespace Rcpp;
#include <cmath>


// [[Rcpp::export]]
double logh(double u, int m, DoubleVector y, DoubleVector eta, double lambda, double sigma)
{
  double f = 0.0;
  for (int j=0; j<m; j++){
    double rit = (y[j] - eta[j] - u) / sigma;
    double arg = (-1.0) * lambda * rit ;
    double Phi = R::pnorm(arg, 0.0, 1.0, 1, 1);
    f += Phi - pow(rit, 2) * 0.5;
  }
  return(f);
}


// [[Rcpp::export]]
DoubleVector getEffectsHN(double lnlambda, double lnsigma, Rcpp::List list_eta, Rcpp::List list_y, int niter, DoubleVector uinit)
{
  int len = list_y.size();
  double sigma = exp(lnsigma);
  double lambda = exp(lnlambda);
  double sigma2 = pow(sigma, 2.0);
  double lambda2 = pow(lambda, 2.0);
  double sigmau = exp(lnsigma + lnlambda) / sqrt(1 + lambda2);
  DoubleVector  out (len);
  // first locate the modes in u
  for(int i=0; i<len; i++){
    DoubleVector etai = list_eta[i];
    DoubleVector yi = list_y[i];
    int ni = yi.size();
    double ui = uinit[i] + sigmau * sqrt(2 / M_PI);
    double H = 0.0;
    for (int k=0; k<niter; k++){
      double g = 0.0;
      double f = 0.0;
      H = 0.0;
      for (int j=0; j<ni; j++){
        double rit = (yi[j] - etai[j] - ui) / sigma;
        double arg = (-1.0) * rit * lambda;
        double phi = R::dnorm(arg, 0.0, 1.0, 1);
        double Phi = R::pnorm(arg, 0.0, 1.0, 1, 1);
        double M = exp(phi - Phi);
        f += Phi - pow(rit, 2) * 0.5;
        g += 1/sigma * (rit + lambda * M);
        H += 1/sigma2 + lambda2 / sigma2 * (arg * M + pow(M, 2.0));
      }
      double db = g / H;
      double uiatt = ui + db;
      double ln = logh(uiatt, ni, yi, etai, lambda, sigma);
      while (ln < f + 0.001 * db * g){
        db *= 0.5;
        uiatt = ui + db;
        ln = logh(uiatt, ni, yi, etai, lambda, sigma);
      }
      ui += db;
    }
     out[i] = ui;
  }
  return(out);
}


// [[Rcpp::export]]
double likHN(double lnlambda, double lnsigma,  Rcpp::List list_eta, Rcpp::List list_y, int niter, DoubleVector uinit,
              DoubleVector ws, DoubleVector z)
{
  int len = list_y.size();
  double sigma = exp(lnsigma);
  double lambda = exp(lnlambda);
  double sigma2 = pow(sigma, 2.0);
  double lambda2 = pow(lambda, 2.0);
  double sigmau = exp(lnsigma + lnlambda) / sqrt(1 + lambda2);
  double ll = 0.0;
  // first locate the modes in u
  for(int i=0; i<len; i++){
    DoubleVector etai = list_eta[i];
    DoubleVector yi = list_y[i];
    int ni = yi.size();
    double H = 0.0;
    double ui = uinit[i] + sigmau * sqrt(2 / M_PI);
    double se = 0.0;
    for (int k=0; k<niter; k++){
      double g = 0.0;
      double f = 0.0;
      H = 0.0;
      for (int j=0; j<ni; j++){
        double rit = (yi[j] - etai[j] - ui) / sigma;
        double arg = (-1.0) * rit * lambda;
        double phi = R::dnorm(arg, 0.0, 1.0, 1);
        double Phi = R::pnorm(arg, 0.0, 1.0, 1, 1);
        double M = exp(phi - Phi);
        f += Phi - pow(rit, 2) * 0.5;
        g += 1/sigma * (rit + lambda * M);
        H += 1/sigma2 + lambda2 / sigma2 * (arg * M + pow(M, 2.0));
      }
      double db = g / H;
      double uiatt = ui + db;
      double ln = logh(uiatt, ni, yi, etai, lambda, sigma);
      while (ln < f + 0.001 * db * g){
        db *= 0.5;
        uiatt = ui + db;
        ln = logh(uiatt, ni, yi, etai, lambda, sigma);
      }
      ui += db;
      se = 1.0 / sqrt(H);
    }
    int nq = ws.size();
    // now computes AGH
    double Ii = 0.0;
    for (int k=0; k<nq; k++){
      double zk = sqrt(2) * se * z[k] + ui;
      double f = 0.0;
      for (int j=0; j<ni; j++){
        double rit = (yi[j] - etai[j] - zk) / sigma;
        double arg = (-1.0) * rit * lambda;
        f -= log(sigma);
        f += R::dnorm(rit, 0.0, 1.0, 1);
        f += R::pnorm(arg, 0.0, 1.0, 1, 1);
      }
     Ii += ws[k] * exp(f);
    }
    Ii *= se * sqrt(2.0);
    ll += log(Ii);
  }
  return(ll);
}

// [[Rcpp::export]]
double logf(double u, int m, DoubleVector y, DoubleVector eta, double sigmav, DoubleVector sigmau)
{
  double f = 0.0;
  for (int j=0; j<m; j++){
    double eit = y[j] - eta[j] - u;
    double arg = (-1.0) * eit / sigmav - sigmav / sigmau[j];
    double Phi = R::pnorm(arg, 0.0, 1.0, 1, 1);
    f += Phi + eit /sigmau[j];
  }
  return(f);
}

// [[Rcpp::export]]
double likEHET(DoubleVector gamma, double lnsigmav, Rcpp::List list_eta, Rcpp::List list_y, Rcpp::List list_z, int niter,
               DoubleVector uinit, DoubleVector ws, DoubleVector z)
{
  int len = list_y.size();
  double sigmav = exp(lnsigmav);
  double sigmav2 = exp(lnsigmav * 2.0);
  double ll = 0.0;
  // first locate the modes in u
  for(int i=0; i<len; i++){
    DoubleVector etai = list_eta[i];
    DoubleVector yi = list_y[i];
    DoubleVector zi = list_z[i];
    int ni = yi.size();
    double ui = uinit[i];
    double sigmaij;
    double se = 0.0;
    double H = 0.0;
    DoubleVector sigmai (ni, 0.0);
    for (int k=0; k<niter; k++){
      double g = 0;
      double f = 0.0;
      H = 0.0;
      for (int j=0; j<ni; j++){
        sigmaij = exp(gamma[0] + gamma[1] * zi[j]);
        if(k==0) sigmai[j] = sigmaij;
        double eit = yi[j] - etai[j] - ui;
        double arg = (-1.0) * eit / sigmav - sigmav / sigmaij;
        double phi = R::dnorm(arg, 0.0, 1.0, 1);
        double Phi = R::pnorm(arg, 0.0, 1.0, 1, 1);
        double M = exp(phi - Phi);
        f += Phi + eit /sigmaij;
        g += 1.0 / sigmav * M - 1.0 / sigmaij;
        H += 1.0 / sigmav2 * (arg * M + pow(M, 2.0));
      }
      double db = g / H;
      double uiatt = ui + db;
      double ln = logf(uiatt, ni, yi, etai, sigmav, sigmai);
      while (ln < f + 0.001 * db * g){
        db *= 0.5;
        uiatt = ui + db;
        ln = logf(uiatt, ni, yi, etai, sigmav, sigmai);
      }
    ui += db;
    se = 1.0 / sqrt(H);
    }
    int nq = ws.size();
    // now computes AGH
    double Ii = 0.0;
    for (int k=0; k<nq; k++){
      double zk = sqrt(2) * se * z[k] + ui;
      double f = 0.0;
      for (int j=0; j<ni; j++){
        double eit = yi[j] - etai[j] - zk;
        double arg = (-1.0) * eit / sigmav - sigmav / sigmai[j];
        f -= (gamma[0] + gamma[1] * zi[j]);
        f += sigmav2 * 0.5 / pow(sigmai[j], 2.0);
        f += R::pnorm(arg, 0.0, 1.0, 1, 1) + eit / sigmai[j];
      }
      Ii += ws[k] * exp(f);
    }
    Ii *= se * sqrt(2.0);
    ll += log(Ii);
  }
  return(ll);
}

// [[Rcpp::export]]
DoubleVector getEffectsEHET(DoubleVector gamma, double lnsigmav, Rcpp::List list_eta, Rcpp::List list_y,
                            Rcpp::List list_z, int niter, DoubleVector uinit)
{
  int len = list_y.size();
  double sigmav = exp(lnsigmav);
  double sigmav2 = exp(lnsigmav * 2.0);
  DoubleVector  out (len);
  // first locate the modes in u
  for(int i=0; i<len; i++){
    DoubleVector etai = list_eta[i];
    DoubleVector yi = list_y[i];
    DoubleVector zi = list_z[i];
    int ni = yi.size();
    double ui = uinit[i];
    double sigmaij;
    double H = 0.0;
    DoubleVector sigmai (ni, 0.0);
    for (int k=0; k<niter; k++){
      double g = 0;
      double f = 0.0;
      H = 0.0;
      for (int j=0; j<ni; j++){
        sigmaij = exp(gamma[0] + gamma[1] * zi[j]);
        if(k==0) sigmai[j] = sigmaij;
        double eit = yi[j] - etai[j] - ui;
        double arg = (-1.0) * eit / sigmav - sigmav / sigmaij;
        double phi = R::dnorm(arg, 0.0, 1.0, 1);
        double Phi = R::pnorm(arg, 0.0, 1.0, 1, 1);
        double M = exp(phi - Phi);
        f += Phi + eit /sigmaij;
        g += 1.0 / sigmav * M - 1.0 / sigmaij;
        H += 1.0 / sigmav2 * (arg * M + pow(M, 2.0));
      }
      double db = g / H;
      double uiatt = ui + db;
      double ln = logf(uiatt, ni, yi, etai, sigmav, sigmai);
      while (ln < f + 0.001 * db * g){
        db *= 0.5;
        uiatt = ui + db;
        ln = logf(uiatt, ni, yi, etai, sigmav, sigmai);
      }
      ui += db;
    }
    out[i] = ui;
  }
  return(out);
}



