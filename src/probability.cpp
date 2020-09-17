
#include <Rcpp.h>
#include <math.h>
#include "probability.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// draw from continuous uniform distribution on interval [0,1)
double runif_0_1() {
  return R::runif(0,1);
}

//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b)
double runif1(double a, double b) {
  return R::runif(a,b);
}

//------------------------------------------------
// draw from Bernoulli(p) distribution
bool rbernoulli1(double p) {
  return R::rbinom(1, p);
}

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(double mean, double sd) {
  return R::rnorm(mean, sd);
}

//------------------------------------------------
// density of univariate normal distribution
double dnorm1(double x, double mean, double sd, bool log_on) {
  return R::dnorm(x, mean, sd, log_on);
}

//------------------------------------------------
// density of log-normal distribution
double dlnorm1(double x, double meanlog, double sdlog, bool log_on) {
  return R::dlnorm(x, meanlog, sdlog, log_on);
}

//------------------------------------------------
// draw from univariate normal distribution and reflect to interval (a,b)
double rnorm1_interval(double mean, double sd, double a, double b) {

  // draw raw value relative to a
  double ret = rnorm1(mean, sd) - a;
  double interval_difference = b - a;
  
  // reflect off boundries at 0 and (b-a)
  if (ret < 0 || ret > interval_difference) {
    // use multiple reflections to bring into range [-(b-a), 2(b-a)]
    double modded = std::fmod(ret, 2*interval_difference);
    
    // use one more reflection to bring into range [0, (b-a)]
    if (modded < 0) {
      modded = -1*modded;
    }
    if (modded > interval_difference) {
      modded = 2*interval_difference - modded;
    }
    ret = modded;
  }

  // no longer relative to a
  ret += a;
  
  // don't let ret equal exactly a or b
  if (ret == a) {
    ret += UNDERFLO;
  } else if (ret == b) {
    ret -= UNDERFLO;
  }

  return ret;
}

//------------------------------------------------
// sample single value from given probability vector (that sums to pSum)
int sample1(vector<double> &p, double pSum) {
  double rand = pSum*runif_0_1();
  double z = 0;
  for (int i=0; i<int(p.size()); i++) {
    z += p[i];
    if (rand<z) {
      return i+1;
    }
  }
  return 0;
}

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal
// probability. Works on positive or negative values of a or b, and works
// irrespective of which of a or b is larger.
int sample2(int a, int b) {
  if (a<b) {
    return floor(runif1(a, b+1));
  } else {
    return floor(runif1(b, a+1));
  }
}

//------------------------------------------------
// sample a given number of values from a vector without replacement (templated
// for different data types). Note, this function re-arranges the original
// vector (passed in by reference), and the result is stored in the first n
// elements.
// sample3
// DEFINED IN HEADER

//------------------------------------------------
// draw from gamma(shape,rate) distribution
double rgamma1(double shape, double rate) {
  double x = R::rgamma(shape, 1/rate);
  
  // check for zero or infinite values (catches bug present in Visual Studio 2010)
  if (x<UNDERFLO) {
    x = UNDERFLO;
  }
  if (x>OVERFLO) {
    x = OVERFLO;
  }
  return x;
}

//------------------------------------------------
// density of gamma(shape,rate) distribution
double dgamma1(double x, double shape, double rate) {
  double y = R::dgamma(x, shape, 1/rate, FALSE);
  
  // check for zero or infinite values (catches bug present in Visual Studio 2010)
  if (y<UNDERFLO) {
    y = UNDERFLO;
  }
  if (y>OVERFLO) {
    y = OVERFLO;
  }
  return y;
}

//------------------------------------------------
// draw from beta(alpha,beta) distribution
double rbeta1(double shape1, double shape2) {
  if (shape1==1 && shape2==1) {
    return runif_0_1();
  }
  return R::rbeta(shape1, shape2);
}

//------------------------------------------------
// probability density of beta(shape1,shape2) distribution
double dbeta1(double x, double shape1, double shape2, bool return_log) {
  return R::dbeta(x, shape1, shape2, return_log);
}

//------------------------------------------------
// draw from dirichlet distribution using vector of shape parameters. Return vector of values.
vector<double> rdirichlet1(vector<double> &shapeVec) {
  // draw a series of gamma random variables
  int n = shapeVec.size();
  vector<double> ret(n);
  double retSum = 0;
  for (int i=0; i<n; i++) {
    ret[i] = rgamma1(shapeVec[i], 1.0);
    retSum += ret[i];
  }
  // divide all by the sum
  double retSumInv = 1.0/retSum;
  for (int i=0; i<n; i++) {
    ret[i] *= retSumInv;
  }
  return(ret);
}

//------------------------------------------------
// draw from dirichlet distribution using bespoke inputs. Outputs are stored in
// x, passed by reference for speed. Shape parameters are equal to alpha+beta,
// where alpha is an integer vector, and beta is a single double.
void rdirichlet2(std::vector<double> &x, std::vector<int> &alpha, double beta) {

  int n = x.size();
  double xSum = 0;
  for (int i=0; i<n; i++) {
    x[i] = rgamma1(alpha[i]+beta, 1.0);
    xSum += x[i];
  }
  double xSumInv = 1.0/xSum;
  for (int i=0; i<n; i++) {
    x[i] *= xSumInv;
    if (x[i] <= UNDERFLO) {
      x[i] = UNDERFLO;
    }
  }
}

//------------------------------------------------
// draw from Poisson(rate) distribution
int rpois1(double rate) {
  return R::rpois(rate);
}

//------------------------------------------------
// probability mass of Poisson(rate) distribution
double dpois1(int n, double rate, bool returnLog) {
  return R::dpois(n,rate,returnLog);
}

//------------------------------------------------
// draw from negative binomial distribution with mean lambda and variance gamma*lambda (gamma must be >1)
int rnbinom1(double lambda, double gamma) {
  return R::rnbinom(lambda/(gamma-1), 1/gamma);
}

//------------------------------------------------
// probability mass of negative binomial distribution with mean lambda and
// variance gamma*lambda (gamma must be >1)
double dnbinom1(int n, double lambda, double gamma, bool returnLog) {
  return R::dnbinom(n, lambda/(gamma-1), 1/gamma, returnLog);
}

//------------------------------------------------
// probability mass of negative binomial distribution with mean and variance 
double dnbinom_mu1(int n, double size, double mean, bool returnLog) {
  return R::dnbinom_mu(n, size, mean, returnLog);
}

//------------------------------------------------
// return closest value to a vector of target values
double closest(std::vector<double> const& vec, double value) {
  auto const it = std::lower_bound(vec.begin(), vec.end(), value);
  if (it == vec.end()) {
    return -1;
  }
  return *it;
}

//------------------------------------------------
// draw from binomial(N,p) distribution
int rbinom1(int N, double p) {
  if (p >= 1.0) {
    return N;
  } else if (p <= 0.0) {
    return 0;
  }
  return R::rbinom(N, p);
}

//------------------------------------------------
// draw from multinomial(N,p) distribution, where p sums to p_sum
void rmultinom1(int N, const vector<double> &p, double p_sum, vector<int> &ret) {
  int k = int(p.size());
  fill(ret.begin(), ret.end(), 0);
  for (int i = 0; i < (k-1); ++i) {
    ret[i] = rbinom1(N, p[i] / p_sum);
    N -= ret[i];
    if (N == 0) {
      break;
    }
    p_sum -= p[i];
  }
  ret[k-1] = N;
}
