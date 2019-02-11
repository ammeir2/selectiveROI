#include <Rcpp.h>
using namespace Rcpp;

const double pi  = 3.141592653589793238463;

///// SLICE SAMPLER FUNCTIONS

void printVec(NumericVector x) {
  for(int i = 0 ; i < x.length() ; i++) {
    Rcpp::Rcout<<x[i]<<" " ;
  }
  Rcpp::Rcout<<"\n" ;
  return ;
}

// target = M %*% x
void innerMatVec(NumericVector target, NumericMatrix M, NumericVector x) {
  for(int i = 0 ; i < M.nrow() ; i++) {
    target[i] = 0.0 ;
    for(int j = 0; j < M.ncol() ; j++) {
      target[i] += M(i, j) * x[j] ;
    }
  }
}

// all(large >= small)
bool allLarge(NumericVector large, NumericVector small) {
  for(int i = 0 ; i < large.length() ; i++) {
    if(large[i] < small[i]) {
      return false ;
    }
  }

  return true ;
}

// target = X * x + Y * y
void addVecTimesConst(NumericVector target, NumericVector X, double x, NumericVector Y, double y) {
  for(int i = 0 ; i < X.length() ; i++) {
    target[i] = X[i] * x + Y[i] * y ;
  }
}

// [[Rcpp::export]]
void sliceSamplerRcpp(NumericMatrix sampMat,
                      NumericVector samp,
                      NumericMatrix chol,
                      NumericVector lth, NumericVector uth) {
  int p = samp.length() ;
  NumericVector v(p) ;
  NumericVector z(p) ;
  double likThreshold ;
  double u ;
  double theta, thMin, thMax ;
  bool keep ;
  NumericVector proposal(p) ;

  for(int m = 0; m < sampMat.nrow() ; m++) {
    z = rnorm(p) ;
    innerMatVec(v, chol, z) ;
    u = runif(1)[0] ;
    likThreshold = std::log(u) ;
    theta = R::runif(0.0, 2.0 * pi) ;
    thMin = theta - 2.0 * pi ;
    thMax = theta ;
    keep = false ;
    int attempts = 0;
    while(!keep) {
      attempts++ ;
      addVecTimesConst(proposal, samp, std::cos(theta), v, std::sin(theta)) ;
      keep = allLarge(proposal, uth) | allLarge(lth, proposal) ;
      if(keep) break ;
      if(theta < 0) {
        thMin = theta ;
      } else {
        thMax = theta ;
      }
      theta = runif(1, thMin, thMax)[0] ;
      if(attempts > 5000) {
        stop("Something went wrong, cant find a sample value that satisfies constraints, please check input values.") ;
      }
    }
    std::copy(proposal.begin(), proposal.end(), samp.begin()) ;
    std::copy(samp.begin(), samp.end(), sampMat.row(m).begin()) ;
    // Rcpp::Rcout<<m<<"\n" ;
    // printVec(proposal) ;
  }

  return ;
}


//////// GIBBS SAMPELR FUNCTIONS -------
double sampleExtreme(double mu, double sd, double lower, double upper) {
  double sign = 1 ;
  double threshold ;
  double proposal ;
  double alpha ;
  double phi ;

  if(std::isinf(lower)) {
    sign = -1 ;
    mu *= sign ;
    threshold = upper * sign ;
  } else {
    threshold = lower ;
  }

  // rescaling
  threshold = (threshold - mu) / sd ;
  alpha = (threshold + std::sqrt(std::pow(threshold, 2) + 4)) / 2 ;

  bool reject = true ;
  int iter = 0;
  while(reject & (iter++ < 10000)) {
    proposal = threshold + R::rexp(1 / alpha) ;
    phi = std::exp(-std::pow(proposal - alpha, 2) / 2) ;
    if(runif(1)[0] < phi) {
      reject = false ;
    }
  }
  //Rcpp::Rcout<<iter<<"\n" ;

  proposal = proposal * sd + mu ;
  return proposal * sign;
}

// [[Rcpp::export]]
double sampleUnivTruncNorm(double mu, double sd, double lower, double upper) {
  double u = runif(1)[0] ;
  double phiB = R::pnorm5(upper, mu, sd, 1, 0) ;
  double phiA = R::pnorm5(lower, mu, sd, 1, 0) ;
  double quantile = u * phiB + phiA * (1 - u) ;
  double sample = R::qnorm(quantile, mu, sd, 1, 0) ;
  //Rcpp::Rcout<<"samp: "<<sample<<"\n" ;

  // Sometimes the sampler doesn't work well, in those cases we use a rejection
  // sampler. They do something similar in the truncnorm package.
  while(sample < lower | sample > upper | std::isinf(std::abs(sample))) {
    // throw std::range_error("Inadmissible value") ;
    sample = sampleExtreme(mu, sd, lower, upper) ;
  }

  return sample ;
}

double computeConditionalMean(NumericVector mu,
                              NumericVector samp,
                              const NumericMatrix precision,
                              int index) {
  double result = 0 ;
  double cor = 0;

  for(int j = 0; j < mu.length() ; j ++) {
    if(j != index) {
      result += precision(index, j) * (samp[j] - mu[j]) ;
      cor += precision(index, j) ;
    }
  }

  result = result / precision(index, index) ;
  //Rcpp::Rcout<<cor<<" "<<result<<" "<< samp[index]<<" ";
  result = mu[index] - result ;

  return result ;
}

// [[Rcpp::export]]
NumericVector sampleTruncNorm(NumericVector sample,
                              NumericVector lower, NumericVector upper,
                              NumericVector mean, NumericMatrix precision,
                              int cycles) {
  double condMean ;
  double condSD ;
  int i, j ;
  int EXTREME_THRESHOLD = 4 ;
  NumericVector samp = clone(sample) ;

  for(i = 0; i < cycles ; i++) {
    for(j = 0; j < samp.length() ; j ++) {
      condMean = computeConditionalMean(mean, samp, precision, j) ;
      //Rcpp::Rcout<<mean[j]<<" "<<condMean<<"\n" ;
      condSD = std::sqrt(1 / precision(j, j)) ;
      if((std::isinf(std::abs(lower[j])) & ((condMean - upper[j]) / condSD) > EXTREME_THRESHOLD) |
         (std::isinf(std::abs(upper[j])) & ((lower[j] - condMean) / condSD) > EXTREME_THRESHOLD)) {
        samp[j] = sampleExtreme(condMean, condSD, lower[j], upper[j]) ;
      } else {
        samp[j] = sampleUnivTruncNorm(condMean, condSD, lower[j], upper[j]) ;
      }
    }
    //Rcpp::Rcout<<samp<<"\n" ;
  }

  return samp ;
}

