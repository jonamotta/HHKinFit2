#ifdef HHKINFIT2
#include "HHFitConstraint4VectorBJet.h"
#include "exceptions/HHCovarianceMatrixException.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint4VectorBJet.h"
#include "HHKinFit2/HHKinFit2/interface/exceptions/HHCovarianceMatrixException.h"
#endif

#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"

#include <Math/SpecFuncMathCore.h>
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFuncMathCore.h>

#include <iostream>
#include <sstream>
#include <cmath>

HHKinFit2::HHFitConstraint4VectorBJet::HHFitConstraint4VectorBJet(HHFitObject* object,
								  double alpha, 
								  double n, 
								  double sigma, 
								  double mean,
								  double beta, 
								  double normalization)
  
  : HHFitConstraint(object)
{
  m_alpha = alpha;
  m_n = n;
  m_sigma = sigma;
  m_mean = mean;
  m_beta = beta;
  m_normalization = normalization;
}

HHKinFit2::HHFitConstraint4VectorBJet::HHFitConstraint4VectorBJet(HHFitObject* object)
  : HHFitConstraint(object)
{
  HHLorentzVector initVec = m_fitobject->getInitial4Vector();
  
  double initVecEt = initVec.Et();
  double initVecEta = fabs(initVec.Eta());  

  if(initVecEta < 1.2)
  {
    if(initVecEt < 25)
    {
      m_alpha = 2.48532;
      m_n = 99.912;
      m_sigma = 0.248993;
      m_mean = 1.35457;
      m_beta = -0.525116;
      m_normalization = 1.16973;
    }
    else if(initVecEt > 25 &&  initVecEt < 30)
    {
      m_alpha = 1.94202;
      m_n = 99.918;
      m_sigma = 0.20096;
      m_mean = 1.24588;
      m_beta = 0.50405;
      m_normalization = 1.41556;
    }
    else if(initVecEt > 30 &&  initVecEt < 40)
    {
      m_alpha = 2.06421;
      m_n = 99.9267;
      m_sigma = 0.170023;
      m_mean = 1.14678;
      m_beta = 0.466032;
      m_normalization = 1.60844;
    }
    else if(initVecEt > 40 &&  initVecEt < 50)
    {
      m_alpha = 1.77141;
      m_n = 99.9939;
      m_sigma = 0.16554;
      m_mean = 1.12466;
      m_beta = -0.591766;
      m_normalization = 1.84852;
    }
    else if(initVecEt > 50 &&  initVecEt < 60)
    {
      m_alpha = 1.67154;
      m_n = 99.9779;
      m_sigma = 0.161855;
      m_mean = 1.1019;
      m_beta = 0.716851;
      m_normalization = 2.03644;
    }
    else if(initVecEt > 60 &&  initVecEt < 80)
    {
      m_alpha = 1.43124;
      m_n = 47.6301;
      m_sigma = 0.149199;
      m_mean = 1.08066;
      m_beta = 0.839028;
      m_normalization = 2.29145;
    }
    else if(initVecEt > 80 &&  initVecEt < 120)
    {
      m_alpha = 1.22277;
      m_n = 14.1823;
      m_sigma = 0.118743;
      m_mean = 1.02656;
      m_beta = 1.04878;
      m_normalization = 2.93557;
    }
    else if(initVecEt > 120 &&  initVecEt < 500)
    {
      m_alpha = 0.70576;
      m_n = 6.70179;
      m_sigma = 0.107425;
      m_mean = 0.993853;
      m_beta = 1.21589;
      m_normalization = 2.75108;
    }
  }
  else
  {
    if(initVecEt < 25)
    {
      m_alpha = 2.04106;
      m_n = 97.4338;
      m_sigma = 0.246799;
      m_mean = 1.31924;
      m_beta = 0.524127;
      m_normalization = 1.17627;
    }
    else if(initVecEt > 25 &&  initVecEt < 30)
    {
      m_alpha = 2.52252;
      m_n = 99.5452;
      m_sigma = 0.226336;
      m_mean = 1.25191;
      m_beta = 0.533287;
      m_normalization = 1.29661;
    }
    else if(initVecEt > 30 &&  initVecEt < 40)
    {
      m_alpha = 2.12896;
      m_n = 99.9998;
      m_sigma = 0.198;
      m_mean = 1.1731;
      m_beta = 0.589841;
      m_normalization = 1.55054;
    }
    else if(initVecEt > 40 &&  initVecEt < 50)
    {
      m_alpha = 1.8613;
      m_n = 99.9795;
      m_sigma = 0.191035;
      m_mean = 1.13761;
      m_beta = 0.680676;
      m_normalization = 1.69989;
    }
    else if(initVecEt > 50 &&  initVecEt < 60)
    {
      m_alpha = 1.6989;
      m_n = 99.7556;
      m_sigma = 0.173905;
      m_mean = 1.09708;
      m_beta = 0.722451;
      m_normalization = 1.90215;
    }
    else if(initVecEt > 60 &&  initVecEt < 80)
    {
      m_alpha = 1.47361;
      m_n = 29.041;
      m_sigma = 0.151756;
      m_mean = 1.07486;
      m_beta = 0.831271;
      m_normalization = 2.24056;
    }
    else if(initVecEt > 80 &&  initVecEt < 120)
    {
      m_alpha = 1.31405;
      m_n = 6.86929;
      m_sigma = 0.124757;
      m_mean = 1.02943;
      m_beta = 1.21758;
      m_normalization = 2.83871;
    }
    else if(initVecEt > 120 &&  initVecEt < 500)
    {
      m_alpha = 0.499206;
      m_n = 5.10681;
      m_sigma = 0.0679968;
      m_mean = 0.988872;
      m_beta = 0.819638;
      m_normalization = 3.26882;
    }
  }
}

double
HHKinFit2::HHFitConstraint4VectorBJet::getChi2() const{
  HHLorentzVector fitVec = m_fitobject->getFit4Vector();
  HHLorentzVector initVec = m_fitobject->getInitial4Vector();
  
  double res = fitVec.E()/initVec.E();
  double cdf = crystBallLikeCDF(res, m_alpha, m_n, m_sigma, m_mean, m_beta, 
					   m_normalization);
  double chi2 = 2*pow(TMath::ErfInverse(2.0*cdf - 1),2);
  
  return chi2;
}

void
HHKinFit2::HHFitConstraint4VectorBJet::printChi2() const{
  HHLorentzVector fitVec = m_fitobject->getFit4Vector();
  HHLorentzVector initVec = m_fitobject->getInitial4Vector();
  
  double res = fitVec.E()/initVec.E();
  std::cout << "Residual is " << res << std::endl;
  double cdf = crystBallLikeCDF(res, m_alpha, m_n, m_sigma, m_mean, m_beta, 
			    m_normalization);
  std::cout << "cdf is " << cdf << std::endl;
  double chi2 = 2*pow(TMath::ErfInverse(2.0*cdf - 1),2);
  std::cout << "Chi2 is " << chi2 << std::endl;

  fitVec.Print();
  initVec.Print();
  std::cout << "Chi2 is: " << chi2 << std::endl;
}

double
HHKinFit2::HHFitConstraint4VectorBJet::getLikelihood() const{
  double likelihood = 1.0;
  std::cout << "WARNING! Likelihood not implemented for BJets!!!" << std::endl;
  return(likelihood);
}

double HHKinFit2::crystBallLikeCDF(double x, double alpha, 
				   double n, double sigma, 
				   double mean,
				   double beta, 
				   double normalization)
{
  if (sigma == 0)   return 0;
  bool useLog = (n == 1.0); 

  double z = (x-mean)/sigma;
  //if (alpha < 0 ) z = -z;
 
  double abs_alpha = std::abs(alpha);
  double abs_beta = std::abs(beta);

  double intgaus = 0.;
  double intpow  = 0.;

  if (z <= -abs_alpha)
  {
    double A = std::pow(n/abs_alpha,n) * std::exp(-0.5 * alpha*alpha);
    double B = n/abs_alpha - abs_alpha;

    if (!useLog) {
      intpow  = (A/(n-1.0) * std::pow(B-z,-n+1))*sigma ;
    }
    else {
      // for n=1 the primitive of 1/x is log(x)
      intpow = -A * std::log( n / abs_alpha ) + A * std::log( B -z );
    }
  }
  else if (z  > - abs_alpha && z < abs_beta)
  {
    double A = std::pow(n/abs_alpha,n) * std::exp(-0.5 * alpha*alpha);
    double B = n/abs_alpha - abs_alpha;

    if (!useLog) {
      intpow  = (A/(n-1.0) * std::pow(B+abs_alpha,-n+1))*sigma ;
    }
    else {
      // for n=1 the primitive of 1/x is log(x)
      intpow = -A * std::log( n / abs_alpha ) + A * std::log( B +abs_alpha );
    }
    
    double xAtAlpha = -(abs_alpha * sigma) + mean;

    intgaus = (ROOT::Math::normal_cdf(x, sigma, mean) - ROOT::Math::normal_cdf(xAtAlpha, sigma, mean))*(std::sqrt( 2.*M_PI)*sigma);
  }
  else
  {
    
    double A = std::pow(n/abs_alpha,n) * std::exp(-0.5 * alpha*alpha);
    double B = n/abs_alpha - abs_alpha;

    if (!useLog) {
      intpow  = (A/(n-1.0) * std::pow(B+abs_alpha,-n+1))*sigma ;
    }
    else {
      // for n=1 the primitive of 1/x is log(x)
      intpow = -A * std::log( n / abs_alpha ) + A * std::log( B +abs_alpha );
    }
    
    
    double xAtAlpha = -(abs_alpha * sigma) + mean;
    double xAtBeta = (abs_beta * sigma) + mean;

    intgaus = (ROOT::Math::normal_cdf(xAtBeta, sigma, mean) - ROOT::Math::normal_cdf(xAtAlpha, sigma, mean))*(std::sqrt( 2.*M_PI)*sigma);
    

    double C = std::pow(n/abs_beta,n) * std::exp(-0.5 * beta*beta);
    double D = n/abs_beta - abs_beta;

    if (!useLog) {
      intpow  += (-C/(n-1.0) * std::pow(D+z,-n+1) + C/(n-1.0) * std::pow(D+abs_beta,-n+1))*sigma ;
    }
    else {
      // for n=1 the primitive of 1/x is log(x)
      intpow += -C * std::log( n / abs_beta ) + C * std::log( D +z );
    }
    
  }

  return normalization*(intgaus + intpow);
}
