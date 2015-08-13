#include <iostream>

#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "HHKinFitMasterHeavyHiggs.h"


using HHKinFit2::HHKinFitMasterHeavyHiggs;

int main(int argc, char* argv[])
{
  TLorentzVector* bjet1   = new TLorentzVector(-65.31,50.89,-97.31,128.13);
  TLorentzVector* bjet2   = new TLorentzVector(280.11,225.64,-517.33,633.32);
  TLorentzVector* tauvis1 = new TLorentzVector(-64.86,-54.66,-178.57,197.69);
  TLorentzVector* tauvis2 = new TLorentzVector(-50.37,-132.73,-444.31,466.44);
  TLorentzVector* met     = new TLorentzVector(-71.73,-122.64,0,0);
  TMatrixD        met_cov(4,4);
  met_cov[0][0]= 528.94;
  met_cov[0][1]= 90.30;
  met_cov[1][0]= met_cov[0][1];
  met_cov[1][1]= 727.47;

  HHKinFitMasterHeavyHiggs heavyhiggsfit(bjet1, bjet2, tauvis1, tauvis2, met, met_cov);
  heavyhiggsfit.addHypo(125,125);
  heavyhiggsfit.doFit();

  return (0);
}
