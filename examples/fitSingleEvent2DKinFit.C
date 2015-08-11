#include <iostream>

#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "HHKinFitMasterHeavyHiggs.h"


using HHKinFit2::HHKinFitMasterHeavyHiggs;

int main(int argc, char* argv[])
{

  TLorentzVector* bjet1   = new TLorentzVector(-46.36, -7.87, 122.31,131.53);
  TLorentzVector* bjet2   = new TLorentzVector(116.81, 15.80, -59.72,132.79);
  TLorentzVector* tauvis1 = new TLorentzVector( -8.26, -21.85, 15.30, 27.93);
  TLorentzVector* tauvis2 = new TLorentzVector(-45.02, 54.19, 35.25, 78.78);
  TLorentzVector* met     = new TLorentzVector(-74.62, 7.71,0,0);
  TMatrixD        met_cov(4,4);
  met_cov[0][0]= 575.57;
  met_cov[0][1]= -46.76;
  met_cov[1][0]= met_cov[0][1];
  met_cov[1][1]= 744.59;

  HHKinFitMasterHeavyHiggs heavyhiggsfit(bjet1, bjet2, tauvis1, tauvis2, met, met_cov);
  heavyhiggsfit.addHypo(125,125);
  heavyhiggsfit.addHypo(125,150);
  heavyhiggsfit.addHypo(125,175);
  heavyhiggsfit.addHypo(135,125);
  heavyhiggsfit.doFit();

  return (0);
}
