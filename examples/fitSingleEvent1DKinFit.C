#include <iostream>
#include <iomanip>

#include "TLorentzVector.h"
#include "TVector2.h"
#include "TMatrixD.h"
#include "HHKinFitMasterSingleHiggs.h"

#include "TGraph.h" 
#include "TFile.h" 

using HHKinFit2::HHKinFitMasterSingleHiggs;

int main(int argc, char* argv[])
{
  TLorentzVector tauvis1(-15.902388,37.502301,43.366285,59.503192);
  TLorentzVector tauvis2(16.824116,-21.066462,7.050932,27.867069);
  TVector2 met(8.283318,11.756746);
  TMatrixD met_cov(4,4);
  met_cov[0][0]= 359.15655;
  met_cov[0][1]= 0;
  met_cov[1][0]= met_cov[0][1];
  met_cov[1][1]= 359.15655;

  HHKinFitMasterSingleHiggs singlehiggsfit(tauvis1,tauvis2,met,met_cov);

  for (int i=50; i<=250; i++)
    singlehiggsfit.addHypo(i);
 
  singlehiggsfit.fit();

  TFile f("out1D.root","RECREATE");
  TGraph gr(0);
  int points=0;
  for (int i=50; i<=250; i++){
    gr.SetPoint(points,i,singlehiggsfit.getChi2(i));
    std::cout << std::setprecision(5) << singlehiggsfit.getChi2(i) << std::endl;
    points++;
  }

  gr.Write();
  f.Close();

  return (0);
}
