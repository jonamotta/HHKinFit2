#include <iostream>

#include "TLorentzVector.h"
#include "TVector2.h"
#include "TMatrixD.h"
#include "HHKinFitMasterSingleHiggs.h"

#include "TGraph.h" 
#include "TFile.h" 

using HHKinFit2::HHKinFitMasterSingleHiggs;

int main(int argc, char* argv[])
{
  TLorentzVector tauvis1(27.0386,-3.37592,10.0049,29.04);
  TLorentzVector tauvis2(-19.7041,7.58246,-1.47997,21.1647);
  TVector2 met(-20.2176, 7.48374);
  TMatrixD met_cov(4,4);
  met_cov[0][0]= 191.604;
  met_cov[0][1]= -3.8281;
  met_cov[1][0]= met_cov[0][1];
  met_cov[1][1]= 156.727;

  HHKinFitMasterSingleHiggs singlehiggsfit(tauvis1,tauvis2,met,met_cov);

  for (int i=50; i<=250; i++)
    singlehiggsfit.addHypo(i);
 
  singlehiggsfit.fit();

  TFile f("out1D.root","RECREATE");
  TGraph gr(0);
  int points=0;
  for (int i=50; i<=250; i++){
    gr.SetPoint(points,i,singlehiggsfit.getChi2(i));
    std::cout << singlehiggsfit.getChi2(i) << std::endl;
    points++;
  }

  gr.Write();
  f.Close();

  return (0);
}
