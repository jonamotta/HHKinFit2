#include <iostream>
#include "HHTauTauEventGenerator.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TArrayD.h"
#include "TH1.h"
#include "TFile.h"
#include <cmath>

using namespace HHKinFit2;

int main(int argc, char* argv[])
{
  TF1 PDF1("PDF1","2*x",0,2);
  TF1 PDF2("PDF2","2*x",0,2);
  TMatrixD covarmatrix(2,2);
  covarmatrix[0][0]=130;
  covarmatrix[0][1]=0;
  covarmatrix[1][0]=0;
  covarmatrix[1][1]=130;
  HHTauTauEventGenerator testgenerator(PDF1,PDF2,covarmatrix);
  for(unsigned int i=0; i<10; i++){
    testgenerator.generateEvent();
    std::cout << "testing the generator methods" << std::endl;
    std::cout << "Vector of the first tau in rest frame of the higgs:" << std::endl;
    HHLorentzVector tau1test(testgenerator.getTau1());
    tau1test.Print();
    std::cout << "Vector of the first tau boosted in lab frame:" << std::endl;
    HHLorentzVector tau1boostedtest(testgenerator.getTau1boosted());
    tau1boostedtest.Print();
    std::cout << "Vector of the second tau in rest frame of the higgs:" << std::endl;
    HHLorentzVector tau2test(testgenerator.getTau2());
    tau2test.Print();
    std::cout << "Vector of the second tau boosted in lab frame:" << std::endl;
    HHLorentzVector tau2boostedtest(testgenerator.getTau2boosted());
    tau2boostedtest.Print();
    std::cout << "teste MET" << std::endl;
    TVectorD MET1(testgenerator.getMET());
    MET1.Print();
    TVectorD testergebnis(testgenerator.getMETwithsigma());
    testgenerator.getMET().Print();
    testgenerator.PrintCovarmatrix();
    testergebnis.Print();
  }
  std::cout << "testing Matrix:" << std::endl;
  std::cout << "Covarmatrix:" << std::endl;
  testgenerator.PrintCovarmatrix();
  std::cout << "L-Matrix:" << std::endl;
  testgenerator.PrintLmatrix();

  return(0);
}
