#include <iostream>
#include "HHTauTauEventGenerator.h"
#include "TF1.h"
#include "TVector2.h"
#include "TH1.h"
#include "TFile.h"
#include <cmath>


int main(int argc, char* argv[])
{
  TF1 PDF1("PDF1","2*x",0,2);
  TF1 PDF2("PDF2","2*x",0,2);
  HHTauTauEventGenerator testgenerator(PDF1,PDF2);
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
    TVector2 MET1(testgenerator.getMET());
    MET1.Print();
  }

  return(0);
}
