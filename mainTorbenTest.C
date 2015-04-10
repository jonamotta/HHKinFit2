#include <iostream>
#include "HHTauTauEventGenerator.h"
#include "TF1.h"
#include "TVector2.h"


int main(int argc, char* argv[])
{
  TF1 *PDF = new TF1("PDF","sin(x)/x",0,10);
  HHTauTauEventGenerator testgenerator(PDF);
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
  testgenerator.generateEvent();
  std::cout << "testing the generator methods" << std::endl;
  std::cout << "Vector of the first tau in rest frame of the higgs:" << std::endl;
  HHLorentzVector tau1test2(testgenerator.getTau1());
  tau1test2.Print();
  std::cout << "Vector of the first tau boosted in lab frame:" << std::endl;
  HHLorentzVector tau1boostedtest2(testgenerator.getTau1boosted());
  tau1boostedtest2.Print();
  std::cout << "Vector of the second tau in rest frame of the higgs:" << std::endl;
  HHLorentzVector tau2test2(testgenerator.getTau2());
  tau2test2.Print();
  std::cout << "Vector of the second tau boosted in lab frame:" << std::endl;
  HHLorentzVector tau2boostedtest2(testgenerator.getTau2boosted());
  tau2boostedtest2.Print();
  std::cout << "teste MET" << std::endl;
  TVector2 MET2=testgenerator.getMET();
  MET2.Print();


  return(0);
}
