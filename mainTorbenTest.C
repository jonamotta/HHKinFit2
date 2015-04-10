#include <iostream>
#include "HHTauTauEventGenerator.h"


int main(int argc, char* argv[])
{
  HHTauTauEventGenerator testgenerator;
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



  return(0);
}
