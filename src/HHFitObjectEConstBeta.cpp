#include "HHFitObjectEConstBeta.h"
#include <cassert>
#include <iostream>


HHFitObjectEConstBeta::HHFitObjectEConstBeta(TLorentzVector const& initial4vector)
  :HHFitObjectE(initial4vector){
  //not yet implemented
  assert(0);
}

TLorentzVector 
HHFitObjectEConstBeta::constrainEtoMinv(Double_t m, TLorentzVector const& other4vector) const{
  TLorentzVector p(0,0,0,0);

  return(p);
}

void
HHFitObjectEConstBeta::print() const{
  std::cout << "---" << std::endl;
  std::cout <<  "energy component fit object with constant beta (=p/E):" << std::endl;
  this->printInitial4Vector();
  this->printFit4Vector();
  this->printCovMatrix();
}
