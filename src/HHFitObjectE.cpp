#include "HHFitObjectE.h"
#include <iostream>
#include <iomanip>


HHKinFit2::HHFitObjectE::HHFitObjectE(HHLorentzVector const& initial4vector)
  :HHFitObject(initial4vector),
   m_upperLimitE(9999999),
   m_lowerLimitE(0){

}

HHKinFit2::HHLorentzVector
HHKinFit2::HHFitObjectE::scaleE(double scale) const{
  return(this->changeE(scale*this->getFit4Vector().E()));
}


void
HHKinFit2::HHFitObjectE::changeEandSave(double E){
  this->setFit4Vector(changeE(E));
}

void
HHKinFit2::HHFitObjectE::scaleEandSave(double E){
  this->setFit4Vector(scaleE(E));
}

void
HHKinFit2::HHFitObjectE::constrainEtoMinvandSave(double m, HHLorentzVector const& other4vector){
  this->setFit4Vector(constrainEtoMinv(m, other4vector));
}

double
HHKinFit2::HHFitObjectE::getE() const{
  return(m_fit4vector.E());
}

double
HHKinFit2::HHFitObjectE::getUpperFitLimitE() const{
  return(m_upperLimitE);
}

double
HHKinFit2::HHFitObjectE::getLowerFitLimitE() const{
  return(m_lowerLimitE);
}

void
HHKinFit2::HHFitObjectE::setUpperFitLimitE(double upperlimit){
  m_upperLimitE = upperlimit;
}

void
HHKinFit2::HHFitObjectE::setUpperFitLimitE(double minv, HHLorentzVector const& other4vectorMin){
  this->setUpperFitLimitE(constrainEtoMinv(minv,other4vectorMin).E());
}

void
HHKinFit2::HHFitObjectE::setLowerFitLimitE(double lowerlimit){
  m_lowerLimitE = lowerlimit;
}

void
HHKinFit2::HHFitObjectE::setLowerFitLimitE(HHLorentzVector const& lowerlimit){
  m_lowerLimitE = lowerlimit.E();
}

void
HHKinFit2::HHFitObjectE::print() const{
  std::cout << "---" << std::endl;
  std::cout << "energy component fit object:" << std::endl;
  this->printInitial4Vector();
  this->printFit4Vector();
  this->printLimits();
  this->printCovMatrix();
}


void
HHKinFit2::HHFitObjectE::printLimits() const{
  std::cout <<  "limits: "
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getLowerFitLimitE()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getUpperFitLimitE()
            << std::endl;
}

