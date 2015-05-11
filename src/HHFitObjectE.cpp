#include "HHFitObjectE.h"
#include <iostream>
#include <iomanip>

#include "exceptions/HHEnergyRangeException.h"
#include "exceptions/HHLimitSettingException.h"


HHKinFit2::HHFitObjectE::HHFitObjectE(HHLorentzVector const& initial4vector)
  :HHFitObject(initial4vector),
   m_upperLimitE(pow(10,10)),
   m_lowerLimitE(0),
   m_initstart(-pow(10,10)),
   m_initprec(-pow(10,10)),
   m_initstep(-pow(10,10)),
   m_initdirection(-pow(10,10)){
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


void
HHKinFit2::HHFitObjectE::setFitLimitsE(double const lowerlimit, double const upperlimit){
  this->setLowerFitLimitE(lowerlimit);
  this->setUpperFitLimitE(upperlimit);
}

void
HHKinFit2::HHFitObjectE::setFitLimitsE(HHLorentzVector const& own4vectorMin, double const minv, HHLorentzVector const& other4vectorMin){
  this->setLowerFitLimitE(own4vectorMin);
  this->setUpperFitLimitE(minv,other4vectorMin);
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
  try{
    this->setUpperFitLimitE(constrainEtoMinv(minv,other4vectorMin).E());
  }
  catch(HHEnergyRangeException const& e){
    throw(HHLimitSettingException(e.what()));
  }
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
HHKinFit2::HHFitObjectE::setCovMatrix(double dE){
  TMatrixD cov(4,4);
  cov(3,3)=dE;
  m_covmatrix=cov;
}

void
HHKinFit2::HHFitObjectE::setInitPrecision(double prec){
  m_initprec=prec;
}

void
HHKinFit2::HHFitObjectE::setInitDirection(double daN){
  m_initdirection=daN;
}

void
HHKinFit2::HHFitObjectE::setInitStepWidth(double h){
  m_initstep=h;
}

void
HHKinFit2::HHFitObjectE::setInitStart(double start){
  m_initstart=start;
}

double
HHKinFit2::HHFitObjectE::getInitPrecision(){
  return(m_initprec);
}

double
HHKinFit2::HHFitObjectE::getInitDirection(){
  return(m_initdirection);
}

double
HHKinFit2::HHFitObjectE::getInitStepWidth(){
  return(m_initstep);
}

double
HHKinFit2::HHFitObjectE::getInitStart(){
  return(m_initstart);
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

