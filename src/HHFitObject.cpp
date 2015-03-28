#include "HHFitObject.h"

#include <iostream>
#include <iomanip>

HHFitObject::HHFitObject()
  :m_fit4vector(TLorentzVector()),
   m_initial4vector(TLorentzVector()),
   m_covmatrix(TMatrixD(4,4)){

}

HHFitObject::HHFitObject(TLorentzVector initial4vector)
  :m_fit4vector(initial4vector),
   m_initial4vector(initial4vector),
   m_covmatrix(TMatrixD(4,4)){

}

TLorentzVector  
HHFitObject::getInitial4Vector(){
  return(m_initial4vector);
}

TLorentzVector
HHFitObject::getFit4Vector(){
  return(m_fit4vector);
}

TMatrixD
HHFitObject::getCovMatrix(){
  return(m_covmatrix);
}

void
HHFitObject::setFit4Vector(TLorentzVector vec){
  this->m_fit4vector=vec;
}

void
HHFitObject::setCovMatrix(TMatrixD covmatrix){
  this->m_covmatrix=covmatrix;
}

void
HHFitObject::reset(){
  this->m_fit4vector=this->m_initial4vector;
}

void
HHFitObject::print(){
  std::cout << "---" << std::endl;
  std::cout <<  "general fit object:" << std::endl;
  this->printInitial4Vector();
  this->printFit4Vector();
  this->printCovMatrix();
}

void
HHFitObject::printInitial4Vector(){
  std::cout <<  "initial vector (px,py,pz,E,m)"
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getInitial4Vector().Px()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getInitial4Vector().Py()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getInitial4Vector().Pz()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getInitial4Vector().E()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getInitial4Vector().M()
            << std::endl;
}

void
HHFitObject::printFit4Vector(){
  std::cout <<  "  final vector (px,py,pz,E,m)"
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getFit4Vector().Px()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getFit4Vector().Py()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getFit4Vector().Pz()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getFit4Vector().E()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getFit4Vector().M()
            << std::endl;
}

void
HHFitObject::printCovMatrix(){
  std::cout <<  "covariance matrix:" << std::endl;
  std::cout   << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(0,0)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(0,1)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(0,2)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(0,3) << std::endl;
  std::cout   << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(1,0)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(1,1)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(1,2)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(1,3) << std::endl;
  std::cout   << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(2,0)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(2,1)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(2,2)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(2,3) << std::endl;
  std::cout   << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(3,0)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(3,1)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(3,2)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(3,3) << std::endl;
}
