#include "HHFitObjectMET.h"
#include <iostream>
#include <iomanip>


HHFitObjectMET::HHFitObjectMET(TVector2 const& v,TVector2 const& vfit)
  :HHFitObject(TLorentzVector(v.X(),v.Y(),0,sqrt(v.X()*v.X()+v.Y()*v.Y()))){
  this->setFit4Vector(TLorentzVector(vfit.X(),vfit.Y(),0,sqrt(vfit.X()*vfit.X()+vfit.Y()*vfit.Y())));
}


void
HHFitObjectMET::setCovMatrix(Double_t xx, Double_t yy, Double_t xy){
  TMatrixD cov(4,4);
  cov(0,0)=xx;
  cov(1,1)=yy;
  cov(0,1)=xy;
  cov(1,0)=xy;
  m_covmatrix=cov;
}


void
HHFitObjectMET::print() const{
  std::cout << "---" << std::endl;
  std::cout << "MET fit object:" << std::endl;
  this->printInitial4Vector();
  this->printFit4Vector();
  this->printCovMatrix();
}


void
HHFitObjectMET::printInitial4Vector() const{
  std::cout <<  "initial vector (px,py,pt)"
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getInitial4Vector().Px()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getInitial4Vector().Py()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getInitial4Vector().Pt()
            << std::endl;
}

void
HHFitObjectMET::printFit4Vector() const{
  std::cout <<  "  final vector (px,py,pt)"
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getFit4Vector().Px()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getFit4Vector().Py()
            << std::setw(10) << std::fixed<< std::setprecision(1) << this->getFit4Vector().Pt()
            << std::endl;
}

void
HHFitObjectMET::printCovMatrix() const{
  std::cout <<  "covariance matrix:" << std::endl;
  std::cout   << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(0,0)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(0,1) << std::endl;
  std::cout   << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(1,0)
              << std::setw(10) << std::fixed<< std::setprecision(2) << this->getCovMatrix()(1,1) << std::endl;
}
