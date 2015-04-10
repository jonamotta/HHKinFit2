#include "HHLorentzVector.h"
#include <cmath>

HHLorentzVector::HHLorentzVector(double x, double y, double z, double t)
  : TLorentzVector(x,y,z,t){
}

HHLorentzVector::HHLorentzVector()
  : TLorentzVector(){
}

void HHLorentzVector::SetEEtaPhiM(double E, double eta, double phi, double m){
  double p = sqrt(E*E-m*m);
  double sinth = sin(2*atan(exp(-eta)));
  double pt = sinth*p;
  this->SetPtEtaPhiE(pt,eta,phi,E);
}

HHLorentzVector HHLorentzVector::operator+(HHLorentzVector const& rhs) const{
	return(HHLorentzVector(rhs.Px()+this->Px(),rhs.Py()+this->Py(),rhs.Pz()+this->Pz(),rhs.E()+this->E()));
}

HHLorentzVector HHLorentzVector::operator-(HHLorentzVector const& rhs) const{
	return(HHLorentzVector(rhs.Px()-this->Px(),rhs.Py()-this->Py(),rhs.Pz()-this->Pz(),rhs.E()-this->E()));
}
