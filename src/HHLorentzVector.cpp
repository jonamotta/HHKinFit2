#include "HHLorentzVector.h"
#include <cmath>

HHLorentzVector::HHLorentzVector(double x, double y, double z, double t)
  : TLorentzVector(x,y,z,t){
}

HHLorentzVector::HHLorentzVector()
  : TLorentzVector(){
}

void HHLorentzVector::SetEEtaPhiM(double E,double Eta,double Phi,double M){
  double p(sqrt(E*E-M*M));
  double sinth(sin(2*atan(exp(-Eta))));
  double pt(sinth*p);
  this->SetPtEtaPhiE(pt,Eta,Phi,E);
}
