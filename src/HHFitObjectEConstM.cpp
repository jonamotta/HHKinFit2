#include "HHFitObjectEConstM.h"
#include <cassert>
#include <iostream>
#include <sstream>
#include <math.h>
#include "exceptions/HHEnergyRangeException.h"
#include "exceptions/HHEnergyConstraintException.h"


HHFitObjectEConstM::HHFitObjectEConstM(HHLorentzVector const& initial4vector)
  :HHFitObjectE(initial4vector){

}


HHLorentzVector
HHFitObjectEConstM::constrainEtoMinv(Double_t m, HHLorentzVector const& pset) const{
  HHLorentzVector pmod = getInitial4Vector();
  HHLorentzVector combined = pset + pmod;

  Double_t Mc = m;
  Double_t M1c = pset.M();
  Double_t M2c = pmod.M();
//  Double_t M = combined.M();

  Double_t C = 0.5*(pow(Mc,2)-pow(M1c,2)-pow(M2c,2));

//  int loopCount = 0;
//  while(fabs(M-Mc) > 0.000001){
//    loopCount++;
    Double_t P1x = pset.Px();
    Double_t P1y = pset.Py();
    Double_t P1z = pset.Pz();
    Double_t P1  = pset.P();
    Double_t E1  = pset.E();

    Double_t P2x = pmod.Px();
    Double_t P2y = pmod.Py();
    Double_t P2z = pmod.Pz();
    Double_t P2  = pmod.P();

    Double_t cosa = (P1x*P2x+P1y*P2y+P1z*P2z)/(P1*P2);
    Double_t E2new = -1;

    if(cosa==0){
      E2new=C/E1;
    }
    else{
      Double_t cp=C/(cosa*P1);
      Double_t dp=E1/(cosa*P1);
      Double_t a=pow(dp,2)-1;
      Double_t b=-2*dp*cp;
      Double_t c=pow(cp,2)+pow(M2c,2);

      if (cosa>0) E2new = (1./(2*a))*(-b+sqrt(pow(b,2)-4*a*c));
      if (cosa<0) E2new = (1./(2*a))*(-b-sqrt(pow(b,2)-4*a*c));
    }

    if (isnan(E2new)||isinf(E2new)){
      std::stringstream msg;
      msg << "problem in constraining energy to inv. mass: minv="<<Mc<< " E(set)="<<E2new;
      throw(HHEnergyConstraintException(msg.str()));
    }

    return(this->changeE(E2new));

//    Double_t P2new = sqrt(pow(E2new,2)-pow(M2c,2));
//    Double_t Pt2new = P2new * sin(2.*atan(exp(-pmod.Eta())));
//    pmod.SetPtEtaPhiE(Pt2new,pmod.Eta(),pmod.Phi(),E2new);
//    combined = pmod + pset;
//    M = combined.M();
//  }
//  std::cout << "needed " << loopCount << " iterations to constrain E on invariant mass." << std::endl;
//  return(pmod);
}

HHLorentzVector
HHFitObjectEConstM::changeE(Double_t E) const{
  HHLorentzVector temp = this->getFit4Vector();

  if(E<temp.M()){
    std::stringstream msg;
    msg << "target energy is smaller than the particle mass: "<<"E(set)="<<E<<" "<<"m="<<temp.M();
    throw(HHEnergyRangeException(msg.str()));
  }
  else if((E<this->getLowerFitLimitE())||(E>this->getUpperFitLimitE())){
    std::stringstream msg;
    msg << "target energy is out of limits: "<<"E(set)="<<E<<" "<<"E(limits)=["<<this->getLowerFitLimitE()<<","<< this->getUpperFitLimitE() << "]";
    throw(HHEnergyRangeException(msg.str()));
  }

  Double_t pnew = sqrt(pow(E,2)-pow(temp.M(),2));
  Double_t ptnew = pnew * sin(2.*atan(exp(-temp.Eta())));
  temp.SetPtEtaPhiE(ptnew,temp.Eta(),temp.Phi(),E);

  return(temp);
}

HHLorentzVector
HHFitObjectEConstM::scaleE(Double_t scale) const{
  return(this->changeE(scale*this->getFit4Vector().E()));
}

void
HHFitObjectEConstM::print() const{
  std::cout << "---" << std::endl;
  std::cout <<  "energy component fit object with constant mass:" << std::endl;
  this->printInitial4Vector();
  this->printFit4Vector();
  this->printLimits();
//  this->printCovMatrix();
}
