#include "HHTauTauEventGenerator.h"
#include <cmath>
#include <iostream>
#include "TVector3.h"
HHTauTauEventGenerator::HHTauTauEventGenerator(TF1* a){
 m_PDF=a;
 m_visfrac1=1;
 m_visfrac2=1;
}
void HHTauTauEventGenerator::generateEvent() {
// Generate lorentzvectors for the two tau in the rest frame of the higgs
  
//Tau 1
  double pi(3.14159265359);
  double mhiggs(125.7);
  double mtau1(1.77682);
  double mtau2(1.77682);
  double cthtau1(m_randomnumber.Uniform(-1,1));
  double etatau1(log(tan(acos(cthtau1)/2)));
  double phitau1(m_randomnumber.Uniform(0,2*pi));
  double etau1((mhiggs*mhiggs-mtau2*mtau2 +mtau1*mtau1)/(2*mhiggs));
  m_tau1.SetEEtaPhiM(etau1, etatau1, phitau1, mtau1);
  //Tau 2
  double etau2((mhiggs*mhiggs-mtau1*mtau1 +mtau2*mtau2)/(2*mhiggs));
  double etatau2(-etatau1);
  double phitau2(phitau1+pi);
  m_tau2.SetEEtaPhiM(etau2, etatau2, phitau2, mtau2);

  //generate lorentzvectors for isr-jet
  double ptisr(fabs(m_randomnumber.Gaus(0,30))); // absolut Value of Pt
  double phiisr(m_randomnumber.Uniform(0,2*pi));
  double etaisr(m_randomnumber.Uniform(-5,5));
  double misr(0); //mass
  m_isr.SetPtEtaPhiM(ptisr,etaisr,phiisr,misr);
  
  //generate lorentzvector for higgs with lorentzvector from isr-jet
  double cthhiggs(m_randomnumber.Uniform(-1,1));
  double etahiggs(log(tan(acos(cthhiggs)/2)));
  double pthiggs(ptisr);
  double phihiggs(phiisr+pi);
  m_higgs.SetPtEtaPhiM(pthiggs,etahiggs,phihiggs,mhiggs);
  
  //boosting the two tau vectors
  double bx(m_higgs.Px()/m_higgs.E());
  double by(m_higgs.Py()/m_higgs.E());
  double bz(m_higgs.Pz()/m_higgs.E());
    HHLorentzVector boost1(m_tau1);
  HHLorentzVector boost2(m_tau2);
  boost1.Boost(bx,by,bz);
  m_tau1boosted=boost1; // why  m_tau1boosted(boost1); dont work?
  boost2.Boost(bx,by,bz);
  m_tau2boosted=boost2;
  //generate lorenzvector for tauvis
  //tau1vis
  HHLorentzVector vis1 =m_tau1boosted;
  vis1.SetPx(m_visfrac1*m_tau1boosted.Px());
  vis1.SetPy(m_visfrac1*m_tau1boosted.Py());
  vis1.SetPz(m_visfrac1*m_tau1boosted.Pz());
  m_tau1vis=vis1;
  //tau2vis
    HHLorentzVector vis2 =m_tau2boosted;
    vis2.SetPx(m_visfrac2*m_tau2boosted.Px());
    vis2.SetPy(m_visfrac2*m_tau2boosted.Py());
    vis2.SetPz(m_visfrac2*m_tau2boosted.Pz());
    m_tau2vis=vis2;
  //find MET
  double METx=-(m_tau1vis.Px()+m_tau2vis.Px()+m_isr.Px());
  std::cout << METx << std::endl;
  double METy=-(m_tau1vis.Py()+m_tau2vis.Py()+m_isr.Py());
  TVector2 met(METx,METy);
  m_MET=met;
  std::cout << "met test" << std::endl;
  m_MET.Print();
  }

HHLorentzVector HHTauTauEventGenerator::getTau1(){
  return(m_tau1);
}

HHLorentzVector HHTauTauEventGenerator::getTau2(){
  return(m_tau2);
}

HHLorentzVector  HHTauTauEventGenerator::getTau1boosted(){
  return(m_tau1boosted);
}

HHLorentzVector  HHTauTauEventGenerator::getTau2boosted(){
  return(m_tau2boosted);
}

HHLorentzVector  HHTauTauEventGenerator::getTau1Vis(){
  return(m_tau1vis);
}

HHLorentzVector  HHTauTauEventGenerator::getTau2Vis(){
  return(m_tau2vis);
}

TVector2  HHTauTauEventGenerator::getMET(){
  return(m_MET);
}





