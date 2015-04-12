#include "HHTauTauEventGenerator.h"
#include <cmath>
#include <iostream>
#include "TVector3.h"
HHTauTauEventGenerator::HHTauTauEventGenerator(TF1 a,TF1 b){
  m_PDF1=a;
  m_PDF2=b;
  m_randomnumber;
  m_tau1;
  m_tau2;
  m_isr;
  m_higgs;
  m_tau1boosted;
  m_tau2boosted;
  m_visfrac1;
  m_visfrac2;
  m_tau1vis;
  m_tau2vis;
  m_MET;
  m_mhiggs=125.7;
  m_mtau1=1.77682;
  m_mtau2=1.77682;
  m_pi=3.14159265359;
  m_misr=0;
}
void HHTauTauEventGenerator::generateEvent() {
  // Generate lorentzvectors for the two tau in the rest frame of the higgs
  
  double cthtau1=m_randomnumber.Uniform(-1,1);
  double etatau1=log(tan(acos(cthtau1)/2));
  double phitau1=m_randomnumber.Uniform(0,2*m_pi);
  double etau1=(m_mhiggs*m_mhiggs-m_mtau2*m_mtau2 +m_mtau1*m_mtau1)/(2*m_mhiggs);
  m_tau1.SetEEtaPhiM(etau1, etatau1, phitau1, m_mtau1);
  //Tau 2

  double etau2=(m_mhiggs*m_mhiggs-m_mtau1*m_mtau1 +m_mtau2*m_mtau2)/(2*m_mhiggs);
  double etatau2=-etatau1;
  double phitau2=phitau1+m_pi;
  m_tau2.SetEEtaPhiM(etau2, etatau2, phitau2, m_mtau2);
  //generate lorentzvectors for isr-jet

  double ptisr=fabs(m_randomnumber.Gaus(0,30)); // absolut Value of Pt
  double phiisr=m_randomnumber.Uniform(0,2*m_pi);
  double etaisr=m_randomnumber.Uniform(-5,5);
  m_isr.SetPtEtaPhiM(ptisr,etaisr,phiisr,m_misr);
  //generate lorentzvector for higgs with lorentzvector from isr-jet

  double cthhiggs=m_randomnumber.Uniform(-1,1);
  double etahiggs=log(tan(acos(cthhiggs)/2));
  double pthiggs=ptisr;
  double phihiggs=phiisr+m_pi;
  m_higgs.SetPtEtaPhiM(pthiggs,etahiggs,phihiggs,m_mhiggs);
  //boosting the two tau vectors

  double bx=m_higgs.Px()/m_higgs.E();
  double by=m_higgs.Py()/m_higgs.E();
  double bz=m_higgs.Pz()/m_higgs.E();
  m_tau1boosted=m_tau1;
  m_tau1boosted.Boost(bx,by,bz);
  m_tau2boosted=m_tau2;
  m_tau2boosted.Boost(bx,by,bz);
  //generate lorenzvector for tauvis
  //generate visfrac

  m_visfrac1=m_PDF1.GetRandom(0,1);
  m_visfrac2=m_PDF2.GetRandom(0,1);
  //tau1vis

  m_tau1vis=m_tau1boosted;
  m_tau1vis.SetEkeepM(m_visfrac1*m_tau1boosted.E());
  //tau2vis

  m_tau2vis=m_tau2boosted;
  m_tau2vis.SetEkeepM(m_visfrac2*m_tau2boosted.E());
  //find MET

  double METx=-(m_tau1vis.Px()+m_tau2vis.Px()+m_isr.Px());
  double METy=-(m_tau1vis.Py()+m_tau2vis.Py()+m_isr.Py());
  TVector2 met(METx,METy);
  m_MET=met;
  // modify MET with a covariance-matrix

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

double HHTauTauEventGenerator::getvisfrac1(){
  return(m_visfrac1);
}

double HHTauTauEventGenerator::getvisfrac2(){
  return(m_visfrac2);
}

void HHTauTauEventGenerator::setMhiggs(double M){
  m_mhiggs=M;
}

void HHTauTauEventGenerator::setMtau1(double M){
  m_mtau1=M;
}

void HHTauTauEventGenerator::setMtau2(double M){
  m_mtau2=M;
}





