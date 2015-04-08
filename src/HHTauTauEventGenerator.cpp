#include "HHTauTauEventGenerator.h"
#include "TRandom.h"
#include <cmath>
#include <iostream>
HHTauTauEventGenerator::HHTauTauEventGenerator()
  : m_tau1(),
    m_tau2() {

// Generate lorentzvectors for the two tau in the rest frame of the higgs
  
//Tau 1
  double pi(3.14159265359);
  double mhiggs(125.7);
  double mtau1(1.77682);
  double mtau2(1.77682);
  TRandom randomnumber;
  double cthtau1(randomnumber.Uniform(-1,1)); // need a generator for random numbers
  double etatau1(log(tan(acos(cthtau1)/2)));
  double phitau1(randomnumber.Uniform(0,2*pi));
  double etau1((mhiggs*mhiggs-mtau2*mtau2 +mtau1*mtau1)/(2*mhiggs));
  m_tau1.SetEEtaPhiM(etau1, etatau1, phitau1, mtau1);
  double etau2((mhiggs*mhiggs-mtau1*mtau1 +mtau2*mtau2)/(2*mhiggs));
  double etatau2(-etatau1);
  double phitau2(phitau1+pi);
  m_tau2.SetEEtaPhiM(etau2, etatau2, phitau2, mtau2);
}

HHLorentzVector HHTauTauEventGenerator::getTau1(){
  return(m_tau1);
}

HHLorentzVector HHTauTauEventGenerator::getTau2(){
  return(m_tau2);
}


