/*
class for generating lorentzvectors for Tau in higgsdecay: Higgs->TauTau

*/

#include "HHLorentzVector.h"

 
class HHTauTauEventGenerator {

public:

HHTauTauEventGenerator();
HHLorentzVector getTau1boosted();
HHLorentzVector getTau2boosted();
 HHLorentzVector getTau1();
 HHLorentzVector getTau2();


private:
 HHLorentzVector m_tau1;
 HHLorentzVector m_tau2;
 HHLorentzVector m_isr;
 HHLorentzVector m_higgs;
 HHLorentzVector m_tau1boosted;
 HHLorentzVector m_tau2boosted;




};
