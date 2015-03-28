#include <iostream>
#include "HHKinFit.h"
#include "HHFitConstraint4Vector.h"
#include "HHFitConstraintEHardM.h"
#include "HHFitConstraint.h"
#include "HHFitConstraintLikelihood.h"
#include "HHFitObjectEConstM.h"
#include "HHFitObjectE.h"
#include "HHFitObjectMET.h"
#include "HHFitObject.h"
#include "HHFitObjectComposite.h"
#include "exceptions/HHEnergyRangeException.h"
#include "exceptions/HHCovarianceMatrixException.h"

int main(int argc, char* argv[])
{
  double mass = 125.4;

  HHFitObjectE* tau1 = new HHFitObjectEConstM(TLorentzVector(58,0,0,60));
  HHFitObjectE* tau2 = new HHFitObjectEConstM(TLorentzVector(0,58,0,60));
  HHFitObjectMET* met = new HHFitObjectMET(TVector2(10,20),TVector2(0,0));
  met->setCovMatrix(100,-100,50);

  std::vector<HHFitObject*> higgsparts;
  higgsparts.push_back(tau1);
  higgsparts.push_back(tau2);
  higgsparts.push_back(met);
  HHFitObject* higgs  = new HHFitObjectComposite(higgsparts);

  HHFitConstraint* invm = new HHFitConstraintEHardM(tau1, tau2, mass);
  HHFitConstraint* balance = new HHFitConstraint4Vector(higgs, true, true, false, false);

  HHKinFit* fit = new HHKinFit();
  fit->addFitObject(tau1);
  fit->addConstraint(invm);
  fit->addConstraint(balance);

  try{
    std::cout << fit->getChi2() << std::endl;
  }
  catch(const HHCovarianceMatrixException& e){
    std::cout << e.what() << std::endl;
    std::cout << "will fix it manually" << std::endl;
    met->setCovMatrix(100,100,0);
    std::cout << fit->getChi2() << std::endl;
  }
  catch(...){
    std::cout << "caught an unexpected exception" << std::endl;
  }
  return (0);
}
