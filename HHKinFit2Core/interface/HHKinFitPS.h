/*
 * class for general fit constraint
 */

#ifndef HHKinFit_
#define HHKinFit_

#include "TGraph.h"
#include "TGraph2D.h"
#include <vector>

#include "HHKinFit2/HHKinFit2Core/interface/HHFitObjectE.h"
#include "HHKinFit2/HHKinFit2Core/interface/HHFitConstraint.h"

namespace HHKinFit2{
  namespace PS{
class HHKinFit {
 public:
  HHKinFit();

  void setPrintLevel(int printlevel=0);

  void fit();
  double getChi2(bool respectLimits=true) const;
  double getL(bool respectLimits=true) const;
  std::vector<HHFitObjectE*> getListOfFitObjects() const;
  std::vector<HHFitConstraint*> getListOfConstraints() const;
  
  void addFitObjectE(HHFitObjectE* fitobject);
  void addConstraint(HHFitConstraint* constraint);
  
  //Returns TGraph for 1D Fit
  TGraph* getChi2Function(int steps);
  ///Returns TGraph2D for 2D Fit
  TGraph2D* getChi2Function(int* steps, double* mins, double* maxs);
  TGraph* getLFunction(int steps);
  int getConvergence() const;

  void printChi2() const;

  int m_loopsNeeded;
 private:
  std::vector<HHFitObjectE*> m_fitobjects;
  std::vector<HHFitConstraint*> m_constraints;

  double m_chi2;
  int m_convergence;
  int m_printlevel;
  int m_maxloops;
};
}
}
#endif /* HHKinFit_ */
