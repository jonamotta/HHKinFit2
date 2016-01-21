/*
 * class for constraint on 4 momentum conservation
 */

#ifndef HHFitConstraint4VectorBJet_
#define HHFitConstraint4VectorBJet_

#ifdef HHKINFIT2
#include "HHFitConstraint.h"
#include "HHFitObject.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHFitConstraint.h"
#include "HHKinFit2/HHKinFit2/interface/HHFitObject.h"
#endif

#include <vector>

namespace HHKinFit2{
  double crystBallLikeCDF(double x, double alpha, double n, double sigma, double mean,
			  double beta, double normalization);

class HHFitConstraint4VectorBJet : public HHFitConstraint {
 public:
  HHFitConstraint4VectorBJet(HHFitObject* object);
  HHFitConstraint4VectorBJet(HHFitObject* object, double alpha, double n, double sigma, 
			     double mean, double beta, double normalization);

  double getChi2() const;
  double getLikelihood() const;

  void printChi2() const;

  double m_alpha;
  double m_n;
  double m_sigma;
  double m_mean;
  double m_beta;
  double m_normalization;  
};
}
#endif /* HHFitConstraint4VectorBJet_ */
