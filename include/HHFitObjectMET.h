/*
 * class for fit objects
 */

#ifndef HHFitObjectMET_
#define HHFitObjectMET_

#include "Rtypes.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TMatrixD.h"
#include "HHFitObject.h"


class HHFitObjectMET : public HHFitObject {
 public:
  HHFitObjectMET(TVector2 initialvector,TVector2 fitvector);
  
  void setCovMatrix(Double_t xx, Double_t yy, Double_t xy);

  virtual void print();
  virtual void printInitial4Vector();
  virtual void printFit4Vector();
  virtual void printCovMatrix();
};

#endif /* HHFitObjectMET_ */
