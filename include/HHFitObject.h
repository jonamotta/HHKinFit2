/*
 * class for fit objects
 */

#ifndef HHFitObject_
#define HHFitObject_

#include "Rtypes.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"

class HHFitObject {
 public:
  HHFitObject();
  HHFitObject(TLorentzVector initial4vector);
  virtual ~HHFitObject(){};

  virtual TLorentzVector getInitial4Vector();
  virtual TLorentzVector getFit4Vector();
  virtual TMatrixD getCovMatrix();

  void setFit4Vector(TLorentzVector vec);
  void setCovMatrix(TMatrixD covmatrix);

  void reset();

  virtual void printInitial4Vector();
  virtual void printFit4Vector();
  virtual void printCovMatrix();

  virtual void print();

 protected:
  TLorentzVector m_fit4vector;
  TLorentzVector m_initial4vector;
  TMatrixD m_covmatrix;
};

#endif /* HHFitObject_ */
