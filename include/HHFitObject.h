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
  HHFitObject(TLorentzVector const& initial4vector);
  virtual ~HHFitObject(){};

  virtual TLorentzVector getInitial4Vector() const;
  virtual TLorentzVector getFit4Vector() const;
  virtual TMatrixD getCovMatrix() const;

  void setFit4Vector(TLorentzVector const& vec);
  void setCovMatrix(TMatrixD const& covmatrix);

  void reset();

  virtual void printInitial4Vector() const;
  virtual void printFit4Vector() const;
  virtual void printCovMatrix() const;

  virtual void print() const;

 protected:
  TLorentzVector m_fit4vector;
  TLorentzVector const m_initial4vector;
  TMatrixD m_covmatrix;
};

#endif /* HHFitObject_ */
