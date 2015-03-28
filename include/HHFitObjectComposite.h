/*
 * class for composite objects
 */

#ifndef HHFitObjectComposite_
#define HHFitObjectComposite_

#include "Rtypes.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"

#include "HHFitObject.h"

class HHFitObjectComposite : public HHFitObject {
 public:
  HHFitObjectComposite(std::vector<HHFitObject*> subobjects);

 protected:
  TLorentzVector getFit4Vector();
  TLorentzVector getInitial4Vector();
  TMatrixD getCovMatrix();

  void setSubobjects(std::vector<HHFitObject*> subobjects);
  void addSubobject(HHFitObject* subobject);

  virtual void print();

 private:
  std::vector<HHFitObject*> m_subobjects;
};

#endif /* HHFitObjectComposite_ */
