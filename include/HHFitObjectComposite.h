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
	  HHFitObjectComposite(std::vector<HHFitObject*> const& subobjects);
	  HHFitObjectComposite(HHFitObject* subobject1, HHFitObject* subobject2, HHFitObject* subobject3);

 protected:
  TLorentzVector getFit4Vector() const;
  TLorentzVector getInitial4Vector() const;
  TMatrixD getCovMatrix() const;

  void setSubobjects(std::vector<HHFitObject*> const& subobjects);
  void addSubobject(HHFitObject* subobject);

  virtual void print() const;

 private:
  std::vector<HHFitObject*> m_subobjects;
};

#endif /* HHFitObjectComposite_ */
