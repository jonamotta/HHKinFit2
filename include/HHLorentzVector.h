#include "TLorentzVector.h"

#ifndef HHLorentzVector_
#define HHLorentzVector_

class HHLorentzVector : public TLorentzVector {

public:
	HHLorentzVector(double x, double y, double z, double t);
	HHLorentzVector();
	void SetEEtaPhiM(double E,double Eta,double Phi,double M);

	HHLorentzVector operator+(HHLorentzVector const& rhs) const;
	HHLorentzVector operator-(HHLorentzVector const& rhs) const;
};

#endif /* HHLorentzVector_ */
