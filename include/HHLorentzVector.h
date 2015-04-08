#include "TLorentzVector.h"

class HHLorentzVector : public TLorentzVector {

public:
	HHLorentzVector(double x, double y, double z, double t);
	HHLorentzVector();
	void SetEEtaPhiM(double E,double Eta,double Phi,double M);
};
