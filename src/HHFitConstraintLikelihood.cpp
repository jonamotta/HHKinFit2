#include "HHFitConstraintLikelihood.h"
#include <cmath>
#include <iostream>


HHKinFit2::HHFitConstraintLikelihood::HHFitConstraintLikelihood(HHFitObject* object1,HHFitObject* object2, TF1* likelihood1,TF1* likelihood2 )
  : HHFitConstraint(object1),
    m_likelihood1(likelihood1),
    m_likelihood2(likelihood2),
    m_likelihoodhisto1(),
    m_likelihoodhisto2(),
    m_likelihoodhisto(),
    mode(1),
    m_object2(object2){

}

HHKinFit2::HHFitConstraintLikelihood::HHFitConstraintLikelihood(HHFitObject* object1,HHFitObject* object2, TH1D* likelihoodhisto1, TH1D* likelihoodhisto2 )
  : HHFitConstraint(object1),
    m_likelihood1(),
    m_likelihood2(),
    m_likelihoodhisto1(likelihoodhisto1),
    m_likelihoodhisto2(likelihoodhisto2),
    m_likelihoodhisto(),
    mode(2),
    m_object2(object2){
}

HHKinFit2::HHFitConstraintLikelihood::HHFitConstraintLikelihood(HHFitObject* object1,HHFitObject* object2, TH2D* likelihoodhisto)
  : HHFitConstraint(object1),
    m_likelihood1(),
    m_likelihood2(),
    m_likelihoodhisto1(),
    m_likelihoodhisto2(),
    m_likelihoodhisto(likelihoodhisto),
    mode(3),
    m_object2(object2){
}
/*Double_t
HHKinFit2::HHFitConstraintLikelihood::getChi2() const{
	double chi2=-2*log(m_likelihood1->Eval(m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E()))-2*log(m_likelihood2->Eval(m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E()));
  return(chi2);
}*

*/

Double_t 
HHKinFit2::HHFitConstraintLikelihood::getChi2() const{
	if(mode==1){
	double chi2=-2*log(m_likelihood1->Eval(m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E()))-2*log(m_likelihood2->Eval(m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E()));
    return(chi2);
	}
	if(mode==2){
  double chi2=-2*log(m_likelihoodhisto1->GetBinContent(m_likelihoodhisto1->FindBin(m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E())))-2*log(m_likelihoodhisto2->GetBinContent(m_likelihoodhisto2->FindBin(m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E())));
  return(chi2);
	}
	if(mode==3){
		double chi2=-2*log(m_likelihoodhisto->GetBinContent(m_likelihoodhisto->FindBin(m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E(),m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E() )));
		  return(chi2);
	}
}
/*double
HHKinFit2::HHFitConstraintLikelihood::getLikelihood() const{
  return(m_likelihood1->Eval(m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E())*m_likelihood2->Eval(m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E()));
}*/

double
HHKinFit2::HHFitConstraintLikelihood::getLikelihood() const{
	if(mode==1){
		return(m_likelihood1->Eval(m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E())*m_likelihood2->Eval(m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E()));
	}
	if(mode==2){
        return(m_likelihoodhisto1->GetBinContent(m_likelihoodhisto1->FindBin(m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E()))*m_likelihoodhisto2->GetBinContent(m_likelihoodhisto2->FindBin(m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E())));
	}
	if(mode==3){
		return(m_likelihoodhisto->GetBinContent(m_likelihoodhisto->FindBin(m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E(),m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E() )));
	}

}
