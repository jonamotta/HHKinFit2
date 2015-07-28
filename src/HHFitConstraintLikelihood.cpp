#include "HHFitConstraintLikelihood.h"
#include <cmath>
#include <iostream>
#include "TSpline.h"

HHKinFit2::HHFitConstraintLikelihood::HHFitConstraintLikelihood(HHFitObject* object1,HHFitObject* object2, TF1* likelihood1,TF1* likelihood2 )
  : HHFitConstraint(object1),
    m_likelihood1(likelihood1),
    m_likelihood2(likelihood2),
    m_likelihoodhisto1(),
    m_likelihoodhisto2(),
    m_likelihoodhisto(),
    m_likelihoodspline1(),
    m_likelihoodspline2(),
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
    m_likelihoodspline1(),
    m_likelihoodspline2(),
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
    m_likelihoodspline1(),
    m_likelihoodspline2(),
    mode(3),
    m_object2(object2){
}

HHKinFit2::HHFitConstraintLikelihood::HHFitConstraintLikelihood(HHFitObject* object1,HHFitObject* object2,TSpline3* likelihoodspline1,  TSpline3* likelihoodspline2)
  : HHFitConstraint(object1),
    m_likelihood1(),
    m_likelihood2(),
    m_likelihoodhisto1(),
    m_likelihoodhisto2(),
    m_likelihoodhisto(),
    m_likelihoodspline1(likelihoodspline1),
    m_likelihoodspline2(likelihoodspline2),
    mode(4),
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
		TSpline3* temp1=new TSpline3(m_likelihoodhisto1);
		TSpline3* temp2=new TSpline3(m_likelihoodhisto2);
		double x1=temp1->Eval(m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E());
		double x2=temp2->Eval(m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E());
		if(x1<0)x1=0.0000001;
		if(x2<0)x2=0.0000001;
		delete(temp1);
		delete(temp2);
  double chi2=-2*log(x1)-2*log(x2);
  return(chi2);
	}
	if(mode==3){
		double chi2=-2*log(m_likelihoodhisto->GetBinContent(m_likelihoodhisto->FindBin(m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E(),m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E() )));
		  return(chi2);
	}
	if(mode==4){
		double x1=m_likelihoodspline1->Eval(m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E());
		double x2=m_likelihoodspline2->Eval(m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E());
		if(x1<0)x1=0.0000001;
	    if(x2<0)x2=0.0000001;
		double chi2=-2*log(x1)-2*log(x2);
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
		TSpline3* temp1=new TSpline3(m_likelihoodhisto1);
		TSpline3* temp2=new TSpline3(m_likelihoodhisto2);
		double x1=temp1->Eval(m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E());
		double x2=temp2->Eval(m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E());
        return(x1*x2);
	}
	if(mode==3){
		return(m_likelihoodhisto->GetBinContent(m_likelihoodhisto->FindBin(m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E(),m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E() )));
	}
	if(mode==4){
		double x1=m_likelihoodspline1->Eval(m_fitobject->getInitial4Vector().E()/m_fitobject->getFit4Vector().E());
	    double x2=m_likelihoodspline2->Eval(m_object2->getInitial4Vector().E()/m_object2->getFit4Vector().E());
	    if(x1<0)x1=0.0000001;
	    if(x2<0)x2=0.0000001;
	    return(x1*x2);
	}

}
