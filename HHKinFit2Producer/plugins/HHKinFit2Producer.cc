// -*- C++ -*-
//
// Package:    HHKinFit2/HHKinFit2Producer
// Class:      HHKinFit2Producer
// 
/**\class HHKinFit2Producer HHKinFit2Producer.cc HHKinFit2/HHKinFit2Producer/plugins/HHKinFit2Producer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Benedikt Vormwald
//         Created:  Tue, 23 Aug 2016 14:30:11 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "HHKinFit2/HHKinFit2/interface/HHKinFitMasterHeavyHiggs.h"

#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TVector2.h>

//
// class declaration
//

class HHKinFit2Producer : public edm::stream::EDProducer<> {
public:
  explicit HHKinFit2Producer(const edm::ParameterSet&);
  ~HHKinFit2Producer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;
  
  edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
  edm::EDGetTokenT<pat::TauCollection> tausToken_;
  edm::EDGetTokenT<pat::JetCollection> jetsToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  bool debug_;
};




HHKinFit2Producer::HHKinFit2Producer(const edm::ParameterSet& iConfig){
  produces<double>("mH").setBranchAlias("mH");
  produces<double>("P").setBranchAlias("P");
  produces<double>("chi2").setBranchAlias("chi2");
  produces<int>("convergence").setBranchAlias("convergence");

  muonsToken_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  tausToken_  = consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"));
  jetsToken_  = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("bjets"));
  metToken_   = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("met"));
}


HHKinFit2Producer::~HHKinFit2Producer(){
}


void
HHKinFit2Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  edm::Handle<pat::MuonCollection> muons;
  edm::Handle<pat::TauCollection> taus;
  edm::Handle<pat::JetCollection> jets;
  edm::Handle<pat::METCollection> met;

  iEvent.getByToken(muonsToken_,muons);
  iEvent.getByToken(tausToken_,taus);
  iEvent.getByToken(jetsToken_,jets);
  iEvent.getByToken(metToken_,met);

  double mH = 0;
  double prob = 0;
  double chi2 = 0;
  int convergence = -10;

  if ((muons->size()!=1)||(taus->size()!=1)){
    std::cout << "WARNING: At least one of the used lepton input collections has more than 1 entry! The fit only uses the first element!" << std::endl;
  }
    
  if (jets->size()>2){
    std::cout << "WARNING: The jet input collections has more than 2 entry! The fit only uses the first two element!" << std::endl;
  }

  if ((jets->size()<2)||(muons->size()<1)||(taus->size()<1)||(met->size()<1)){
    std::cout << "ERROR: At least one input collection has not enough entries! Event skipped!" << std::endl;
  }
  else{

    TLorentzVector bjet1(jets->at(0).px(),jets->at(0).py(),jets->at(0).pz(),jets->at(0).energy());
    TLorentzVector bjet2(jets->at(1).px(),jets->at(1).py(),jets->at(1).pz(),jets->at(1).energy());
    TLorentzVector tau1(muons->at(0).px(),muons->at(0).py(),muons->at(0).pz(),muons->at(0).energy());
    TLorentzVector tau2(taus->at(0).px(),taus->at(0).py(),taus->at(0).pz(),taus->at(0).energy());

    TVector2 met_vec(met->at(0).px(),met->at(0).py());
    TMatrixD met_cov(2,2);
    met_cov[0][0]=100;
    met_cov[0][1]=0;
    met_cov[1][0]=0;
    met_cov[1][1]=100;

    HHKinFit2::HHKinFitMasterHeavyHiggs* kinFit = new HHKinFit2::HHKinFitMasterHeavyHiggs(bjet1,bjet2,tau1,tau2,met_vec,met_cov);

    kinFit->useAdvancedBJetChi2(true);
    kinFit->addHypo(125, 125);
    kinFit->fit();
    
    // TLorentzVector h1 = kinFit->getFittedH1(125, 125);
    // TLorentzVector h2 = kinFit->getFittedH2(125, 125);
    // TLorentzVector b1 = kinFit->getFittedBJet1(125, 125);
    // TLorentzVector b2 = kinFit->getFittedBJet2(125, 125);
    // TLorentzVector tau1 = kinFit->getFittedTau1(125, 125);
    // TLorentzVector tau2 = kinFit->getFittedTau2(125, 125);
    
    mH           = kinFit->getMH(125, 125);
    chi2         = kinFit->getChi2(125, 125);
    convergence  = kinFit->getConvergence(125, 125);
    prob         = kinFit->getFitProb(125, 125);
  }

  std::auto_ptr<double> pmH(new double(mH));
  std::auto_ptr<double> pprob(new double(prob));
  std::auto_ptr<double> pchi2(new double(chi2));
  std::auto_ptr<int> pconvergence(new int(convergence));

  iEvent.put(std::move(pmH),  "mH");
  iEvent.put(std::move(pprob),   "P");
  iEvent.put(std::move(pchi2),"chi2");
  iEvent.put(std::move(pconvergence),"convergence");
}


void
HHKinFit2Producer::beginStream(edm::StreamID){
}


void
HHKinFit2Producer::endStream() {
}


void
HHKinFit2Producer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("muons")->setComment("muon input collection");
  desc.add<edm::InputTag>("taus")->setComment("tau input collection");
  desc.add<edm::InputTag>("bjets")->setComment("jet input collection");
  desc.add<edm::InputTag>("met")->setComment("met input collection");
  desc.add<bool>("debug", false)->setComment("flag to switch on debug output");

  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(HHKinFit2Producer);
