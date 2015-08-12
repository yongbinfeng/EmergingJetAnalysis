// -*- C++ -*-
//
// Package:    EmergingAnalysis/WJetFilter
// Class:      WJetFilter
// 
/**\class WJetFilter WJetFilter.cc EmergingAnalysis/WJetFilter/plugins/WJetFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Young Ho Shin
//         Created:  Mon, 10 Aug 2015 23:22:24 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TVector3.h>
#include <TMath.h>
#include <Math/VectorUtil.h>

// For creating/writing histograms
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"

// Data formats
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

// Namespace shorthands
using std::string;
using std::vector;

//
// class declaration
//

class WJetFilter : public edm::EDFilter {
  public:
    explicit WJetFilter(const edm::ParameterSet&);
    ~WJetFilter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() override;
    virtual bool filter(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------
    string alias_; // Alias suffix for all products
    // Retrieve once per event
    edm::EDGetTokenT< pat::MuonCollection > muonCollectionToken_;
    edm::EDGetTokenT< pat::ElectronCollection > electronCollectionToken_;
    edm::EDGetTokenT< pat::JetCollection > jetCollectionToken_;
    edm::EDGetTokenT< pat::METCollection > metCollectionToken_;
    double minPtMuon_;
    double minPtElectron_;
    double minPtMET_;
    double maxDeltaPhi_;
    double minPtSelectedJet_;
    double maxPtAdditionalJets_;

    std::unordered_map<string, TH1*> histoMap1D_;
    std::unordered_map<string, TH2*> histoMap2D_;
};

//
// constants, enums and typedefs
//
/// Lorentz vector
typedef reco::Candidate::LorentzVector RecoLorentzVector;
/// Lorentz vector
typedef math::XYZTLorentzVector LorentzVector;
/// Lorentz vector
typedef math::PtEtaPhiMLorentzVector PolarLorentzVector;

//
// static data member definitions
//

//
// constructors and destructor
//
WJetFilter::WJetFilter(const edm::ParameterSet& iConfig) :
    minPtMuon_           (  iConfig.getParameter<double>("minPtMuon") ),
    minPtElectron_       (  iConfig.getParameter<double>("minPtElectron") ),
    minPtMET_            (  iConfig.getParameter<double>("minPtMET") ),
    maxDeltaPhi_         (  iConfig.getParameter<double>("maxDeltaPhi") ),
    minPtSelectedJet_    (  iConfig.getParameter<double>("minPtSelectedJet") ),
    maxPtAdditionalJets_ (  iConfig.getParameter<double>("maxPtAdditionalJets") )
{
   //now do what ever initialization is needed

  alias_ = iConfig.getParameter<string>("@module_label");

  muonCollectionToken_ = consumes< pat::MuonCollection > (iConfig.getParameter<edm::InputTag>("srcMuons"));
  electronCollectionToken_ = consumes< pat::ElectronCollection > (iConfig.getParameter<edm::InputTag>("srcElectrons"));
  jetCollectionToken_ = consumes< pat::JetCollection > (iConfig.getParameter<edm::InputTag>("srcJets"));
  metCollectionToken_ = consumes< pat::METCollection > (iConfig.getParameter<edm::InputTag>("srcMET"));

  //Register products
  // produces< bool > ("zValidity")  . setBranchAlias( string("zValidity_").append(alias_) );
  // produces< PolarLorentzVector > ("zP4")  . setBranchAlias( string("zP4_").append(alias_) );
  //
  // produces< pat::JetCollection > ("jetSelected")  . setBranchAlias( string("jetSelected_").append(alias_) );
  // produces< vector<double> > ("deltaR")  . setBranchAlias( string("deltaR_").append(alias_) );
  // produces< vector<double> > ("deltaPhi")  . setBranchAlias( string("deltaPhi_").append(alias_) );
  // produces< bool > ("eventPassed")  . setBranchAlias( string("eventPassed_").append(alias_) );

  edm::Service<TFileService> fs;
  string name;
  name="mT_e"  ; histoMap1D_.emplace( name , fs->make<TH1D>(name.c_str() , name.c_str() , 100 , 0., 200.) );
  name="mT_mu" ; histoMap1D_.emplace( name , fs->make<TH1D>(name.c_str() , name.c_str() , 100 , 0., 200.) );
}


WJetFilter::~WJetFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
WJetFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  bool eventPassed = false;
  edm::Handle< pat::MuonCollection > muonCollection;
  iEvent.getByToken(muonCollectionToken_, muonCollection);
  edm::Handle< pat::ElectronCollection > electronCollection;
  iEvent.getByToken(electronCollectionToken_, electronCollection);
  edm::Handle< pat::JetCollection > jetCollection;
  iEvent.getByToken(jetCollectionToken_, jetCollection);
  edm::Handle< pat::METCollection > metCollection;
  iEvent.getByToken(metCollectionToken_, metCollection);
  auto met = (metCollection->at(0));
  if ( met.pt() < minPtMET_ ) return eventPassed;

  ////////////////////////////////////////////////////////////
  // Calculate mT with electrons/muons
  ////////////////////////////////////////////////////////////

  // Construct zP4 from electrons
  {
    // Select good leptons
    auto leptons = electronCollection.product();
    int nGoodLepton = 0;
    pat::ElectronCollection goodLeptons;
    // Lepton-specific criteria
    for ( auto it = leptons->begin(); it != leptons->end(); it++ ) {
      auto lepton = *it;
      if ( lepton.pt() < minPtElectron_ ) continue;
      if ( abs(lepton.eta()) > 2.5 ) continue;
      if ( lepton.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-medium") < 1 ) continue;
      nGoodLepton++;
      goodLeptons.push_back(lepton);
    }

    if ( goodLeptons.size() >= 1 ) {
      float dPhi = ROOT::Math::VectorUtil::DeltaPhi( met.p4(), goodLeptons[0].p4() );
      double mT = TMath::Sqrt( 2. * goodLeptons[0].pt() * met.pt() * ( 1 - TMath::Cos(dPhi) ) );
      LogTrace("ZJetFilter") << "mT_e: " << mT;
      histoMap1D_["mT_e"]->Fill(mT);
    }
  }

  {
    // Select good leptons
    auto leptons = muonCollection.product();
    int nGoodLepton = 0;
    pat::MuonCollection goodLeptons;
    // Lepton-specific criteria
    for ( auto it = leptons->begin(); it != leptons->end(); it++ ) {
      auto lepton = *it;
      if ( lepton.pt() < minPtElectron_ ) continue;
      if ( abs(lepton.eta()) > 2.5 ) continue;
      if ( !lepton.muonID("AllGlobalMuons") ) continue;
      nGoodLepton++;
      goodLeptons.push_back(lepton);
    }

    if ( goodLeptons.size() >= 1 ) {
      float dPhi = ROOT::Math::VectorUtil::DeltaPhi( met.p4(), goodLeptons[0].p4() );
      double mT = TMath::Sqrt( 2. * goodLeptons[0].pt() * met.pt() * ( 1 - TMath::Cos(dPhi) ) );
      LogTrace("ZJetFilter") << "mT_mu " << mT;
      histoMap1D_["mT_mu"]->Fill(mT);
    }
  }

  // LogTrace("ZJetFilter") << "Printing electron pts.";
  // for ( auto it = electrons->begin(); it != electrons->end(); it++ ) {
  //   auto electron = *it;
  //   LogTrace("ZJetFilter") << electron.pt();
  // }


  // std::auto_ptr< bool > _zValidity( new bool (zValidity) );
  // iEvent.put(_zValidity, "zValidity");
  // std::auto_ptr< PolarLorentzVector > _zP4( new PolarLorentzVector (zP4) );
  // iEvent.put(_zP4, "zP4");

  // LogTrace("ZJetFilter") << "eventPassed: " << eventPassed;

  // std::auto_ptr< reco::PFJetCollection > _jetToSave( new reco::PFJetCollection (jetToSave) );
  // iEvent.put(_jetToSave, "jetSelected");
  // std::auto_ptr< vector<double> > _deltaRs( new vector<double> (deltaRs) );
  // iEvent.put(_deltaRs, "deltaR");
  // std::auto_ptr< vector<double> > _deltaPhis( new vector<double> (deltaPhis) );
  // iEvent.put(_deltaPhis, "deltaPhi");
  // std::auto_ptr< bool > _eventPassed( new bool (eventPassed) );
  // iEvent.put(_eventPassed, "eventPassed");

  return eventPassed;
}

// ------------ method called once each job just before starting event loop  ------------
void 
WJetFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
WJetFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
WJetFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
WJetFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
WJetFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
WJetFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
WJetFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(WJetFilter);
