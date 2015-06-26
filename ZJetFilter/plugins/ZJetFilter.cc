// -*- C++ -*-
//
// Package:    EmergingJetAnalysis/ZJetFilter
// Class:      ZJetFilter
// 
/**\class ZJetFilter ZJetFilter.cc EmergingJetAnalysis/ZJetFilter/plugins/ZJetFilter.cc

 Description: Selects Z+Jet events where there is exactly one hard jet back to back with the Z boson.

 Implementation:
     // Requires Z boson 4-vector to be present in the event. (reco::Candidate::LorentzVector)
     Saves selected hard jet.
*/
//
// Original Author:  Young Ho Shin
//         Created:  Mon, 08 Jun 2015 15:35:15 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TVector3.h>
#include <Math/VectorUtil.h>

// Data formats
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"

// Namespace shorthands
using std::string;
using std::vector;

//
// class declaration
//

class ZJetFilter : public edm::EDFilter {
   public:
      explicit ZJetFilter(const edm::ParameterSet&);
      ~ZJetFilter();

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
			// Retrieve once 
			string alias_; // Alias suffix for all products
			// Retrieve once per event
			edm::EDGetTokenT< edm::View<reco::Candidate> > muonCollectionToken_;
			edm::EDGetTokenT< pat::ElectronCollection > electronCollectionToken_;
			edm::EDGetTokenT< reco::PFJetCollection > jetCollectionToken_;
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
ZJetFilter::ZJetFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

	alias_ = iConfig.getParameter<string>("@module_label");

	muonCollectionToken_ = consumes< edm::View<reco::Candidate> > (iConfig.getParameter<edm::InputTag>("srcMuons"));
	electronCollectionToken_ = consumes< pat::ElectronCollection > (iConfig.getParameter<edm::InputTag>("srcElectrons"));
	jetCollectionToken_ = consumes< reco::PFJetCollection > (iConfig.getParameter<edm::InputTag>("srcJets"));

  //Register products
	produces< bool > ("zValidity")  . setBranchAlias( string("zValidity_").append(alias_) );
	produces< PolarLorentzVector > ("zP4")  . setBranchAlias( string("zP4_").append(alias_) );

	produces< reco::PFJetCollection > ("jetSelected")  . setBranchAlias( string("jetSelected_").append(alias_) );
	produces< vector<double> > ("deltaR")  . setBranchAlias( string("deltaR_").append(alias_) );
	produces< vector<double> > ("deltaPhi")  . setBranchAlias( string("deltaPhi_").append(alias_) );
	produces< bool > ("eventPassed")  . setBranchAlias( string("eventPassed_").append(alias_) );
}


ZJetFilter::~ZJetFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ZJetFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

	edm::Handle< edm::View<reco::Candidate> > muonCollection;
	iEvent.getByToken(muonCollectionToken_, muonCollection);
	edm::Handle< pat::ElectronCollection > electronCollection;
	iEvent.getByToken(electronCollectionToken_, electronCollection);
	edm::Handle< reco::PFJetCollection > jetCollection;
	iEvent.getByToken(jetCollectionToken_, jetCollection);

  ////////////////////////////////////////////////////////////
  // Calculate zP4 with electrons/muons
  ////////////////////////////////////////////////////////////
	bool zValidity = false;
	PolarLorentzVector zP4;

	auto electrons = electronCollection.product();
	// int electronCount = 0;
	// for (auto electron = electrons->begin(); electron != electrons->end(); electron++) {
	// 	if ( electron->electronID("eidLoose") > 0 ) electronCount++;
	// }
	if ( electrons->size() >= 2 ) {
		RecoLorentzVector zP4_;
		zP4_ = electrons->at(0).p4() + electrons->at(1).p4();
		if ( 80 < zP4_.mass() && zP4_.mass() < 100 ) {
			zP4 = zP4_;
			zValidity = true;
		}
	}

	auto muons = muonCollection.product();
	if ( muons->size() >= 2 ) {
		RecoLorentzVector zP4_;
		zP4_ = muons->at(0).p4() + muons->at(1).p4();
		if ( 80 < zP4_.mass() && zP4_.mass() < 100 ) {
			zP4 = zP4_;
			zValidity = true;
		}
	}

	std::auto_ptr< bool > _zValidity( new bool (zValidity) );
	iEvent.put(_zValidity, "zValidity");
	std::auto_ptr< PolarLorentzVector > _zP4( new PolarLorentzVector (zP4) );
	iEvent.put(_zP4, "zP4");

  ////////////////////////////////////////////////////////////
  // Check if there is a jet opposite to Z
  // Veto presence of additional hard jets
  ////////////////////////////////////////////////////////////
  bool eventPassed = false;
  reco::PFJetCollection jetToSave;
  vector<double> deltaRs;
  vector<double> deltaPhis;
  if (zValidity) {
    auto jets = jetCollection.product();
    reco::PFJetCollection::const_iterator jetSelected = jets->end();
    // Loop over jets to find first jet that is within dR cut of -zP4
    int iJet = 0;
    for ( auto jet = jets->begin(); jet!= jets->end(); jet++ ) {
      float dR = ROOT::Math::VectorUtil::DeltaR( jet->p4(), - zP4 );
      float dPhi = ROOT::Math::VectorUtil::DeltaPhi( jet->p4(), -zP4 );
      // std::cout << "dR: " << dR << " \n";
      // std::cout << "jet->pt(): " << jet->pt() << " \n";
      // std::cout << iJet << "\t" << dR << " \t" << jet->pt() << " \n";
      if ( dPhi < 0.4 && jet->pt() > 50 ) {
        eventPassed = true;
        jetSelected = jet;
        jetToSave.push_back(*jet);
        deltaRs.push_back(dR);
        deltaPhis.push_back(dPhi);
        break;
      }
      iJet++;
    }
    // Veto presence of additional hard jets
    for ( auto jet = jets->begin(); jet!= jets->end(); jet++ ) {
      if ( jet == jetSelected ) {
        continue;
      }
      else if ( jet->pt() > 50 ) {
        eventPassed = false;
        break;
      }
    }
  }

	std::auto_ptr< reco::PFJetCollection > _jetToSave( new reco::PFJetCollection (jetToSave) );
	iEvent.put(_jetToSave, "jetSelected");
	std::auto_ptr< vector<double> > _deltaRs( new vector<double> (deltaRs) );
	iEvent.put(_deltaRs, "deltaR");
	std::auto_ptr< vector<double> > _deltaPhis( new vector<double> (deltaPhis) );
	iEvent.put(_deltaPhis, "deltaPhi");
	std::auto_ptr< bool > _eventPassed( new bool (eventPassed) );
	iEvent.put(_eventPassed, "eventPassed");

	return eventPassed;
}

// ------------ method called once each job just before starting event loop  ------------
void 
ZJetFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZJetFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
ZJetFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ZJetFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ZJetFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ZJetFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZJetFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ZJetFilter);
