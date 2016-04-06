// -*- C++ -*-
//
// Package:    EmergingJetAnalysis/EmJetAnalyzer
// Class:      EmJetAnalyzer
//
/**\class EmJetAnalyzer EmJetAnalyzer.cc EmergingJetAnalysis/EmJetAnalyzer/plugins/EmJetAnalyzer.cc

 Description: Analyzer for Emerging Jet analysis, supercedes EmergingJetAnalyzer

 Implementation:
     Re-write of EmergingJetAnalyzer to be more modular.
*/
//
// Original Author:  Young Ho Shin
//         Created:  Tue, 05 Apr 2016 19:37:25 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "EmergingJetAnalysis/EmJetAnalyzer/interface/OutputTree.h"
#include "EmergingJetAnalysis/EmJetAnalyzer/interface/EmJetEvent.h"

//
// class declaration
//

class EmJetAnalyzer : public edm::EDFilter {
  public:
    explicit EmJetAnalyzer(const edm::ParameterSet&);
    ~EmJetAnalyzer();

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
    emjet:: Event  event_  ; // Current event
    emjet:: Jet    jet_    ; // Current jet
    emjet:: Track  track_  ; // Current track
    emjet:: Vertex vertex_ ; // Current vertex
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EmJetAnalyzer::EmJetAnalyzer(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed

}


EmJetAnalyzer::~EmJetAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
EmJetAnalyzer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
  return true;
}

// ------------ method called once each job just before starting event loop  ------------
void
EmJetAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
EmJetAnalyzer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
EmJetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
EmJetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
EmJetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
EmJetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EmJetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(EmJetAnalyzer);
