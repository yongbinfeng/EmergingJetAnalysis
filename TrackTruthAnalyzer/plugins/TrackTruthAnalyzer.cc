// -*- C++ -*-
//
// Package:    EmergingJetAnalysis/TrackTruthAnalyzer
// Class:      TrackTruthAnalyzer
//
/**\class TrackTruthAnalyzer TrackTruthAnalyzer.cc EmergingJetAnalysis/TrackTruthAnalyzer/plugins/TrackTruthAnalyzer.cc

   Description: Analyzer to evaluate tracking/vertexing efficiency, etc for Emerging Jet analysis

   Implementation:
   Heavily borrowed from SimTracker/TrackAssociation/test/testTrackAssociator.cc
*/
//
// Original Author:  Young Ho Shin
//         Created:  Mon, 07 Mar 2016 19:08:10 GMT
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

#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

//
// class declaration
//

class TrackTruthAnalyzer : public edm::EDFilter {
public:
  explicit TrackTruthAnalyzer(const edm::ParameterSet&);
  ~TrackTruthAnalyzer();

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
  edm::InputTag tracksTag, tpTag;
  edm::EDGetTokenT<reco::RecoToSimCollection> recSimToken;
  edm::EDGetTokenT<reco::SimToRecoCollection> simRecToken;
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
TrackTruthAnalyzer::TrackTruthAnalyzer(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  tracksTag = iConfig.getParameter< edm::InputTag >("tracksTag");
  tpTag = iConfig.getParameter< edm::InputTag >("tpTag");
  recSimToken = consumes<reco::RecoToSimCollection>( iConfig.getParameter< edm::InputTag >("recSim") );
  simRecToken = consumes<reco::SimToRecoCollection>( iConfig.getParameter< edm::InputTag >("simRec") );
}


TrackTruthAnalyzer::~TrackTruthAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TrackTruthAnalyzer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using reco::Track;

  Handle<View<Track> > trackCollectionH;
  iEvent.getByLabel(tracksTag,trackCollectionH);
  const View<Track>&  tC = *(trackCollectionH.product());

  Handle<TrackingParticleCollection> TPCollectionH;
  iEvent.getByLabel(tpTag,TPCollectionH);

  Handle<reco::RecoToSimCollection > recToSimCollection;
  iEvent.getByToken(recSimToken, recToSimCollection);
  Handle<reco::SimToRecoCollection > simToRecCollection;
  iEvent.getByToken(simRecToken, simToRecCollection);

  // for(View<Track>::size_type i=0; i<trackCollection->size(); ++i){
  //   RefToBase<Track> track(trackCollection, i);
  //   std::vector<std::pair<TrackingParticleRef, double> > tp;
  //   if(recSimColl.find(track) != recSimColl.end()){
  //     tp = recSimColl[track];
  //     if (tp.size()!=0) {
  //       TrackingParticleRef tpr = tp.begin()->first;
  //       double associationQuality = tp.begin()->second;
  //         }
  //   }
  // }
  std::cout << "Number of tracking particles" << TPCollectionH->size() << std::endl;
  for ( unsigned int itp = 0; itp != TPCollectionH->size(); itp++ ) {
    std::cout << "Running on tp index: " << itp << std::endl;
    int nMatch = 0;
    const TrackingParticleRef tpRef(TPCollectionH, itp);
    if ( simToRecCollection->find(tpRef) != simToRecCollection->end() ) {
      const auto vec_rt_quality = (*simToRecCollection)[tpRef];
      std::cout << "Number of matched recoTracks: " << vec_rt_quality.size() << std::endl;
      for (const auto rt_quality : vec_rt_quality) {
        auto rt = rt_quality.first;
        auto quality = rt_quality.second;
        std::cout << "    SIM-RECO association quality: " << quality << std::endl;
        nMatch++;
      }
    }
  }

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

// ------------ method called once each job just before starting iEvent loop  ------------
void
TrackTruthAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the iEvent loop  ------------
void
TrackTruthAnalyzer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
  void
  TrackTruthAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void
  TrackTruthAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void
  TrackTruthAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void
  TrackTruthAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackTruthAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TrackTruthAnalyzer);
