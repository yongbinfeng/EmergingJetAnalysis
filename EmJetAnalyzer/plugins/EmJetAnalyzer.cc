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
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometrySurface/interface/Line.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "RecoVertex/AdaptiveVertexFinder/interface/AdaptiveVertexReconstructor.h"
#include "RecoVertex/TrimmedKalmanVertexFinder/interface/KalmanTrimmedVertexFinder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/JetReco/interface/GenJet.h"

// track association
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"

#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"

#include "EmergingJetAnalysis/EmJetAnalyzer/interface/OutputTree.h"
#include "EmergingJetAnalysis/EmJetAnalyzer/interface/EmJetEvent.h"

#ifndef OUTPUT
#define OUTPUT(x) std::cout<<#x << ": " << x << std::endl
#endif

//
// class declaration
//

using namespace emjet;

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
    bool isData_;

    edm::Service<TFileService> fs;
    edm::EDGetTokenT< reco::PFJetCollection > jetCollectionToken_;

    edm::ParameterSet         m_trackParameterSet;
    TrackDetectorAssociator   m_trackAssociator;
    TrackAssociatorParameters m_trackParameters;

    emjet::OutputTree otree_ ; // OutputTree object
    TTree* t_tree;

    emjet:: Event  event_  ; // Current event
    emjet:: Jet    jet_    ; // Current jet
    emjet:: Track  track_  ; // Current track
    emjet:: Vertex vertex_ ; // Current vertex
    int jet_index_    ; // Current jet index
    int track_index_  ; // Current track index
    int vertex_index_ ; // Current vertex index

    // Retrieve once per event
    // Intermediate objects used for calculations
    edm::Handle<reco::VertexCollection> primary_verticesH_;

    edm::ESHandle<TransientTrackBuilder> transienttrackbuilderH_;
    std::vector<reco::TransientTrack> generalTracks_;
    edm::Handle<reco::GenParticleCollection> genParticlesH_;
    const reco::BeamSpot* theBeamSpot_;
    reco::VertexCollection selectedSecondaryVertices_;
    std::vector<TransientVertex> avrVertices_;
    edm::Handle<reco::GenJetCollection> genJets_;
    std::vector<GlobalPoint> dt_points_;
    std::vector<GlobalPoint> csc_points_;
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
EmJetAnalyzer::EmJetAnalyzer(const edm::ParameterSet& iConfig):
  // event_  (new Event  ()),
  // jet_    (new Jet    ()),
  // track_  (new Track  ()),
  // vertex_ (new Vertex ())
  event_  (),
  jet_    (),
  track_  (),
  vertex_ ()
{
  OUTPUT(jet_.track_vector.size());
  // OUTPUT(jet_.vertex_vector.size());
  OUTPUT(event_.jet_vector.size());
  // Config-independent initialization
  {
    t_tree           = fs->make<TTree>("emJetTree","emJetTree");
    otree_.Branch(t_tree);
  }

  // Config-dependent initialization
  {
    isData_ = iConfig.getUntrackedParameter<bool>("isData");
    m_trackParameterSet = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
    if (isData_) {
      std::cout << "running on data" << std::endl;
    } else {
      std::cout << "running on MC"   << std::endl;
    }

    edm::ConsumesCollector iC = consumesCollector();
    m_trackParameters.loadParameters( m_trackParameterSet, iC );
    m_trackAssociator.useDefaultPropagator();

    jetCollectionToken_ = consumes< reco::PFJetCollection > (iConfig.getParameter<edm::InputTag>("srcJets"));
  }
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
  // Reset output tree to default values
  otree_.Init();
  OUTPUT(jet_.track_vector.size());
  // OUTPUT(jet_.vertex_vector.size());
  OUTPUT(event_.jet_vector.size());
  Event a;
  event_ = a;
  // Reset Event variables
  std::cout<<"1" << std::endl;
  vertex_.Init();
  std::cout<<"2" << std::endl;
  jet_.Init();
  std::cout<<"3" << std::endl;
  event_.Init();

  event_.run   = iEvent.id().run();
  event_.event = iEvent.id().event();
  event_.lumi  = iEvent.id().luminosityBlock();
  event_.bx    = iEvent.bunchCrossing();

  // Calculate Vertex info :EVENTLEVEL:
  {
    iEvent.getByLabel("offlinePrimaryVerticesWithBS", primary_verticesH_);
    const reco::Vertex& primary_vertex = primary_verticesH_->at(0);
    if (!isData_) { // :MCONLY: Add true number of interactions
      edm::Handle<std::vector<PileupSummaryInfo> > PileupInfo;
      iEvent.getByLabel("addPileupInfo", PileupInfo);
      for (auto const& puInfo : *PileupInfo) {
        int bx = puInfo.getBunchCrossing();
        if (bx == 0) {
          event_.nTrueInt = puInfo.getTrueNumInteractions();
        }
      }
    }
    for (auto ipv = primary_verticesH_->begin(); ipv != primary_verticesH_->end(); ++ipv) {
      event_.nVtx ++;
      if ( (ipv->isFake()) || (ipv->ndof() <= 4.) || (ipv->position().Rho() > 2.0) || (fabs(ipv->position().Z()) > 24.0) ) continue; // :CUT: Primary vertex cut for counting
      event_.nGoodVtx++;
    }
  }

  // Calculate MET :EVENTLEVEL:
  {
    edm::Handle<reco::PFMETCollection> pfmet;
    iEvent.getByLabel("pfMet",pfmet);
    event_.met_pt = pfmet->begin()->pt();
    event_.met_phi = pfmet->begin()->phi();
  }

  // Write current Event to OutputTree
  WriteEventToOutput(event_, &otree_);
  // Write OutputTree to TTree
  t_tree->Fill();

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
