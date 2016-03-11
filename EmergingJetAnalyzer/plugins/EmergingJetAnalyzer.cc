// -*- C++ -*-
//
// Package:    MyAnalysis/EmergingJetAnalyzer
// Class:      EmergingJetAnalyzer
//
/**\class EmergingJetAnalyzer EmergingJetAnalyzer.cc MyAnalysis/EmergingJetAnalyzer/plugins/EmergingJetAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Ted Ritchie Kolberg
//         Created:  Tue, 27 Jan 2015 12:30:51 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

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

#include "EmergingJetAnalysis/EmergingJetAnalyzer/interface/OutputTree.h"
//
// class declaration
//

class EmergingJetAnalyzer : public edm::EDAnalyzer {
  public:
    explicit EmergingJetAnalyzer(const edm::ParameterSet&);
    ~EmergingJetAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ---------- helper functions ---------------------------
    // Take a single PFJet and add to output tree
    void fillSingleJet(const reco::PFJet&);
    // Select all secondary vertices passing certain criteria, from a given vertexCollection
    reco::VertexCollection selectSecondaryVertices (edm::Handle<reco::VertexCollection>) const;
    // Calculate sum of pt-squares of all tracks, for a given vertex
    double calculatePt2Sum (const reco::Vertex&) const;

    TrackDetectorAssociator   m_trackAssociator;
    TrackAssociatorParameters m_trackParameters;
    edm::ParameterSet         m_trackParameterSet;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------
    edm::Service<TFileService> fs;
    edm::EDGetTokenT< reco::PFJetCollection > jetCollectionToken_;

    emjet::OutputTree otree;

    // Retrieve once per event
    // Intermediate objects used for calculations
    edm::ESHandle<TransientTrackBuilder> transienttrackbuilder_;
    edm::Handle<reco::VertexCollection> primary_vertices_;
    std::vector<reco::TransientTrack> generalTracks_;
    edm::Handle<reco::GenParticleCollection> genParticlesH_;
    const reco::BeamSpot* theBeamSpot_;
    reco::VertexCollection selectedSecondaryVertices_;
    edm::Handle<reco::GenJetCollection> genJets_;
    std::vector<GlobalPoint> dt_points_;
    std::vector<GlobalPoint> csc_points_;


    TH1F * h_dr_jet_track;

    TH1F * h_genHT;

    TH2F * h_ptWtPixHits_vs_chargeFraction;
    TH2F * h_ptWtPixHits_vs_emFraction;
    TH1F * h_ptWtD0;
    TH1F * h_ptWtD0_trackFrac05;
    TH2F * h_jetEtaPhi;
    TH1F * h_nVtxInEvent;
    TH1F * h_rVtxGenTracks;
    TH2F * h_jetPtVsIndex;
    TH1F * h_eventSumLogIp;
    TH1F * h_jetSumLogIp;
    TH1F * h_eventHT;
    TH2F * h_rVsLogSig;
    TH2F * h_rVsZ;
    TH2F * h_xVsY;
    TH1F * h_jetTrackSumLogIp;
    TH1F * h_vertexMass;
    TH1F * h_vertexPt;
    TH1F * h_trackFrac;
    TH1F * h_trackLogD0;
    TH1F * h_jet_cef;
    TH1F * h_jet_nef;
    TH1F * h_jet_chf;
    TH1F * h_jet_nhf;
    TH1F * h_jet_phf;
    TH1F * h_jet_cef_lt5prtk;
    TH1F * h_jet_nef_lt5prtk;
    TH1F * h_jet_chf_lt5prtk;
    TH1F * h_jet_nhf_lt5prtk;
    TH1F * h_jet_phf_lt5prtk;
    TH1F * h_jet_cef_trackFrac05;
    TH1F * h_jet_nef_trackFrac05;
    TH1F * h_jet_chf_trackFrac05;
    TH1F * h_jet_nhf_trackFrac05;

    TH2F * h_log_rz_vertex;

    TH1F * h_pTxIPxSig;

    TH1F * h_nPV;
    TH1F * h_nSV;
    TH1F * h_nSV_disp;

    TH1F * h_nSV_perjet;

    TH1F * h_rVtx;
    TH1F * h_rSV;
    TH1F * h_chi2Vtx;
    TH1F * h_nTracksVtx;
    TH1F * h_massVtx;
    TH1F * h_ptVtx;

    TH1F * h_jet_phf_trackFrac05;

    TH1F * h_jetIpTracks;
    TH1F * h_jetMatchedVertices;
    TH1F * h_vertexDrJet;
    TH1F * h_jetDispMuons;
    TH1F * h_muonDrJet;

    TH1F * h_eventIpTracks;
    TH1F * h_eventIpTracks0_1;
    TH1F * h_eventIpTracks1_0;
    TH1F * h_eventIpTracks10_;
    TH1F * h_eventVertices;
    TH1F * h_eventDispMuons;

    TH1F * h_jetNeutralFraction;

    TH1F * h_sumEt;

    TH1F * h_logTrackDz;
    TH1F * h_logTrackDxy;
    TH1F * h_logTrackDz_genMatched;
    TH1F * h_logTrackDxy_genMatched;


    TH1F * h_medianLogIpSig;
    TH1F * h_medianLogIpSig_2highest;
    TH1F * h_medianIp;

    TH1F * h_Lxy_perjet;

    TH1F * h_pTdispTracks;

    TH1F * h_deltaEtaJetTrack;
    TH1F * h_deltaPhiJetTrack;

    TH1F * h_jetPromptTracks;
    TH1F * h_jetPromptTracks_2lowest;
    TH1F * h_jetDispTracks;

    TH1F * h_rDecayGen;

    TH1F * h_tagsPerEvent;
    TH1F * h_median_rPCA;
    TH1F * h_rPCA;

    TH1F * h_MET;
    TH1F * h_dPhiMetJet;

    TH2F * h_jetPCAdist;

    TH1F * h_DTsegments;
    TH1F * h_CSCsegments;
    TH1F * h_RPChits;

    TH1F * h_dr_disp_jet;

    TH1F * h_pca_dr_jet;
    TH1F * h_pca_r;

    TH1F * h_sum_ipSig;

    TH1F * h_nMissingHitsPerJet;

    TH1F * h_dsz_perJet;

    TH1F * h_nDarkPions_genJet;

    TH1F * h_nJets_gt1piD;
    TH1F * h_nMatchedGen;

    TH2F * h_dz_vs_dxy;
    TH2F * h_dz_vs_dxy_genMatched;

    TH1F * h_3Dip;
    TH1F * h_3Dip_genMatched;
    TH1F * h_3DipSig;
    TH1F * h_3DipSig_genMatched;

    TH1F * h_trackPtInJet;
    TH1F * h_trackPtInJet_genMatched;

    TH1F * h_medianSVradius;

    bool isData_;

    TTree * t_tree;

    float t_met;

    float jet0_pt;
    float jet0_eta;
    int jet0_promptTracks;
    int jet0_dispTracks;
    int jet0_nSV;
    float jet0_cef;
    float jet0_nef;
    float jet0_chf;
    float jet0_nhf;
    float jet0_phf;
    float jet0_medianLogIpSig;
    int   jet0_missHits;
    int   jet0_muonHits;

    float jet1_pt;
    float jet1_eta;
    int jet1_promptTracks;
    int jet1_dispTracks;
    int jet1_nSV;
    float jet1_cef;
    float jet1_nef;
    float jet1_chf;
    float jet1_nhf;
    float jet1_phf;
    float jet1_medianLogIpSig;
    int   jet1_missHits;
    int   jet1_muonHits;

    float jet2_pt;
    float jet2_eta;
    int jet2_promptTracks;
    int jet2_dispTracks;
    int jet2_nSV;
    float jet2_cef;
    float jet2_nef;
    float jet2_chf;
    float jet2_nhf;
    float jet2_phf;
    float jet2_medianLogIpSig;
    int   jet2_missHits;
    int   jet2_muonHits;

    float jet3_pt;
    float jet3_eta;
    int jet3_promptTracks;
    int jet3_dispTracks;
    int jet3_nSV;
    float jet3_cef;
    float jet3_nef;
    float jet3_chf;
    float jet3_nhf;
    float jet3_phf;
    float jet3_medianLogIpSig;
    int   jet3_missHits;
    int   jet3_muonHits;
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
EmergingJetAnalyzer::EmergingJetAnalyzer(const edm::ParameterSet& iConfig) :
  m_trackParameterSet(iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters")),
  isData_(iConfig.getUntrackedParameter<bool>("isData"))
{
  //now do what ever initialization is needed

  if (isData_) {
    std::cout << "running on data" << std::endl;
  } else {
    std::cout << "running on MC"   << std::endl;
  }

  edm::ConsumesCollector iC = consumesCollector();
  m_trackParameters.loadParameters( m_trackParameterSet, iC );
  m_trackAssociator.useDefaultPropagator();

  jetCollectionToken_ = consumes< reco::PFJetCollection > (iConfig.getParameter<edm::InputTag>("srcJets"));

  //    h_ptWtPixHits_vs_chargeFraction = fs->make<TH2F>( "h_ptWtPixHits_vs_chargeFraction"  , ";p_{T} wt. PIX hits;charge fraction", 10,  0., 4., 10, 0., 1. );
  //    h_ptWtPixHits_vs_emFraction = fs->make<TH2F>( "h_ptWtPixHits_vs_emFraction"  , ";p_{T} wt. PIX hits;EM fraction", 10,  0., 4., 10, 0., 1. );

  //    h_ptWtD0 = fs->make<TH1F>( "h_ptWtD0",";p_{T} weighted d0 [cm];" , 25, -3., 3.);
  //    h_ptWtD0_trackFrac05 = fs->make<TH1F>( "h_ptWtD0_trackFrac05",";p_{T} weighted d0 [cm];" , 25, -3., 3.);
  //    h_jetEtaPhi = fs->make<TH2F>( "h_jetEtaPhi" , ";jet #eta;jet #phi;" , 25, -2.5, 2.5, 25, -TMath::Pi(), TMath::Pi() );

  h_genHT      = fs->make<TH1F> ("h_genHT",";gen HT [GeV];",40,0.,4000.);

  h_medianSVradius = fs->make<TH1F> ("h_medianSVradius",";median SV radius;",100,0.,100.);

  h_dz_vs_dxy  = fs->make<TH2F> ("h_dz_vs_dxy", ";log(dxy);log(dz)",40,-7.,5.,40,-7.,5.);
  h_dz_vs_dxy_genMatched = fs->make<TH2F> ("h_dz_vs_dxy_genMatched", ";log(dxy);log(dz)",40,-7.,5.,40,-7.,5.);

  h_dsz_perJet = fs->make<TH1F> ("h_dsz_perJet",";dsz per track;",40,0.,100.);

  h_DTsegments = fs->make<TH1F> ("h_DTsegments",";DT segments per jet;",40,-0.5,199.5);
  h_CSCsegments = fs->make<TH1F> ("h_CSCsegments",";CSC segments per jet;",40,-0.5,199.5);
  h_RPChits = fs->make<TH1F> ("h_RPChits",";RPC hits per jet;",40,-0.5,199.5);
  h_rPCA = fs->make<TH1F> ("h_rPCA",";r_{PCA} [cm];",40,0.,100.);

  h_trackPtInJet         = fs->make<TH1F> ("h_trackPtInJet",";track p_{T} [GeV];",40,0.,100.);
  h_trackPtInJet_genMatched = fs->make<TH1F> ("h_trackPtInJet_genMatched",";track p_{T} [GeV];",40,0.,100.);

  h_median_rPCA=fs->make<TH1F> ("h_median_rPCA",";median r_{PCA} (2 highest);",40,0.,20.);
  h_pTxIPxSig = fs->make<TH1F> ("h_pTxIPxSig",";pT #times IP #times IPsig;", 40.,-4.,10.);
  h_log_rz_vertex = fs->make<TH2F> ("h_log_rz_vertex",";log(z);log(r)", 100, -5., 5., 100, -5.,5.);

  h_deltaEtaJetTrack = fs->make<TH1F> ("h_deltaEtaJetTrack", ";d#eta(jet, trk);", 40, 0.,1.5);
  h_deltaPhiJetTrack = fs->make<TH1F> ("h_deltaPhiJetTrack", ";d#phi(jet, trk);", 40, 0.,1.5);

  h_pTdispTracks = fs->make<TH1F> ("h_pTdispTracks", ";p_{T} of displaced tracks;",40,0.,200.);
  h_dr_jet_track = fs->make<TH1F> ("h_dr_jet_track", ";dR(jet, track);",20,0.,TMath::Pi()/2);

  h_nDarkPions_genJet = fs->make<TH1F> ("h_nDarkPions_genJet", ";N(#pi_{D}) per genJet;",20,-0.5,19.5);

  h_pca_dr_jet   = fs->make<TH1F> ("h_pca_dr_jet",   ";DR of PCA to jet;",20.,0.,1.5);
  h_pca_r        = fs->make<TH1F> ("h_pca_r",        ";r of PCA to jet [cm];",20,0.,150.);

  h_nMatchedGen = fs->make<TH1F>( "h_nMatchedGen",  ";N_{GEN matched};", 5, -0.5,4.5);
  h_nVtxInEvent = fs->make<TH1F>( "h_nVtxInEvent",  ";nVtx;",101,-0.5,100.5);
  h_rVtx      = fs->make<TH1F>("h_rVtx",";vertex r [cm];",101,-0.5,100.5);
  h_rSV      = fs->make<TH1F>("h_rSV",";vertex r [cm];",101,-0.5,100.5);
  h_chi2Vtx   = fs->make<TH1F>("h_chi2Vtx",";vertex #chi^{2};",16,-0.5,16);
  h_nTracksVtx= fs->make<TH1F>("h_nTracksVtx",";N tracks/vertex;",10,-0.5,9.5);
  h_massVtx   = fs->make<TH1F>("h_massVtx",";mass [GeV];",20,0.,10.);
  h_ptVtx     = fs->make<TH1F>("h_ptVtx",";pT [GeV];",50,0.,100.);

  h_nJets_gt1piD = fs->make<TH1F>("h_nJets_gt1piD",";N(jets) with > 1 #pi_{D};",5,-0.5,4.5);

  h_nPV       = fs->make<TH1F>("h_nPV",";N_{PV};",50,-0.5,49.5);
  h_nSV       = fs->make<TH1F>("h_nSV",";N_{SV};",50,-0.5,49.5);
  h_nSV_disp       = fs->make<TH1F>("h_nSV_disp",";N_{SV};",50,-0.5,49.5);

  h_3Dip = fs->make<TH1F>("h_3Dip",";log(3D IP);",40,-8.,3.);
  h_3Dip_genMatched = fs->make<TH1F>("h_3Dip_genMatched",";log(3D IP);",40,-8.,3.);
  h_3DipSig = fs->make<TH1F>("h_3DipSig",";3D IP significance;",40,0.,40.);
  h_3DipSig_genMatched = fs->make<TH1F>("h_3DipSig_genMatched",";3D IP significance;",40,0.,40.);

  h_MET       = fs->make<TH1F>("h_MET",";MET;",50,0.,400.);
  h_dPhiMetJet= fs->make<TH1F>("h_dPhiMetJet",";#Delta#phi(MET, jet);",50,-TMath::Pi()/2.,TMath::Pi()/2.);

  h_jetPromptTracks = fs->make<TH1F>("h_jetPromptTracks", ";N prompt tracks;", 40, -0.5, 39.5);
  h_jetPromptTracks_2lowest = fs->make<TH1F>("h_jetPromptTracks_2lowest", ";N prompt tracks;", 40, -0.5, 39.5);
  h_jetDispTracks = fs->make<TH1F>("h_jetDispTracks", ";N prompt tracks;", 40, -0.5, 39.5);

  //    h_rVtxGenTracks = fs->make<TH1F>( "h_rVtxGenTracks", ";vertex r [cm];",300.,0.,150.);
  h_tagsPerEvent = fs->make<TH1F>("h_tagsPerEvent",";N_{tags};",10,-0.5,9.5);

  h_jetPtVsIndex = fs->make<TH2F>( "h_jetPtVsIndex", ";p_{T} [GeV];index", 100,0.,1000.,4,-0.5,3.5);
  h_jetIpTracks  = fs->make<TH1F>("h_jetIpTracks",";N disp trks;",20,-0.5,19.5);
  h_jetMatchedVertices  = fs->make<TH1F>("h_jetMatchedVertices",";N disp vtx;",20,-0.5,19.5);
  h_jetDispMuons  = fs->make<TH1F>("h_jetDispMuons",";N disp STA #mu;",20,-0.5,19.5);
  h_vertexDrJet   = fs->make<TH1F>("h_vertexDrJet",";#DeltaR(track, jet);",20.,0.,1.5);
  h_muonDrJet   = fs->make<TH1F>("h_muonDrJet",";#DeltaR(track, jet);",20.,0.,1.5);
  h_eventIpTracks  = fs->make<TH1F>("h_eventIpTracks",";N disp trks;",50,-0.5,99.5);
  h_eventIpTracks0_1  = fs->make<TH1F>("h_eventIpTracks0_1",";N disp trks;",20,-0.5,19.5);
  h_eventIpTracks1_0  = fs->make<TH1F>("h_eventIpTracks1_0",";N disp trks;",20,-0.5,19.5);
  h_eventIpTracks10_  = fs->make<TH1F>("h_eventIpTracks10_",";N disp trks;",20,-0.5,19.5);
  h_eventVertices  = fs->make<TH1F>("h_eventVertices",";N disp vtx;",20,-0.5,19.5);
  h_eventDispMuons  = fs->make<TH1F>("h_eventDispMuons",";N disp STA #mu;",20,-0.5,19.5);

  h_dr_disp_jet = fs->make<TH1F>("h_dr_disp_jet",";#Delta R (genParticle, nearest jet);",20,0.,1.5);

  h_eventSumLogIp = fs->make<TH1F>( "h_eventSumLogIp", ";#Sigma(ln[ipSig]);",100,3.,11.);
  //    h_jetSumLogIp = fs->make<TH1F>( "h_jetSumLogIp", ";#Sigma(ln[ipSig]);",100,0.,100.);

  //    h_jetTrackSumLogIp = fs->make<TH1F>("h_jetTrackSumLogIp",";#Sigma(ln[ipSig]);",100,0.,100.);

  //    h_eventHT = fs->make<TH1F>("h_eventHT", ";H_{T} [GeV];",100.,0.,3000.);

  //    h_rVsLogSig = fs->make<TH2F>("h_rVsLogSig", ";r [cm];significance", 150,0.,150,20,-5.,15.);

  //    h_rVsZ = fs->make<TH2F>("h_rVsZ", ";z [cm];r [cm]",150,0.,150.,100,0.,100.);
  //    h_xVsY = fs->make<TH2F>("h_xVsY", ";x [cm];y [cm]",200,-100.,100.,200,-100.,100.);

  //    h_vertexMass = fs->make<TH1F>("h_vertexMass",";M [GeV];",20,0.,20.);
  //    h_vertexPt   = fs->make<TH1F>("h_vertexPt",  ";pT [GeV];",20,0.,20.);

  h_trackFrac  = fs->make<TH1F>("h_trackFrac", ";track fraction;",20,0.,1.);
  //    h_trackLogD0 = fs->make<TH1F>("h_trackLogD0",";log(|d0|);",50,-3.,3.);

  h_jet_cef    = fs->make<TH1F>("h_jet_cef",";charged EM fraction;",20,0.,1.);
  h_jet_nef    = fs->make<TH1F>("h_jet_nef",";neutral EM fraction;",20,0.,1.);
  h_jet_chf    = fs->make<TH1F>("h_jet_chf",";charged had fraction;",20,0.,1.);
  h_jet_nhf    = fs->make<TH1F>("h_jet_nhf",";neutral had fraction;",20,0.,1.);
  h_jet_phf    = fs->make<TH1F>("h_jet_phf",";photon fraction;",20,0.,1.);
  h_jet_cef_lt5prtk     = fs->make<TH1F>("h_jet_cef_lt5prtk",";charged EM fraction;",20,0.,1.);
  h_jet_nef_lt5prtk     = fs->make<TH1F>("h_jet_nef_lt5prtk",";neutral EM fraction;",20,0.,1.);
  h_jet_chf_lt5prtk     = fs->make<TH1F>("h_jet_chf_lt5prtk",";charged had fraction;",20,0.,1.);
  h_jet_nhf_lt5prtk     = fs->make<TH1F>("h_jet_nhf_lt5prtk",";neutral had fraction;",20,0.,1.);
  h_jet_phf_lt5prtk     = fs->make<TH1F>("h_jet_phf_lt5prtk",";photon fraction;",20,0.,1.);
  h_jetNeutralFraction =  fs->make<TH1F>("h_jetNeutralFraction",";neutral fraction;",20,0.,1.);
  h_Lxy_perjet= fs->make<TH1F>("h_Lxy_perjet",";log(#Sigma L_{xy});",50,-3.,9.);

  h_nSV_perjet = fs->make<TH1F>("h_nSV_perjet",";N_{SV} per jet;",20,-0.5,19.5);

  h_sumEt      = fs->make<TH1F>("h_sumEt",";#Sigma(E_{T}) [GeV];",130,0.,13000.);

  h_logTrackDxy   = fs->make<TH1F>("h_logTrackDxy",";log(dxy);",100.,-7.,5.);
  h_logTrackDz    = fs->make<TH1F>("h_logTrackDz", ";log(dz);", 100.,-7.,5.);
  h_logTrackDxy_genMatched   = fs->make<TH1F>("h_logTrackDxy_genMatched",";log(dxy);",100.,-7.,5.);
  h_logTrackDz_genMatched    = fs->make<TH1F>("h_logTrackDz_genMatched", ";log(dz);", 100.,-7.,5.);

  h_medianLogIpSig = fs->make<TH1F>("h_medianLogIpSig",";median log(ipSig);",20.,-1.,8.);
  h_medianLogIpSig_2highest = fs->make<TH1F>("h_medianLogIpSig_2highest",";median log(ipSig);",20.,-1.,8.);
  h_medianIp       = fs->make<TH1F>("h_medianIp",      ";median IP;",        100.,-3.,3.);

  h_rDecayGen      = fs->make<TH1F>("h_rDecayGen",   ";gen decay r [cm];" , 100, 0., 500.);

  h_jetPCAdist     = fs->make<TH2F>("h_jetPCAdist",  ";PCA to jet axis [cm];", 50, 0.,0.5,50,0.,100.);

  h_sum_ipSig      = fs->make<TH1F>("h_sum_ipSig",   ";sum(IP sig);", 20, 0., 250.);

  h_nMissingHitsPerJet = fs->make<TH1F>("h_nMissingHitsPerJet",  ";missing hits;", 20, 0., 40.);


  //    h_jet_cef_trackFrac05     = fs->make<TH1F>("h_jet_cef_trackFrac05",";charged EM fraction;",20,0.,1.);
  //    h_jet_nef_trackFrac05     = fs->make<TH1F>("h_jet_nef_trackFrac05",";neutral EM fraction;",20,0.,1.);
  //    h_jet_chf_trackFrac05     = fs->make<TH1F>("h_jet_chf_trackFrac05",";charged had fraction;",20,0.,1.);
  //    h_jet_nhf_trackFrac05     = fs->make<TH1F>("h_jet_nhf_trackFrac05",";neutral had fraction;",20,0.,1.);
  //    h_jet_phf_trackFrac05     = fs->make<TH1F>("h_jet_phf_trackFrac05",";photon fraction;",20,0.,1.);

  t_tree           = fs->make<TTree>("emergingJetsTree","emergingJetsTree");
  otree.Branch(t_tree);

  t_tree->Branch("jet0_pt", &jet0_pt);
  t_tree->Branch("jet0_eta", &jet0_eta);
  t_tree->Branch("jet0_promptTracks", &jet0_promptTracks);
  t_tree->Branch("jet0_dispTracks", &jet0_dispTracks);
  t_tree->Branch("jet0_nSV", &jet0_nSV);
  t_tree->Branch("jet0_cef", &jet0_cef);
  t_tree->Branch("jet0_nef", &jet0_nef);
  t_tree->Branch("jet0_chf", &jet0_chf);
  t_tree->Branch("jet0_nhf", &jet0_nhf);
  t_tree->Branch("jet0_phf", &jet0_phf);
  t_tree->Branch("jet0_medianLogIpSig", &jet0_medianLogIpSig);
  t_tree->Branch("jet0_missHits", &jet0_missHits);
  t_tree->Branch("jet0_muonHits", &jet0_muonHits);

  t_tree->Branch("jet1_pt", &jet1_pt);
  t_tree->Branch("jet1_eta", &jet1_eta);
  t_tree->Branch("jet1_promptTracks", &jet1_promptTracks);
  t_tree->Branch("jet1_dispTracks", &jet1_dispTracks);
  t_tree->Branch("jet1_nSV", &jet1_nSV);
  t_tree->Branch("jet1_cef", &jet1_cef);
  t_tree->Branch("jet1_nef", &jet1_nef);
  t_tree->Branch("jet1_chf", &jet1_chf);
  t_tree->Branch("jet1_nhf", &jet1_nhf);
  t_tree->Branch("jet1_phf", &jet1_phf);
  t_tree->Branch("jet1_medianLogIpSig", &jet1_medianLogIpSig);
  t_tree->Branch("jet1_missHits", &jet1_missHits);
  t_tree->Branch("jet1_muonHits", &jet1_muonHits);

  t_tree->Branch("jet2_pt", &jet2_pt);
  t_tree->Branch("jet2_eta", &jet2_eta);
  t_tree->Branch("jet2_promptTracks", &jet2_promptTracks);
  t_tree->Branch("jet2_dispTracks", &jet2_dispTracks);
  t_tree->Branch("jet2_nSV", &jet2_nSV);
  t_tree->Branch("jet2_cef", &jet2_cef);
  t_tree->Branch("jet2_nef", &jet2_nef);
  t_tree->Branch("jet2_chf", &jet2_chf);
  t_tree->Branch("jet2_nhf", &jet2_nhf);
  t_tree->Branch("jet2_phf", &jet2_phf);
  t_tree->Branch("jet2_medianLogIpSig", &jet2_medianLogIpSig);
  t_tree->Branch("jet2_missHits", &jet2_missHits);
  t_tree->Branch("jet2_muonHits", &jet2_muonHits);

  t_tree->Branch("jet3_pt", &jet3_pt);
  t_tree->Branch("jet3_eta", &jet3_eta);
  t_tree->Branch("jet3_promptTracks", &jet3_promptTracks);
  t_tree->Branch("jet3_dispTracks", &jet3_dispTracks);
  t_tree->Branch("jet3_nSV", &jet3_nSV);
  t_tree->Branch("jet3_cef", &jet3_cef);
  t_tree->Branch("jet3_nef", &jet3_nef);
  t_tree->Branch("jet3_chf", &jet3_chf);
  t_tree->Branch("jet3_nhf", &jet3_nhf);
  t_tree->Branch("jet3_phf", &jet3_phf);
  t_tree->Branch("jet3_medianLogIpSig", &jet3_medianLogIpSig);
  t_tree->Branch("jet3_missHits", &jet3_missHits);
  t_tree->Branch("jet3_muonHits", &jet3_muonHits);
}


EmergingJetAnalyzer::~EmergingJetAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
  void
EmergingJetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  // Reset output tree to default values
  otree.Init();

  //   float ipCut = 0.05;

  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  // edm::ESHandle<TransientTrackBuilder> transienttrackbuilder_;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transienttrackbuilder_);

  edm::Handle<reco::PFJetCollection> pfjetH;
  // iEvent.getByLabel("ak4PFJetsCHS", pfjetH);
  iEvent.getByToken(jetCollectionToken_, pfjetH);
  auto selectedJets = *pfjetH;

  Handle<reco::BeamSpot> theBeamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", theBeamSpotHandle);
  theBeamSpot_ = theBeamSpotHandle.product();

  //   edm::Handle<reco::CaloJetCollection> caloJetH;
  //   iEvent.getByLabel("ak5CaloJets", caloJetH);

  edm::Handle<reco::TrackCollection> genTrackH;
  iEvent.getByLabel("generalTracks", genTrackH);
  generalTracks_ = transienttrackbuilder_->build(genTrackH);

  edm::Handle<reco::TrackCollection> sdmH;
  iEvent.getByLabel("displacedStandAloneMuons",sdmH);
  std::vector<reco::TransientTrack> standaloneDisplacedMuons;
  standaloneDisplacedMuons = transienttrackbuilder_->build(sdmH);

  edm::Handle<reco::PFMETCollection> pfmet;
  iEvent.getByLabel("pfMet",pfmet);

  iEvent.getByLabel("offlinePrimaryVerticesWithBS",primary_vertices_);
  const reco::Vertex& primary_vertex = primary_vertices_->at(0);

  edm::Handle<reco::VertexCollection> secondary_vertices;
  iEvent.getByLabel("inclusiveSecondaryVertices",secondary_vertices);

  ////////////////////////////////////////////////////////////
  // Reconstruct vertices from scratch
  ////////////////////////////////////////////////////////////
  if (1)
  {
    TStopwatch timer;
    timer.Start();
    std::cout << "Running KalmanTrimmedVertexFinder" << std::endl;
    KalmanTrimmedVertexFinder finder;
    finder.setPtCut(0.0);
    finder.setTrackCompatibilityCut(0.05);
    finder.setTrackCompatibilityToSV(0.01);
    finder.setVertexFitProbabilityCut(0.01);
    finder.setMaxNbOfVertices(0);
    vector<TransientVertex> vertices = finder.vertices ( generalTracks_ );
    std::cout << "Number of KTVF vertices: " << vertices.size() << std::endl;
    int nSV = 0;
    for (TransientVertex vertex: vertices) {
      auto vtx = reco::Vertex(vertex);
      double Lxy = 0;
      float dx = primary_vertex.position().x() - vtx.position().x();
      float dy = primary_vertex.position().y() - vtx.position().y();
      Lxy = TMath::Sqrt( dx*dx + dy*dy );
      if (Lxy>1) nSV++;
    }
    std::cout << "Number of displaced KTVF vertices: " << nSV << std::endl;
    std::cout << "Time elapsed:" << timer.RealTime() << std::endl;

    timer.Start();
    std::cout << "Running AdaptiveVertexReconstructor" << std::endl;
    AdaptiveVertexReconstructor avr (2.0, 6.0, 0.5, true );
    std::vector<TransientVertex> theVertices = avr.vertices(generalTracks_);
    std::cout << "Number of AVR vertices: " << theVertices.size() << std::endl;
    nSV = 0;
    for (TransientVertex vertex: theVertices) {
      auto vtx = reco::Vertex(vertex);
      double Lxy = 0;
      float dx = primary_vertex.position().x() - vtx.position().x();
      float dy = primary_vertex.position().y() - vtx.position().y();
      Lxy = TMath::Sqrt( dx*dx + dy*dy );
      if (Lxy>1) nSV++;
    }
    std::cout << "Number of displaced AVR vertices: " << nSV << std::endl;
    std::cout << "Time elapsed:" << timer.RealTime() << std::endl;

    nSV = 0;
    for (size_t ivx = 0; ivx < secondary_vertices->size(); ++ivx) {
      auto vtx = secondary_vertices->at(ivx);
      double Lxy = 0;
      float dx = primary_vertex.position().x() - vtx.position().x();
      float dy = primary_vertex.position().y() - vtx.position().y();
      Lxy = TMath::Sqrt( dx*dx + dy*dy );
      if (Lxy>1) nSV++;
    }
    std::cout << "Number of displaced vertices in inclusiveSecondaryVertices: " << nSV << std::endl;
  }


  // Fill event level GEN quantities
  if (!isData_) {
    edm::Handle<std::vector<reco::GenMET> > genMetH;
    iEvent.getByLabel("genMetTrue", genMetH);

    iEvent.getByLabel("genParticles", genParticlesH_);
    edm::Handle<std::vector<reco::GenJet> >  genJetH;
    iEvent.getByLabel("ak4GenJets",   genJetH);

  }

  otree.run   = iEvent.id().run();
  otree.event = iEvent.id().event();
  otree.lumi  = iEvent.id().luminosityBlock();
  otree.bx    = iEvent.bunchCrossing();

  // Do stuff with muon system
  {
    Handle<DTRecSegment4DCollection> segments4D;
    iEvent.getByLabel("dt4DSegments",segments4D);

    GlobalPoint segpos;
    dt_points_.clear();
    for (auto iseg = segments4D->begin(); iseg != segments4D->end(); ++iseg) {
      const GeomDet* geomDet = theTrackingGeometry->idToDet((*iseg).geographicalId());
      segpos = geomDet->toGlobal((*iseg).localPosition());
      dt_points_.push_back(segpos);
    }

    Handle<CSCSegmentCollection> csc_segments;
    iEvent.getByLabel("cscSegments",csc_segments);

    csc_points_.clear();
    for (auto icsc = csc_segments->begin(); icsc != csc_segments->end(); ++icsc) {
      const GeomDet* geomDet = theTrackingGeometry->idToDet((*icsc).geographicalId());
      segpos = geomDet->toGlobal((*icsc).localPosition());
      csc_points_.push_back(segpos);
    }

    Handle<RPCRecHitCollection> rpc_hits;
    iEvent.getByLabel("rpcRecHits",rpc_hits);
  }

  otree.met_pt = pfmet->begin()->pt();
  otree.met_phi = pfmet->begin()->phi();

  // Select secondary vertices
  selectedSecondaryVertices_ = selectSecondaryVertices(secondary_vertices);
  // Loop over secondary vertices
  for (auto vtx : selectedSecondaryVertices_) {
    float Lxy = 0.;
    float mass = 0.;
    float pt2sum = 0.;
    float dx = primary_vertex.position().x() - vtx.position().x();
    float dy = primary_vertex.position().y() - vtx.position().y();
    Lxy = TMath::Sqrt( dx*dx + dy*dy );
    mass = vtx.p4().mass();
    pt2sum = calculatePt2Sum(vtx);
    // bool hasRefittedTracks = vtx.hasRefittedTracks();
    // std::cout << "pt2sum: " << pt2sum << std::endl;
    // std::cout << "hasRefittedTracks: " << hasRefittedTracks << std::endl;
    otree.vertex_x      .push_back( vtx.x()                   );
    otree.vertex_y      .push_back( vtx.y()                   );
    otree.vertex_z      .push_back( vtx.z()                   );
    otree.vertex_xError .push_back( vtx.xError()              );
    otree.vertex_yError .push_back( vtx.yError()              );
    otree.vertex_zError .push_back( vtx.zError()              );
    otree.vertex_Lxy    .push_back( Lxy                       );
    otree.vertex_mass   .push_back( mass                      );
    otree.vertex_chi2   .push_back( vtx.chi2()                );
    otree.vertex_ndof   .push_back( vtx.ndof()                );
    otree.vertex_pt2sum .push_back( pt2sum                    );
  }

  ////////////////////////////////////////////////////////////
  // selectedJet Loop begin
  ////////////////////////////////////////////////////////////
  int ijet=0;
  for ( reco::PFJetCollection::const_iterator jet = selectedJets.begin(); jet != selectedJets.end(); ++jet ) {
    fillSingleJet(*jet);
    ++ijet;
  }
  ////////////////////////////////////////////////////////////
  // selectedJet Loop end
  ////////////////////////////////////////////////////////////

  t_tree->Fill();

  return;



}


// ------------ method called once each job just before starting event loop  ------------
  void
EmergingJetAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
  void
EmergingJetAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
/*
   void
   EmergingJetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void
   EmergingJetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void
   EmergingJetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void
   EmergingJetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EmergingJetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
EmergingJetAnalyzer::fillSingleJet(const reco::PFJet& jet) {
  // Shared objects
  const reco::Vertex& primary_vertex = primary_vertices_->at(0);
  TLorentzVector jetVector; jetVector.SetPtEtaPhiM(jet.pt(),jet.eta(),jet.phi(),0.);
  // std::cout << "jet.pt(): " << jet.pt() << std::endl;
  const float maxSigPromptTrack = 3.;
  const float minSigDispTrack = 3.;

  // Calculate nPromptTracks
  int nPromptTracks = 0;
  {
    reco::TrackRefVector trackRefs = jet.getTrackRefs();
    // Loop over tracks belonging to jet and calculate nPromptTracks
    for (reco::TrackRefVector::iterator ijt = trackRefs.begin(); ijt != trackRefs.end(); ++ijt) {
      reco::TransientTrack itk = transienttrackbuilder_->build(*ijt);
      if (itk.track().pt() < 1.) continue;
      auto d3d_ipv = IPTools::absoluteImpactParameter3D(itk, primary_vertex);
      if (d3d_ipv.second.significance() < maxSigPromptTrack) nPromptTracks++;
      //       std::cout << "Track with value, significance " << dxy_ipv.second.value() << "\t" << dxy_ipv.second.significance() << std::endl;
      //           if (dxy_ipv.second.value() > ipCut) continue;
    }
  }

  // Calculate nDispTracks, medianLogIpSig
  int nDispTracks   = 0;
  float medianLogIpSig = 0.;
  int nTags = 0;
  const float logTagCut    = 2.;
  // per track variables, reuse for efficiency
  std::vector<int>   vec_algo;
  std::vector<int>   vec_originalAlgo;
  std::vector<int>   vec_nHits;
  std::vector<int>   vec_nMissInnerHits;
  std::vector<float> vec_ipXY;
  std::vector<float> vec_ipZ;
  std::vector<float> vec_ipXYSig;
  int itrack = 0;
  {
    vec_algo           .clear();
    vec_originalAlgo   .clear();
    vec_nHits          .clear();
    vec_nMissInnerHits .clear();
    vec_ipXY           .clear();
    vec_ipZ            .clear();
    vec_ipXYSig        .clear();

    std::vector<float> ipVector;
    int ip0_1 = 0.;
    int ip1_0 = 0.;
    int ip10_ = 0.;
    float pTxIPxSig = 0.;
    std::vector<reco::TransientTrack> ipTracks;
    float trackFrac = 0.;

    math::XYZVector jetMomentum = jet.momentum();
    GlobalVector direction(jetMomentum.x(), jetMomentum.y(), jetMomentum.z());

    int misshits = 0;
    float tot_dsz = 0.;
    std::vector<float> r_pca;
    std::vector<float> logIpSig;
    float sum_Lxy = 0.;
    // std::cout << "generalTracks_.size(): " << generalTracks_.size() << std::endl;
    for (std::vector<reco::TransientTrack>::iterator itk = generalTracks_.begin(); itk != generalTracks_.end(); ++itk) {

      // Skip tracks with pt<1
      if (itk->track().pt() < 1.) continue;
      itrack++;

      auto dxy_ipv = IPTools::absoluteTransverseImpactParameter(*itk, primary_vertex);
      float ipXY = fabs(dxy_ipv.second.value());
      float ipXYSig = fabs(dxy_ipv.second.significance());

      //       std::cout << "Track with value, significance " << dxy_ipv.second.value() << "\t" << dxy_ipv.second.significance() << std::endl;

      TrajectoryStateOnSurface pca = IPTools::closestApproachToJet(itk->impactPointState(), primary_vertex,
          direction, itk->field());

      // Skip tracks with invalid point-of-closest-approach
      GlobalPoint closestPoint;
      if (pca.isValid()) {
        closestPoint = pca.globalPosition();
      } else {
        continue;
      }

      // Skip tracks if point-of-closest-approach has -nan or nan x/y/z coordinates
      if ( ! ( std::isfinite(closestPoint.x()) && std::isfinite(closestPoint.y()) && std::isfinite(closestPoint.z()) ) )
        continue;

      TLorentzVector trackVector;
      trackVector.SetPxPyPzE(
          closestPoint.x() - primary_vertex.position().x(),
          closestPoint.y() - primary_vertex.position().y(),
          closestPoint.z() - primary_vertex.position().z(),
          itk->track().p());

      // Skip tracks with deltaR > 0.4 w.r.t. current jet
      float deltaR = trackVector.DeltaR(jetVector);
      // if (itrack==1) std::cout << "deltaR: " << deltaR << std::endl;
      if (deltaR > 0.4) continue;

      int algo = itk->track().algo();
      int origAlgo = itk->track().originalAlgo();
      int nHits = itk->numberOfValidHits();
      int nMissInnerHits = itk->hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_INNER_HITS);
      // std::cout << "nHits:" << nHits << std::endl;
      // std::cout << "nMissInnerHits:" << nMissInnerHits << std::endl;
      vec_algo           . push_back ( algo           ) ;
      vec_originalAlgo   . push_back ( origAlgo       ) ;
      vec_nHits          . push_back ( nHits          ) ;
      vec_nMissInnerHits . push_back ( nMissInnerHits ) ;
      vec_ipXY           . push_back ( ipXY           ) ;
      vec_ipXYSig        . push_back ( ipXYSig        ) ;
      ipVector.push_back(fabs(dxy_ipv.second.value()));
      logIpSig.push_back(TMath::Log(fabs(dxy_ipv.second.significance())));


      misshits += itk->track().numberOfLostHits();
      tot_dsz += itk->track().dsz(theBeamSpot_->position());

      if (pca.isValid()) {
        r_pca.push_back(closestPoint.perp());
      }

      //           if (dxy_ipv.second.value() < ipCut) continue;
      if (dxy_ipv.second.significance() < minSigDispTrack) continue;

      // logSum += dxy_ipv.second.significance();

      Measurement1D ip2d = dxy_ipv.second;
      float r = 100*3.3*itk->track().pt()/3.8;
      float gLxy = ip2d.value()/sin(itk->track().phi()-jet.phi())*(1-2.5*fabs(ip2d.value())/r);

      sum_Lxy += gLxy;

      //           pTxIPxSig += itk->track().pt() * dxy_ipv.second.value() * dxy_ipv.second.significance() / jet->pt();
      pTxIPxSig += itk->track().pt() * dxy_ipv.second.significance() ;
      //           --promptTracks;
      ++nDispTracks;
      trackFrac += itk->track().pt();

      ipTracks.push_back(*itk);

      if (dxy_ipv.second.value() > 0.1) ip0_1++;
      if (dxy_ipv.second.value() > 1.0) ip1_0++;
      if (dxy_ipv.second.value() > 10.) ip10_++;

    }
    trackFrac /= jet.pt();

    std::sort(logIpSig.begin(),logIpSig.end());
    std::sort(ipVector.begin(),   ipVector.end());
    std::sort(r_pca.begin(), r_pca.end());
    if (logIpSig.size() != 0) {
      medianLogIpSig = logIpSig[logIpSig.size()/2];
      if (logIpSig[logIpSig.size()/2] > logTagCut) {
        ++nTags;
      }
    }
  }
  // std::cout << "itrack: " << itrack << std::endl;
  // std::cout << "vec_ipXY.size(): " << vec_ipXY.size() << std::endl;

  // Calculate jet alpha_max
  double alpha_max = 0.;
  {
    reco::TrackRefVector trackRefs = jet.getTrackRefs();
    // Loop over all tracks and calculate scalar pt-sum of all tracks in current jet
    double jet_pt_sum = 0.;
    for (reco::TrackRefVector::iterator ijt = trackRefs.begin(); ijt != trackRefs.end(); ++ijt) {
      jet_pt_sum += (*ijt)->pt();
    } // End of track loop

    auto ipv_chosen = primary_vertices_->end(); // iterator to chosen primary vertex
    double max_vertex_pt_sum = 0.; // scalar pt contribution of vertex to jet
    // Loop over all PVs and choose the one with highest scalar pt contribution to jet
    for (auto ipv = primary_vertices_->begin(); ipv != primary_vertices_->end(); ++ipv) {
      double vertex_pt_sum = 0.; // scalar pt contribution of vertex to jet
      // std::cout << "New vertex\n";
      for (reco::TrackRefVector::iterator ijt = trackRefs.begin(); ijt != trackRefs.end(); ++ijt) {
        double trackWeight = ipv->trackWeight(*ijt);
        // std::cout << "Track weight: " << trackWeight << std::endl;
        if (trackWeight > 0) vertex_pt_sum += (*ijt)->pt();
      } // End of track loop
      if (vertex_pt_sum > max_vertex_pt_sum) {
        max_vertex_pt_sum = vertex_pt_sum;
        ipv_chosen = ipv;
        // Calculate alpha
        alpha_max = vertex_pt_sum / jet_pt_sum;
      }
    } // End of vertex loop
  }
  // std::cout<< "Jet alpha_max: " << alpha_max << std::endl;

  // Count number of dark pions
  int nDarkPions = 0;
  double minDist = 9999.;
  {
    if (!isData_) {
      for (auto gp = genParticlesH_->begin(); gp != genParticlesH_->end(); ++gp) {
        if (fabs(gp->pdgId()) != 4900111) continue;
        //                   std::cout << "genParticle pdgId, status = " << gp->status() << "\t" << gp->pdgId() << std::endl;
        //                   std::cout << "             mass = " << gp->mass() << std::endl;
        //                   std::cout << "               pT = " << gp->pt()   << std::endl;
        //                   std::cout << "      r, eta, phi = " << gp->vertex().r() << "\t" << gp->vertex().eta() << "\t" << gp->vertex().phi() << std::endl;
        TLorentzVector gpVector;
        gpVector.SetPtEtaPhiM(gp->pt(),gp->eta(),gp->phi(),gp->mass());
        double dist = jetVector.DeltaR(gpVector);
        if (dist < 0.4) {
          nDarkPions++;
        }
        if (dist < minDist) minDist = dist;
      }
    }
  }

  bool matched = false;
  // Find out if there is a GenJet with deltaR < 0.2 relative to current jet
  for ( auto igj = genJets_->begin(); igj != genJets_->end(); ++igj) {
    TLorentzVector gjVector;
    gjVector.SetPtEtaPhiM(igj->pt(), igj->eta(), igj->phi(), 0.);
    if (gjVector.DeltaR(jetVector) < 0.2) {
      matched = true;
      break;
    }
  }

  int dtHits = 0;
  // Count DT hits with position vector that has deltaR < 0.5 relative to current jet
  for (std::vector<GlobalPoint>::iterator idt = dt_points_.begin(); idt != dt_points_.end(); ++idt) {
    TLorentzVector pointVector;
    pointVector.SetPtEtaPhiM(idt->perp(), idt->eta(), idt->phi(), 0.);
    if (pointVector.DeltaR(jetVector) < 0.5) ++dtHits;
  }

  int cscHits = 0;
  // Count CSC hits with position vector that has deltaR < 0.5 relative to current jet
  for (std::vector<GlobalPoint>::iterator icp = csc_points_.begin(); icp != csc_points_.end(); ++icp) {
    TLorentzVector pointVector;
    pointVector.SetPtEtaPhiM(icp->perp(), icp->eta(), icp->phi(), 0.);
    if (pointVector.DeltaR(jetVector) < 0.5) ++cscHits;
  }

  // Calculate number of vertices within jet cone, and median of displacement
  int matchedVertices = 0;
  float medianSVradius = 0.;
  {
    std::vector<float> radiusVector;
    // Loop over a given vertex vector/collection
    for (reco::VertexCollection::iterator ivtx = selectedSecondaryVertices_.begin(); ivtx != selectedSecondaryVertices_.end(); ++ivtx) {
      TLorentzVector vertexPosition;  vertexPosition.SetPtEtaPhiM(ivtx->position().r(),ivtx->position().eta(),ivtx->position().phi(),0.);
      if (vertexPosition.DeltaR(jetVector) < 0.4) ++matchedVertices;
      if (ivtx->normalizedChi2() > 15.) continue;
      radiusVector.push_back(ivtx->position().r());
      TLorentzVector cand;
      // Build TransientTrack vector from refitted tracks
      std::vector<reco::TransientTrack> transRefitTracks;
      for (size_t itrack = 0; itrack < ivtx->refittedTracks().size(); itrack++) {
        transRefitTracks.push_back(transienttrackbuilder_->build(ivtx->refittedTracks()[itrack]));
      }
    }
    std::sort(radiusVector.begin(), radiusVector.end());
    if (radiusVector.size() != 0) {
      medianSVradius = radiusVector.at(radiusVector.size()/2);
    }
  }


  otree.jets_pt             .push_back( jet.pt()                          );
  otree.jets_eta            .push_back( jet.eta()                         );
  otree.jets_phi            .push_back( jet.phi()                         );
  otree.jets_cef            .push_back( jet.chargedEmEnergyFraction()     );
  otree.jets_nef            .push_back( jet.neutralEmEnergyFraction()     );
  otree.jets_chf            .push_back( jet.chargedHadronEnergyFraction() );
  otree.jets_nhf            .push_back( jet.neutralHadronEnergyFraction() );
  otree.jets_phf            .push_back( jet.photonEnergyFraction()        );
  otree.jets_nPromptTracks  .push_back( nPromptTracks                     );
  otree.jets_nDispTracks    .push_back( nDispTracks                       );
  otree.jets_nSV            .push_back( matchedVertices                   );
  otree.jets_medianLogIpSig .push_back( medianLogIpSig                    );
  // otree.jets_missHits       .push_back( misshits                           );
  // otree.jets_muonHits       .push_back( dtHits+cscHits                     );
  otree.jets_alphaMax       .push_back( alpha_max                         );
  otree.jets_nDarkPions     .push_back( nDarkPions                        );
  otree.jets_minDRDarkPion  .push_back( minDist                           );
  otree.tracks_algo           .push_back ( vec_algo           ) ;
  otree.tracks_originalAlgo   .push_back ( vec_originalAlgo   ) ;
  otree.tracks_nHits          .push_back ( vec_nHits          ) ;
  otree.tracks_nMissInnerHits .push_back ( vec_nMissInnerHits ) ;
  otree.tracks_ipXY           .push_back ( vec_ipXY           ) ;
  otree.tracks_ipXYSig        .push_back ( vec_ipXYSig        ) ;


}

reco::VertexCollection
EmergingJetAnalyzer::selectSecondaryVertices (edm::Handle<reco::VertexCollection> secondary_vertices) const {
  const reco::Vertex& primary_vertex = primary_vertices_->at(0);
  const int minDispSv = 0.1; // Minimum transverse displacement for a secondary vertex to be selected
  reco::VertexCollection selectedVertices;
  for (size_t ivx = 0; ivx < secondary_vertices->size(); ++ivx) {
    float dx = primary_vertex.position().x() - secondary_vertices->at(ivx).position().x();
    float dy = primary_vertex.position().y() - secondary_vertices->at(ivx).position().y();
    // If displacement from primary_vertex is greater than minDispSv
    if (TMath::Sqrt( dx*dx + dy*dy ) > minDispSv) {
        selectedVertices.push_back(secondary_vertices->at(ivx));
      }
  }
  return selectedVertices;
}

double
EmergingJetAnalyzer::calculatePt2Sum (const reco::Vertex& vertex) const {
  // Modified from reco::Vertex::p4()
  double sum = 0.;
  double pt2 = 0.;

  if(vertex.hasRefittedTracks()) {
    for(std::vector<reco::Track>::const_iterator iter = vertex.refittedTracks().begin();
        iter != vertex.refittedTracks().end(); ++iter) {
      pt2 = iter->pt();
      sum += pt2;
    }
  }
  else
    {
      for(std::vector<reco::TrackBaseRef>::const_iterator iter = vertex.tracks_begin();
          iter != vertex.tracks_end(); iter++) {
        pt2 = (*iter)->pt();
        sum += pt2;
      }
    }
  return sum;
}

//define this as a plug-in
DEFINE_FWK_MODULE(EmergingJetAnalyzer);
