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

#include "EmergingJetAnalysis/EmergingJetAnalyzer/interface/OutputTree.h"
//
// class declaration
//

class EmergingJetAnalyzer : public edm::EDAnalyzer {
  public:
    explicit EmergingJetAnalyzer(const edm::ParameterSet&);
    ~EmergingJetAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    OutputTree otree;


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

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
    int t_run, t_event, t_ls, t_bx;

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
  t_tree->Branch("run",&t_run);
  t_tree->Branch("event",&t_event);
  t_tree->Branch("LS",&t_ls);
  t_tree->Branch("BX",&t_bx);

  t_tree->Branch("met",&t_met);

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

  //   float ipCut = 0.05;


  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  edm::Handle<reco::PFJetCollection> pfjetH;
  // iEvent.getByLabel("ak4PFJetsCHS", pfjetH);
  iEvent.getByToken(jetCollectionToken_, pfjetH);
  auto selectedJets = *pfjetH;

  Handle<reco::BeamSpot> theBeamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", theBeamSpotHandle);
  const reco::BeamSpot* theBeamSpot = theBeamSpotHandle.product();

  //   edm::Handle<reco::CaloJetCollection> caloJetH;
  //   iEvent.getByLabel("ak5CaloJets", caloJetH);

  edm::Handle<reco::TrackCollection> genTrackH;
  iEvent.getByLabel("generalTracks", genTrackH);

  std::vector<reco::TransientTrack> generalTracks;
  generalTracks = theB->build(genTrackH);

  edm::Handle<reco::TrackCollection> sdmH;
  iEvent.getByLabel("displacedStandAloneMuons",sdmH);

  std::vector<reco::TransientTrack> standaloneDisplacedMuons;
  standaloneDisplacedMuons = theB->build(sdmH);

  edm::Handle<reco::PFMETCollection> pfmet;
  iEvent.getByLabel("pfMet",pfmet);

  edm::Handle<reco::VertexCollection> primary_vertices;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS",primary_vertices);
  const reco::Vertex& primary_vertex = primary_vertices->at(0);

  edm::Handle<reco::VertexCollection> secondary_vertices;
  iEvent.getByLabel("inclusiveSecondaryVertices",secondary_vertices);



  //   /*
  //   */

  float dx, dy;
  int disp_SV = 0;
  reco::VertexCollection theVertices;
  //   for (reco::VertexCollection::iterator ivx = secondary_vertices->begin(); ivx != secondary_vertices->end(); ++ivx) {
  for (size_t ivx = 0; ivx < secondary_vertices->size(); ++ivx) {
    //        std::cout << "ivx " << ivx << std::endl;
    dx = primary_vertex.position().x() - secondary_vertices->at(ivx).position().x();
    dy = primary_vertex.position().y() - secondary_vertices->at(ivx).position().y();
    //        dy = ivx->position.y() - primary_vertex.position().y();
    //        std::cout << "TMath::Sqrt( dx*dx + dy*dy ) " << TMath::Sqrt( dx*dx + dy*dy ) << std::endl;
    h_rSV->Fill(TMath::Sqrt( dx*dx + dy*dy ));
    if (TMath::Sqrt( dx*dx + dy*dy ) > 0.1) ++disp_SV;
    theVertices.push_back(secondary_vertices->at(ivx));
  }
  h_nVtxInEvent->Fill(theVertices.size());
  // testing some event selection criteria here
  //   if ( (pfmet->begin()->et() < 150.) && (disp_SV < 10.) ) return;


  std::vector<reco::GenJet> emergingGenJets;
  if (!isData_) {
    edm::Handle<std::vector<reco::GenMET> > genMetH;
    iEvent.getByLabel("genMetTrue", genMetH);
    h_genHT->Fill(genMetH->at(0).sumEt());

    edm::Handle<reco::GenParticleCollection> genParticlesH;
    iEvent.getByLabel("genParticles", genParticlesH);


    for (reco::GenParticleCollection::const_iterator gp = genParticlesH->begin(); gp != genParticlesH->end(); ++gp) {
      if (gp->numberOfMothers() < 1) continue;
      if (fabs(gp->mother()->pdgId()) != 4900111) continue;
      //                   std::cout << "genParticle pdgId, status = " << gp->status() << "\t" << gp->pdgId() << std::endl;
      //                   std::cout << "             mass = " << gp->mass() << std::endl;
      //                   std::cout << "               pT = " << gp->pt()   << std::endl;
      //                   std::cout << "      r, eta, phi = " << gp->vertex().r() << "\t" << gp->vertex().eta() << "\t" << gp->vertex().phi() << std::endl;
      h_rDecayGen->Fill(gp->vertex().r());

      // find closest jet
      float minDR = 999.;
      reco::PFJetCollection::iterator bestJet;
      for (reco::PFJetCollection::iterator ij = selectedJets.begin(); ij != selectedJets.end(); ++ij) {
        TLorentzVector tlv_jet;
        tlv_jet.SetPtEtaPhiM(ij->pt(),ij->eta(),ij->phi(),0.);
        TLorentzVector tlv_gen;
        tlv_gen.SetPtEtaPhiM(gp->pt(),gp->eta(),gp->phi(),gp->mass());

        float dist = tlv_jet.DeltaR(tlv_gen);
        if (dist < minDR) {
          bestJet = ij;
          minDR   = dist;
        }
      }
      h_dr_disp_jet->Fill(minDR);

      //           std::cout << "   Nearest jet " << std::endl;
      //           std::cout << "      dr pt eta phi " << minDR << "\t" << bestJet->pt() << "\t" << bestJet->eta() << "\t" << bestJet->eta() << "\t" << bestJet->phi() << std::endl;
    }

    edm::Handle<std::vector<reco::GenJet> >  genJetH;
    iEvent.getByLabel("ak4GenJets",   genJetH);

    int nPiDJets = 0;

    for (std::vector<reco::GenJet>::const_iterator igj = genJetH->begin(); igj != genJetH->end(); igj++) {
      if (igj->pt() < 50.) continue;
      //           std::cout << "New GenJet -------------" << std::endl; 
      int nDarkPions = 0;
      std::vector <const reco::GenParticle*> constituents = igj->getGenConstituents();
      for (std::vector<const reco::GenParticle*>::iterator ic = constituents.begin(); ic != constituents.end(); ++ic) {
        //                std::cout << "\t" << (*ic)->pdgId() << "\t" << (*ic)->mother()->pdgId() << std::endl;
        if ((*ic)->mother() != NULL) {  
          if (TMath::Abs((*ic)->mother()->pdgId()) == 4900111) nDarkPions++;
        }
      } 

      h_nDarkPions_genJet->Fill(nDarkPions);
      if (nDarkPions > 1) {
        ++nPiDJets;
        emergingGenJets.push_back(*igj);
      }


    }

    h_nJets_gt1piD->Fill(nPiDJets);

  }

  t_run   = iEvent.id().run();
  t_event = iEvent.id().event();
  t_ls    = iEvent.id().luminosityBlock();
  t_bx    = iEvent.bunchCrossing();

  int muonRecHits = 0;

  Handle<DTRecSegment4DCollection> segments4D;
  iEvent.getByLabel("dt4DSegments",segments4D);

  GlobalPoint segpos;
  DTRecSegment4DCollection::const_iterator iseg;

  std::vector<GlobalPoint> dt_points;

  for (iseg = segments4D->begin(); iseg != segments4D->end(); ++iseg) {

    const GeomDet* geomDet = theTrackingGeometry->idToDet((*iseg).geographicalId());
    segpos = geomDet->toGlobal((*iseg).localPosition());
    dt_points.push_back(segpos);

  }

  muonRecHits += segments4D->size();
  //   h_DTsegments->Fill(muonRecHits);

  Handle<CSCSegmentCollection> csc_segments;
  iEvent.getByLabel("cscSegments",csc_segments);

  CSCSegmentCollection::const_iterator icsc;
  std::vector<GlobalPoint> csc_points;

  for (icsc = csc_segments->begin(); icsc != csc_segments->end(); ++icsc) {

    const GeomDet* geomDet = theTrackingGeometry->idToDet((*icsc).geographicalId());
    segpos = geomDet->toGlobal((*icsc).localPosition());
    csc_points.push_back(segpos);

  }
  //   h_CSCsegments->Fill(csc_segments->size());

  Handle<RPCRecHitCollection> rpc_hits;
  iEvent.getByLabel("rpcRecHits",rpc_hits);
  //   h_RPChits->Fill(rpc_hits->size());

  h_nPV->Fill(primary_vertices->size());
  h_nSV->Fill(secondary_vertices->size());


  h_nSV_disp->Fill(disp_SV);

  h_MET->Fill(pfmet->begin()->et());
  t_met = pfmet->begin()->et();
  float dphi = 999.;
  for (reco::PFJetCollection::iterator ijet = selectedJets.begin(); ijet != selectedJets.end(); ++ijet) {
    float dphi_i = ijet->phi() - pfmet->begin()->phi();
    while (dphi_i > TMath::Pi()) dphi_i -= 2*TMath::Pi();
    while (dphi_i < -TMath::Pi()) dphi_i +=2*TMath::Pi();
    if (TMath::Abs(dphi_i) < TMath::Abs(dphi)) dphi = dphi_i;
  }
  h_dPhiMetJet->Fill(dphi);

  /*
     edm::Handle<reco::GenParticleCollection> genParticlesH;
     iEvent.getByLabel("genParticles", genParticlesH);

  //   std::vector<math::XYZPoint> decayVertices;

  for (reco::GenParticleCollection::const_iterator gp = genParticlesH->begin(); gp != genParticlesH->end(); ++gp) {
  if (fabs(gp->pdgId()) != 4900111) continue;
  std::cout << "genParticle pdgId = " << gp->pdgId() << std::endl;
  std::cout << "             mass = " << gp->mass() << std::endl;
  std::cout << "               pT = " << gp->pt()   << std::endl;

  math::XYZPoint decayVertex = gp->daughter(0)->vertex();
  //        decayVertices.push_back(decayVertex);
  std::cout << "         vertex r = " << decayVertex.r() << std::endl;
  std::cout << "            x,y,z = "
  << "\t" << decayVertex.x()
  << "\t" << decayVertex.y()
  << "\t" << decayVertex.z()
  <<std::endl;
  std::cout << "          daughters" << std::endl;
  for (size_t id = 0; id < gp->numberOfDaughters(); ++id) {

  std::cout << "          daughter " << id << " " << gp->daughter(id)->pdgId() << std::endl;
  }

  }
  //
  */

  //   for (std::vector<reco::TransientTrack>::iterator itk = generalTracks.begin(); itk < generalTracks.end(); ++itk) {
  //   }


  // try my own vertex reco


  //   return;

  int ijet = 0;
  int nTags = 0;
  float logTagCut    = 2.;
  float logSum = 0.;
  std::vector<float> vec_medianIpSig;
  std::vector<float> vec_promptTracks;
  std::vector<float> vec_medianRpca;
  int nMatchedGen = 0;
  for ( reco::PFJetCollection::const_iterator jet = selectedJets.begin(); jet != selectedJets.end(); ++jet ) {


    //       if (jet->pt() < 50.) continue;
    //       if (fabs(jet->eta()) > 2.5) continue;

    h_jetPtVsIndex->Fill(jet->pt(), ijet);
    TLorentzVector jetVector; jetVector.SetPtEtaPhiM(jet->pt(),jet->eta(),jet->phi(),0.);

    bool matched = false;

    for ( std::vector<reco::GenJet>::iterator igj = emergingGenJets.begin(); igj != emergingGenJets.end() && matched == false; ++igj) {
      TLorentzVector gjVector;
      gjVector.SetPtEtaPhiM(igj->pt(), igj->eta(), igj->phi(), 0.);
      if (gjVector.DeltaR(jetVector) < 0.2) matched = true;
    }

    if (matched) nMatchedGen++;
    //       std::cout << "Jet pt,eta,phi " << jet->pt() << "\t" << jet->eta() << "\t" << jet->phi() << std::endl;

    //std::cout << "PFJet pt,eta,phi" << "\t" << jet->pt() << "\t" << jet->eta() << "\t" << jet->phi() << std::endl;
    //std::cout << "\t neutralEM:  " << jet->neutralEmEnergyFraction() << std::endl;
    //std::cout << "\tneutralHad:  " << jet->neutralHadronEnergyFraction() << std::endl;
    float trackFrac = 0.;
    int matchedTracks = 0;

    int dtHits = 0;
    for (std::vector<GlobalPoint>::iterator idt = dt_points.begin(); idt != dt_points.end(); ++idt) {
      TLorentzVector pointVector;
      pointVector.SetPtEtaPhiM(idt->perp(), idt->eta(), idt->phi(), 0.);
      if (pointVector.DeltaR(jetVector) < 0.5) ++dtHits;
    }
    h_DTsegments->Fill(dtHits);

    int cscHits = 0;
    for (std::vector<GlobalPoint>::iterator icp = csc_points.begin(); icp != csc_points.end(); ++icp) {
      TLorentzVector pointVector;
      pointVector.SetPtEtaPhiM(icp->perp(), icp->eta(), icp->phi(), 0.);
      if (pointVector.DeltaR(jetVector) < 0.5) ++cscHits;
    }
    h_CSCsegments->Fill(cscHits);


    //       int pfTracks = jet->getTrackRefs().size();
    int promptTracks = 0;
    for (reco::TrackRefVector::iterator ijt = jet->getTrackRefs().begin(); ijt != jet->getTrackRefs().end(); ++ijt) {
      reco::TransientTrack itk = theB->build(*ijt);
      if (itk.track().pt() < 1.) continue;
      auto d3d_ipv = IPTools::absoluteImpactParameter3D(itk, primary_vertex);
      if (d3d_ipv.second.significance() < 3.) promptTracks++;
      //       std::cout << "Track with value, significance " << dxy_ipv.second.value() << "\t" << dxy_ipv.second.significance() << std::endl;
      //           if (dxy_ipv.second.value() > ipCut) continue;
      h_3Dip->Fill(TMath::Log(d3d_ipv.second.value()));
      h_3DipSig->Fill(d3d_ipv.second.significance());
      if (matched) {
        h_3Dip_genMatched->Fill(TMath::Log(d3d_ipv.second.value()));
        h_3DipSig_genMatched->Fill(d3d_ipv.second.significance());
      }
      //           if (dxy_ipv.second.significance() > ipSig) continue;

      //           trackFrac += itk.track().pt();
    } 
    // now loop over the displaced tracks, checking which match the jet geometrically
    std::vector<float> ipVector;
    float ipSig = 3.;
    int ip0_1 = 0.;
    int ip1_0 = 0.;
    int ip10_ = 0.;
    int dispTracks   = 0;
    float pTxIPxSig = 0.;
    std::vector<reco::TransientTrack> ipTracks;

    math::XYZVector jetMomentum = jet->momentum();
    GlobalVector direction(jetMomentum.x(), jetMomentum.y(), jetMomentum.z());

    int misshits = 0;
    float tot_dsz = 0.;
    std::vector<float> r_pca;
    std::vector<float> logIpSig;
    float sum_Lxy = 0.;
    for (std::vector<reco::TransientTrack>::iterator itk = generalTracks.begin(); itk != generalTracks.end(); ++itk) {

      if (itk->track().pt() < 1.) continue;

      //           if (itk->track().extra().isNull()) continue;

      //           TLorentzVector trackVector; 
      //           trackVector.SetPxPyPzE(
      //                   itk->innermostMeasurementState().globalPosition().x() - primary_vertex.position().x(), 
      //                   itk->innermostMeasurementState().globalPosition().y() - primary_vertex.position().y(), 
      //                   itk->innermostMeasurementState().globalPosition().z() - primary_vertex.position().z(), 
      //                   itk->track().p());
      //           if (trackVector.DeltaR(jetVector) > 0.5) continue;


      //           ++promptTracks;
      //           TrackDetMatchInfo info = m_trackAssociator.associate(iEvent, iSetup, itk->track(), m_trackParameters);
      //           TLorentzVector impactVector; impactVector.SetPtEtaPhiM(itk->track().pt(),info.trkGlobPosAtEcal.eta(),info.trkGlobPosAtEcal.phi(),0.);

      //           std::cout << "Impact point at ECAL eta, phi, dR " 
      //               << info.trkGlobPosAtEcal.eta() << "\t" 
      //               << info.trkGlobPosAtEcal.phi() << "\t" 
      //               << jetVector.DeltaR(impactVector) << std::endl;
      //           if (jetVector.DeltaR(impactVector) > 0.4) continue;

      auto dxy_ipv = IPTools::absoluteTransverseImpactParameter(*itk, primary_vertex);

      //       std::cout << "Track with value, significance " << dxy_ipv.second.value() << "\t" << dxy_ipv.second.significance() << std::endl;


      TrajectoryStateOnSurface pca = IPTools::closestApproachToJet(itk->impactPointState(), primary_vertex, 
          direction, itk->field());

      GlobalPoint closestPoint;
      if (pca.isValid()) {
        closestPoint = pca.globalPosition();
      } else {
        continue;
      }

      TLorentzVector trackVector;
      trackVector.SetPxPyPzE(
          closestPoint.x() - primary_vertex.position().x(),
          closestPoint.y() - primary_vertex.position().y(),
          closestPoint.z() - primary_vertex.position().z(),
          itk->track().p());
      if (pca.isValid()) {
        h_jetPCAdist->Fill( trackVector.DeltaR(jetVector) , closestPoint.perp() );
      }
      h_pca_dr_jet->Fill(trackVector.DeltaR(jetVector));
      h_pca_r->Fill(closestPoint.perp());
      if (trackVector.DeltaR(jetVector) > 0.4) continue;
      h_trackPtInJet->Fill(itk->track().pt());
      if (matched) h_trackPtInJet_genMatched->Fill(itk->track().pt());
      ipVector.push_back(fabs(dxy_ipv.second.value()));
      logIpSig.push_back(TMath::Log(fabs(dxy_ipv.second.significance())));

      misshits += itk->track().numberOfLostHits();
      tot_dsz += itk->track().dsz(theBeamSpot->position());

      if (pca.isValid()) {
        r_pca.push_back(closestPoint.perp()); 
        h_rPCA->Fill(closestPoint.perp());
      }

      //           if (dxy_ipv.second.value() < ipCut) continue;
      if (dxy_ipv.second.significance() < ipSig) continue;

      h_dz_vs_dxy->Fill(TMath::Log(fabs(dxy_ipv.second.value())),TMath::Log(fabs(itk->track().dz() - primary_vertex.position().z())));
      if (matched) { 
        h_dz_vs_dxy_genMatched->Fill(TMath::Log(fabs(dxy_ipv.second.value())),TMath::Log(fabs(itk->track().dz() - primary_vertex.position().z()))); 
        h_logTrackDz_genMatched->Fill(TMath::Log(fabs(itk->track().dz() - primary_vertex.position().z())));
        h_logTrackDxy_genMatched->Fill(TMath::Log(fabs(dxy_ipv.second.value())));
      }
      h_logTrackDz->Fill(TMath::Log(fabs(itk->track().dz() - primary_vertex.position().z())));
      h_logTrackDxy->Fill(TMath::Log(fabs(dxy_ipv.second.value())));
      logSum += dxy_ipv.second.significance();

      Measurement1D ip2d = dxy_ipv.second;
      float r = 100*3.3*itk->track().pt()/3.8;
      float gLxy = ip2d.value()/sin(itk->track().phi()-jet->phi())*(1-2.5*fabs(ip2d.value())/r);

      sum_Lxy += gLxy;

      //           pTxIPxSig += itk->track().pt() * dxy_ipv.second.value() * dxy_ipv.second.significance() / jet->pt();
      pTxIPxSig += itk->track().pt() * dxy_ipv.second.significance() ;
      //           --promptTracks;
      ++dispTracks;
      //           h_dr_jet_track->Fill(jetVector.DeltaR(impactVector));
      //           h_deltaEtaJetTrack->Fill(fabs(jetVector.Eta() - impactVector.Eta()));
      //           h_deltaPhiJetTrack->Fill(fabs(jetVector.DeltaPhi(impactVector)));
      trackFrac += itk->track().pt();
      //       if (itk->track().dxy(*theBeamSpot) < ipCut && itk->track().dxy(*theBeamSpot) > -ipCut) continue;
      //       if (itk->track().dxy(*theBeamSpot) / itk->track().dxyError() < ipSig) continue;
      //       std::cout << "\tq dxy " << itk->track().charge() << "\t" << itk->track().dxy(*theBeamSpot) << std::endl;
      //       std::cout << "\tImpact parameter significance " << itk->track().dxy(*theBeamSpot) / itk->track().dxyError() << std::endl;

      ipTracks.push_back(*itk);

      //           std::cout << "dxy_ipv.second.value() " << dxy_ipv.second.value() << std::endl;

      h_pTdispTracks->Fill(itk->track().pt());
      if (dxy_ipv.second.value() > 0.1) ip0_1++; 
      if (dxy_ipv.second.value() > 1.0) ip1_0++; 
      if (dxy_ipv.second.value() > 10.) ip10_++; 

      ++matchedTracks;
    }

    h_Lxy_perjet->Fill(TMath::Log(sum_Lxy));

    h_dsz_perJet->Fill(tot_dsz);
    h_nMissingHitsPerJet->Fill(misshits);
    h_jetPromptTracks->Fill(promptTracks);
    vec_promptTracks.push_back(promptTracks);
    h_jetDispTracks->Fill(dispTracks);
    h_eventIpTracks0_1->Fill(ip0_1);
    h_eventIpTracks1_0->Fill(ip1_0);
    h_eventIpTracks10_->Fill(ip10_);

    std::sort(logIpSig.begin(),logIpSig.end());
    std::sort(ipVector.begin(),   ipVector.end());
    std::sort(r_pca.begin(), r_pca.end());
    float medianIpSig = 0.;
    if (logIpSig.size() != 0) {
      h_medianLogIpSig->Fill(logIpSig[logIpSig.size()/2]);
      vec_medianIpSig.push_back(logIpSig[logIpSig.size()/2]);
      medianIpSig = logIpSig[logIpSig.size()/2];
      if (r_pca.size() > 0) vec_medianRpca.push_back(r_pca[r_pca.size()/2]);
      h_medianIp->Fill(ipVector[ipVector.size()/2]);
      if (logIpSig[logIpSig.size()/2] > logTagCut) {
        ++nTags;
      }
    }
    trackFrac /= jet->pt();
    h_trackFrac->Fill(trackFrac);

    h_pTxIPxSig->Fill(TMath::Log(pTxIPxSig));

    h_jetIpTracks->Fill(matchedTracks);


    //       AdaptiveVertexReconstructor avr (2.0, 6.0, 0.5, true );
    //       std::vector<TransientVertex> theVertices = avr.vertices(ipTracks);
    std::vector<TLorentzVector>  vertexVectors;
    int matchedVertices = 0;
    std::vector<float> radiusVector;
    for (reco::VertexCollection::iterator ivtx = theVertices.begin(); ivtx != theVertices.end(); ++ivtx) {
      TLorentzVector vertexPosition;  vertexPosition.SetPtEtaPhiM(ivtx->position().r(),ivtx->position().eta(),ivtx->position().phi(),0.);
      h_vertexDrJet->Fill(vertexPosition.DeltaR(jetVector));
      if (vertexPosition.DeltaR(jetVector) < 0.15) ++matchedVertices;
      if (ivtx->normalizedChi2() > 15.) continue;
      h_rVtx->Fill(ivtx->position().r());
      radiusVector.push_back(ivtx->position().r());
      h_chi2Vtx->Fill(ivtx->normalizedChi2());
      h_nTracksVtx->Fill(ivtx->refittedTracks().size());
      TLorentzVector cand;
      // loop over the tracks
      //           std::vector<reco::Track> vecRefittedTracks = ivtx->refittedTracks();
      std::vector<reco::TransientTrack> transRefitTracks;
      for (size_t itrack = 0; itrack < ivtx->refittedTracks().size(); itrack++) {
        transRefitTracks.push_back(theB->build(ivtx->refittedTracks()[itrack]));
      }
      for (std::vector<reco::TransientTrack>::const_iterator itk = transRefitTracks.begin(); itk != transRefitTracks.end(); ++itk) {
        TrajectoryStateClosestToPoint trajectory = itk->trajectoryStateClosestToPoint(GlobalPoint(ivtx->position().x(),ivtx->position().y(),ivtx->position().z()));
        //            int charge = itk->track().charge();
        GlobalVector momentum = trajectory.momentum();
        TLorentzVector trackVector;  trackVector.SetPtEtaPhiM(momentum.perp(), momentum.eta(), momentum.phi(), 0.13957);
        cand += trackVector;
      }
      vertexVectors.push_back(cand);
      h_massVtx->Fill(cand.M());
      h_ptVtx->Fill(cand.Pt());
      h_log_rz_vertex->Fill(TMath::Log(TMath::Abs(ivtx->position().z())),TMath::Log(TMath::Abs(ivtx->position().r())));
    }
    float medianSVradius = 0.;
    std::sort(radiusVector.begin(), radiusVector.end());
    if (radiusVector.size() != 0) {
      medianSVradius = radiusVector.at(radiusVector.size()/2);
    }
    h_medianSVradius->Fill(medianSVradius);
    h_jetMatchedVertices->Fill(matchedVertices);

    // displaced standalone muons
    int matchedMuons = 0;
    //       for (std::vector<reco::TransientTrack>::iterator itk = standaloneDisplacedMuons.begin(); itk < standaloneDisplacedMuons.end(); ++itk) {

    //           if (!itk->innermostMeasurementState().isValid()) continue;
    //           TLorentzVector trackVector; trackVector.SetPtEtaPhiM(itk->innermostMeasurementState().globalPosition().perp(), itk->innermostMeasurementState().globalPosition().eta(), itk->innermostMeasurementState().globalPosition().phi(), 0.13957);
    //           h_muonDrJet->Fill(trackVector.DeltaR(jetVector));
    //           if (trackVector.DeltaR(jetVector) < 0.5) ++matchedMuons;

    //       }
    h_jetDispMuons->Fill(matchedMuons);


    h_jet_cef->Fill(jet->chargedEmEnergyFraction());
    h_jet_nef->Fill(jet->neutralEmEnergyFraction());
    h_jet_chf->Fill(jet->chargedHadronEnergyFraction());
    if (jet->chargedHadronEnergyFraction() < 0.3) nTags++;
    h_jet_nhf->Fill(jet->neutralHadronEnergyFraction());
    h_jet_phf->Fill(jet->photonEnergyFraction());

    if (promptTracks < 5) {
      h_jet_cef_lt5prtk->Fill(jet->chargedEmEnergyFraction());
      h_jet_nef_lt5prtk->Fill(jet->neutralEmEnergyFraction());
      h_jet_chf_lt5prtk->Fill(jet->chargedHadronEnergyFraction());
      h_jet_nhf_lt5prtk->Fill(jet->neutralHadronEnergyFraction());
      h_jet_phf_lt5prtk->Fill(jet->photonEnergyFraction());
    }

    h_jetNeutralFraction->Fill(jet->neutralHadronEnergyFraction() + jet->neutralEmEnergyFraction());
    if (ijet == 0) {
      jet0_pt = jet->pt();
      jet0_eta = jet->eta();
      jet0_promptTracks = promptTracks;
      jet0_dispTracks   = dispTracks;
      jet0_nSV          = matchedVertices;
      jet0_cef          = jet->chargedEmEnergyFraction();
      jet0_nef          = jet->neutralEmEnergyFraction();
      jet0_chf          = jet->chargedHadronEnergyFraction();
      jet0_nhf          = jet->neutralHadronEnergyFraction();
      jet0_phf          = jet->photonEnergyFraction();
      jet0_medianLogIpSig = medianIpSig;
      jet0_missHits     = misshits;
      jet0_muonHits     = dtHits+cscHits;

    } else if (ijet == 1) {
      jet1_pt = jet->pt();
      jet1_eta = jet->eta();
      jet1_promptTracks = promptTracks;
      jet1_dispTracks   = dispTracks;
      jet1_nSV          = matchedVertices;
      jet1_cef          = jet->chargedEmEnergyFraction();
      jet1_nef          = jet->neutralEmEnergyFraction();
      jet1_chf          = jet->chargedHadronEnergyFraction();
      jet1_nhf          = jet->neutralHadronEnergyFraction();
      jet1_phf          = jet->photonEnergyFraction();
      jet1_medianLogIpSig = medianIpSig;
      jet1_missHits     = misshits;
      jet1_muonHits     = dtHits+cscHits;

    } else if (ijet == 2) {
      jet2_pt = jet->pt();
      jet2_eta = jet->eta();
      jet2_promptTracks = promptTracks;
      jet2_dispTracks   = dispTracks;
      jet2_nSV          = matchedVertices;
      jet2_cef          = jet->chargedEmEnergyFraction();
      jet2_nef          = jet->neutralEmEnergyFraction();
      jet2_chf          = jet->chargedHadronEnergyFraction();
      jet2_nhf          = jet->neutralHadronEnergyFraction();
      jet2_phf          = jet->photonEnergyFraction();
      jet2_medianLogIpSig = medianIpSig;
      jet2_missHits     = misshits;
      jet2_muonHits     = dtHits+cscHits;

    } else if (ijet == 3) {
      jet3_pt = jet->pt();
      jet3_eta = jet->eta();
      jet3_promptTracks = promptTracks;
      jet3_dispTracks   = dispTracks;
      jet3_nSV          = matchedVertices;
      jet3_cef          = jet->chargedEmEnergyFraction();
      jet3_nef          = jet->neutralEmEnergyFraction();
      jet3_chf          = jet->chargedHadronEnergyFraction();
      jet3_nhf          = jet->neutralHadronEnergyFraction();
      jet3_phf          = jet->photonEnergyFraction();
      jet3_medianLogIpSig = medianIpSig;
      jet3_missHits     = misshits;
      jet3_muonHits     = dtHits+cscHits;

    }

    ++ijet;

  }

  h_nMatchedGen->Fill(nMatchedGen);

  std::sort(vec_medianIpSig.begin(), vec_medianIpSig.end());
  //   h_medianLogIpSig_2highest->Fill(vec_medianIpSig[0] + vec_medianIpSig[1]);
  h_medianLogIpSig_2highest->Fill(vec_medianIpSig[3] + vec_medianIpSig[2]);

  std::sort(vec_promptTracks.begin(), vec_promptTracks.end());
  //   h_jetPromptTracks_2lowest->Fill(vec_promptTracks[vec_promptTracks.size()-1]+vec_promptTracks[vec_promptTracks.size()-2]);
  h_jetPromptTracks_2lowest->Fill(vec_promptTracks[0]+vec_promptTracks[1]);

  std::sort(vec_medianRpca.begin(), vec_medianRpca.end());
  //   for (std::vector<float>::iterator ir = vec_medianRpca.begin(); ir != vec_medianRpca.end(); ++ir) {
  //       std::cout << *ir << std::endl;
  //   }
  h_median_rPCA->Fill(vec_medianRpca[2]);  if (vec_medianRpca[2] > 0.5) nTags++;
  h_median_rPCA->Fill(vec_medianRpca[3]);  if (vec_medianRpca[3] > 0.5) nTags++;

  // check out the standalone muons

  int dispMuons = 0;
  for (std::vector<reco::TransientTrack>::iterator itk = standaloneDisplacedMuons.begin(); itk < standaloneDisplacedMuons.end(); ++itk) {
    //        std::cout << "Displaced muon pt eta phi dxy " 
    //            << "\t" << itk->track().pt()
    //            << "\t" << itk->track().eta()
    //            << "\t" << itk->track().phi()
    //            << "\t" << itk->track().dxy()
    //            << std::endl;
    ++dispMuons;
  }


  //   h_eventIpTracks->Fill(ipTracks.size());
  //   h_eventVertices->Fill(theVertices.size());
  h_eventDispMuons->Fill(dispMuons);
  h_tagsPerEvent->Fill(nTags);
  h_eventSumLogIp->Fill(TMath::Log(logSum));

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

//define this as a plug-in
DEFINE_FWK_MODULE(EmergingJetAnalyzer);
