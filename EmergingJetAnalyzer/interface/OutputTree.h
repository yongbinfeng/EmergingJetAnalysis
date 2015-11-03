#ifndef EmergingJetAnalysis_EmergingJetAnalyzer_OutputTree_h
#define EmergingJetAnalysis_EmergingJetAnalyzer_OutputTree_h

#include <vector>

#include "TTree.h"

using std::vector;

namespace emjet
{
  class OutputTree {
    public:
      OutputTree() { Init(); }
      void Init();
      void Branch(TTree* tree);
      // per event
      int run;
      int lumi;
      int event;
      int bx;
      float  met_pt;
      float  met_phi;
      // per jet
      vector<float>  jets_pt;
      vector<float>  jets_eta;
      vector<float>  jets_phi;
      vector<int>    jets_promptTracks;
      vector<int>    jets_dispTracks;
      vector<int>    jets_nSV;
      vector<float>  jets_cef;
      vector<float>  jets_nef;
      vector<float>  jets_chf;
      vector<float>  jets_nhf;
      vector<float>  jets_phf;
      vector<float>  jets_medianLogIpSig;
      vector<int>    jets_missHits;
      vector<int>    jets_muonHits;
  };
}

void
emjet::OutputTree::Init() {
  // per event
  run   = 0;
  lumi  = 0;
  event = 0;
  bx    = 0;
  met_pt  = 0.;
  met_phi = 0.;
  // per jet
  jets_pt             .clear();
  jets_eta            .clear();
  jets_phi            .clear();
  jets_promptTracks   .clear();
  jets_dispTracks     .clear();
  jets_nSV            .clear();
  jets_cef            .clear();
  jets_nef            .clear();
  jets_chf            .clear();
  jets_nhf            .clear();
  jets_phf            .clear();
  jets_medianLogIpSig .clear();
  jets_missHits       .clear();
  jets_muonHits       .clear();
}

void
emjet::OutputTree::Branch(TTree* tree) {
  // per event
  tree->Branch("run"     , &run     ) ;
  tree->Branch("lumi"    , &lumi    ) ;
  tree->Branch("event"   , &event   ) ;
  tree->Branch("bx"      , &bx      ) ;
  tree->Branch("met_pt"  , &met_pt  ) ;
  tree->Branch("met_phi" , &met_phi ) ;
  // per jet
  tree->Branch("jets_pt"             , &jets_pt             ) ;
  tree->Branch("jets_eta"            , &jets_eta            ) ;
  tree->Branch("jets_phi"            , &jets_phi            ) ;
  tree->Branch("jets_promptTracks"   , &jets_promptTracks   ) ;
  tree->Branch("jets_dispTracks"     , &jets_dispTracks     ) ;
  tree->Branch("jets_nSV"            , &jets_nSV            ) ;
  tree->Branch("jets_cef"            , &jets_cef            ) ;
  tree->Branch("jets_nef"            , &jets_nef            ) ;
  tree->Branch("jets_chf"            , &jets_chf            ) ;
  tree->Branch("jets_nhf"            , &jets_nhf            ) ;
  tree->Branch("jets_phf"            , &jets_phf            ) ;
  tree->Branch("jets_medianLogIpSig" , &jets_medianLogIpSig ) ;
  tree->Branch("jets_missHits"       , &jets_missHits       ) ;
  tree->Branch("jets_muonHits"       , &jets_muonHits       ) ;
}

#endif
