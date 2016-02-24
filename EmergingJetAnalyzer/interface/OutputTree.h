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
      vector<float>  jets_cef;
      vector<float>  jets_nef;
      vector<float>  jets_chf;
      vector<float>  jets_nhf;
      vector<float>  jets_phf;
      vector<int>    jets_nPromptTracks;
      vector<int>    jets_nDispTracks;
      vector<int>    jets_nSV;
      vector<float>  jets_medianLogIpSig;
      vector<int>    jets_missHits;
      vector<int>    jets_muonHits;
      vector<float>  jets_alphaMax;
      vector<int>    jets_nDarkPions;
      vector<float>  jets_minDRDarkPion;
      // per track (for a given jet)
      vector< vector<int> >   tracks_nHits;
      vector< vector<int> >   tracks_nMissInnerHits;
      vector< vector<float> > tracks_ipXY;
      vector< vector<float> > tracks_ipZ;
      vector< vector<float> > tracks_ipXYSig;
      // per secondary vertex
      vector<float>  vertex_x;
      vector<float>  vertex_y;
      vector<float>  vertex_z;
      vector<float>  vertex_xError;
      vector<float>  vertex_yError;
      vector<float>  vertex_zError;
      vector<float>  vertex_Lxy;
      vector<float>  vertex_mass;
      vector<float>  vertex_chi2;
      vector<float>  vertex_ndof;
      vector<float>  vertex_pt2sum;
      // vector< vector<float> > vector_tracks_weight;
      // vector< vector<float> > vector_tracks_pt;
      // vector< vector<float> > vector_tracks_eta;
      // vector< vector<float> > vector_tracks_phi;
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
 jets_pt               .clear();
 jets_eta              .clear();
 jets_phi              .clear();
 jets_cef              .clear();
 jets_nef              .clear();
 jets_chf              .clear();
 jets_nhf              .clear();
 jets_phf              .clear();
 jets_nPromptTracks    .clear();
 jets_nDispTracks      .clear();
 jets_nSV              .clear();
 jets_medianLogIpSig   .clear();
 jets_missHits         .clear();
 jets_muonHits         .clear();
 jets_alphaMax         .clear();
 jets_nDarkPions       .clear();
 jets_minDRDarkPion    .clear();
 // per track
 tracks_nHits          .clear();
 tracks_nMissInnerHits .clear();
 tracks_ipXY           .clear();
 tracks_ipZ            .clear();
 tracks_ipXYSig        .clear();
 // per secondary vertex
 vertex_x              .clear();
 vertex_y              .clear();
 vertex_z              .clear();
 vertex_xError         .clear();
 vertex_yError         .clear();
 vertex_zError         .clear();
 vertex_Lxy            .clear();
 vertex_mass           .clear();
 vertex_chi2           .clear();
 vertex_ndof           .clear();
 vertex_pt2sum         .clear();
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
  tree->Branch("jets_cef"            , &jets_cef            ) ;
  tree->Branch("jets_nef"            , &jets_nef            ) ;
  tree->Branch("jets_chf"            , &jets_chf            ) ;
  tree->Branch("jets_nhf"            , &jets_nhf            ) ;
  tree->Branch("jets_phf"            , &jets_phf            ) ;
  tree->Branch("jets_nPromptTracks"  , &jets_nPromptTracks  ) ;
  tree->Branch("jets_nDispTracks"    , &jets_nDispTracks    ) ;
  tree->Branch("jets_nSV"            , &jets_nSV            ) ;
  tree->Branch("jets_medianLogIpSig" , &jets_medianLogIpSig ) ;
  tree->Branch("jets_missHits"       , &jets_missHits       ) ;
  tree->Branch("jets_muonHits"       , &jets_muonHits       ) ;
  tree->Branch("jets_alphaMax"       , &jets_alphaMax       ) ;
  tree->Branch("jets_nDarkPions"     , &jets_nDarkPions     ) ;
  tree->Branch("jets_minDRDarkPion"  , &jets_minDRDarkPion  ) ;
  // per track
  tree->Branch("tracks_nHits"          , &tracks_nHits          ) ;
  tree->Branch("tracks_nMissInnerHits" , &tracks_nMissInnerHits ) ;
  tree->Branch("tracks_ipXY"           , &tracks_ipXY           ) ;
  tree->Branch("tracks_ipZ"            , &tracks_ipZ            ) ;
  tree->Branch("tracks_ipXYSig"        , &tracks_ipXYSig        ) ;
  // per secondary vertex
  tree->Branch("vertex_x"      ,   &vertex_x      ) ;
  tree->Branch("vertex_y"      ,   &vertex_y      ) ;
  tree->Branch("vertex_z"      ,   &vertex_z      ) ;
  tree->Branch("vertex_xError" ,   &vertex_xError ) ;
  tree->Branch("vertex_yError" ,   &vertex_yError ) ;
  tree->Branch("vertex_zError" ,   &vertex_zError ) ;
  tree->Branch("vertex_Lxy"    ,   &vertex_Lxy    ) ;
  tree->Branch("vertex_mass"   ,   &vertex_mass   ) ;
  tree->Branch("vertex_chi2"   ,   &vertex_chi2   ) ;
  tree->Branch("vertex_ndof"   ,   &vertex_ndof   ) ;
  tree->Branch("vertex_pt2sum" ,   &vertex_pt2sum ) ;
}

#endif
