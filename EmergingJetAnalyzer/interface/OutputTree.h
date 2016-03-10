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
      vector< vector<int> >   tracks_algo;
      vector< vector<int> >   tracks_originalAlgo;
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
 tracks_algo           .clear();
 tracks_originalAlgo   .clear();
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
  #define BRANCH(tree, branch) (tree)->Branch(#branch, &branch);
  // per event
  BRANCH(tree, run);
  BRANCH(tree, lumi);
  BRANCH(tree, event);
  BRANCH(tree, bx);
  BRANCH(tree, met_pt);
  BRANCH(tree, met_phi);
  // per jet
  BRANCH(tree, jets_pt);
  BRANCH(tree, jets_eta);
  BRANCH(tree, jets_phi);
  BRANCH(tree, jets_cef);
  BRANCH(tree, jets_nef);
  BRANCH(tree, jets_chf);
  BRANCH(tree, jets_nhf);
  BRANCH(tree, jets_phf);
  BRANCH(tree, jets_nPromptTracks);
  BRANCH(tree, jets_nDispTracks);
  BRANCH(tree, jets_nSV);
  BRANCH(tree, jets_medianLogIpSig);
  BRANCH(tree, jets_missHits);
  BRANCH(tree, jets_muonHits);
  BRANCH(tree, jets_alphaMax);
  BRANCH(tree, jets_nDarkPions);
  BRANCH(tree, jets_minDRDarkPion);
  // per track
  BRANCH(tree, tracks_algo);
  BRANCH(tree, tracks_originalAlgo);
  BRANCH(tree, tracks_nHits);
  BRANCH(tree, tracks_nMissInnerHits);
  BRANCH(tree, tracks_ipXY);
  BRANCH(tree, tracks_ipZ);
  BRANCH(tree, tracks_ipXYSig);
  // per secondary vertex
  BRANCH(tree, vertex_x);
  BRANCH(tree, vertex_y);
  BRANCH(tree, vertex_z);
  BRANCH(tree, vertex_xError);
  BRANCH(tree, vertex_yError);
  BRANCH(tree, vertex_zError);
  BRANCH(tree, vertex_Lxy);
  BRANCH(tree, vertex_mass);
  BRANCH(tree, vertex_chi2);
  BRANCH(tree, vertex_ndof);
  BRANCH(tree, vertex_pt2sum);
}

#endif
