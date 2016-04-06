#ifndef EmergingJetAnalysis_EmergingJetAnalyzer_EmJetEvent_h
#define EmergingJetAnalysis_EmergingJetAnalyzer_EmJetEvent_h

// Class describing the content of output for EmergingJetAnalyzer

#include <vector>
#include <functional>
#include "EmergingJetAnalysis/EmergingJetAnalyzer/interface/OutputTree.h"
#define DEFAULTVALUE -1

using std::vector;

namespace emjet
{
  class Track;
  class Vertex;
  class Jet;
  class Event {
  public:
    Event() {
      //[[[cog
      //template_string = "$name = DEFAULTVALUE;"
      //import vars_EmJetAnalyzer as m
      //for vardict in m.event_vardicts: m.replaceSingleLine(template_string, vardict)
      //]]]
      //[[[end]]]

      jet_vector.clear();
    };
    //[[[cog
    //template_string = "$type $name;"
    //import vars_EmJetAnalyzer as m
    //for vardict in m.event_vardicts: m.replaceSingleLine(template_string, vardict)
    //]]]
    //[[[end]]]

    vector<Jet> jet_vector;
  };
  class Jet {
  public:
    Jet(){
      //[[[cog
      //template_string = "$name = DEFAULTVALUE;"
      //import vars_EmJetAnalyzer as m
      //for vardict in m.jet_vardicts: m.replaceSingleLine(template_string, vardict)
      //]]]
      //[[[end]]]

      track_vector.clear();
      vertex_vector.clear();
    }
    //[[[cog
    //template_string = "$type $name;"
    //import vars_EmJetAnalyzer as m
    //for vardict in m.jet_vardicts: m.replaceSingleLine(template_string, vardict)
    //]]]
    //[[[end]]]
    vector<Track>    track_vector;
    vector<Vertex>   vertex_vector;
  };
  class Track {
  public:
    Track() {
      //[[[cog
      //template_string = "$name = DEFAULTVALUE;"
      //import vars_EmJetAnalyzer as m
      //for vardict in m.jet_track_vardicts: m.replaceSingleLine(template_string, vardict)
      //]]]
      //[[[end]]]
    }
    //[[[cog
    //template_string = "$type $name;"
    //import vars_EmJetAnalyzer as m
    //for vardict in m.jet_track_vardicts: m.replaceSingleLine(template_string, vardict)
    //]]]
    //[[[end]]]
  };
  class Vertex {
  public:
    Vertex() {
      //[[[cog
      //template_string = "$name = DEFAULTVALUE;"
      //import vars_EmJetAnalyzer as m
      //for vardict in m.jet_vertex_vardicts: m.replaceSingleLine(template_string, vardict)
      //]]]
      //[[[end]]]
    }
    //[[[cog
    //template_string = "$type $name;"
    //import vars_EmJetAnalyzer as m
    //for vardict in m.jet_vertex_vardicts: m.replaceSingleLine(template_string, vardict)
    //]]]
    //[[[end]]]
  };
}

// Turn vector of objects, into vector of member variable by calling func(object)
// Output is placed in provided vector reference
template <typename Object, typename T>
void
vectorize(const vector<Object>& input, std::function<T (const Object &)> func, vector<T>& output)
{
  output.clear();
  output.reserve(input.size()); // Doesn't reduce capacity
  for (const auto& obj : input) {
    output.push_back( func(obj) );
  }
}

// Version of vectorize() that returns a vector object
template <typename Object, typename T>
vector<T>
vectorize_new(const vector<Object>& input, std::function<T (const Object &)> func)
{
  vector<T> output;
  vectorize(input, func, output);
  return output;
}

using emjet::Track;
using emjet::Vertex;
using emjet::Jet;
using emjet::Event;

void
WriteEventToOutput(Event event, emjet::OutputTree* otree)
{
  otree->Init(); // Reset all values and clear all vectors
  // Event-level variables, e.g. int, float, etc.
  {
    //[[[cog
    //template_string = "otree->$name = event.$name;"
    //import vars_EmJetAnalyzer as m
    //for vardict in m.event_vardicts: m.replaceSingleLine(template_string, vardict)
    //]]]
    otree->run                  = event.run                 ;
    otree->lumi                 = event.lumi                ;
    otree->event                = event.event               ;
    otree->bx                   = event.bx                  ;
    otree->nVtx                 = event.nVtx                ;
    otree->nGoodVtx             = event.nGoodVtx            ;
    otree->nTrueInt             = event.nTrueInt            ;
    otree->met_pt               = event.met_pt              ;
    otree->met_phi              = event.met_phi             ;
    //[[[end]]]
    otree->run      = event.run      ;
  }
  // Jet-level variables, e.g. vector<int>, vector<float>, etc.
  {
    //[[[cog
    //template_string = "vectorize<Jet, $cpptype>(event.jet_vector, [](const emjet::Jet& obj ){return obj.$name;}, otree->jet_$name);"
    //import vars_EmJetAnalyzer as m
    //for vardict in m.jet_vardicts: m.replaceSingleLine(template_string, vardict)
    //]]]
    vectorize<Jet, int   >(event.jet_vector, [](const emjet::Jet& obj ){return obj.index               ;}, otree->jet_index               );
    vectorize<Jet, int   >(event.jet_vector, [](const emjet::Jet& obj ){return obj.source              ;}, otree->jet_source              );
    vectorize<Jet, float >(event.jet_vector, [](const emjet::Jet& obj ){return obj.pt                  ;}, otree->jet_pt                  );
    vectorize<Jet, float >(event.jet_vector, [](const emjet::Jet& obj ){return obj.eta                 ;}, otree->jet_eta                 );
    vectorize<Jet, float >(event.jet_vector, [](const emjet::Jet& obj ){return obj.phi                 ;}, otree->jet_phi                 );
    vectorize<Jet, float >(event.jet_vector, [](const emjet::Jet& obj ){return obj.cef                 ;}, otree->jet_cef                 );
    vectorize<Jet, float >(event.jet_vector, [](const emjet::Jet& obj ){return obj.nef                 ;}, otree->jet_nef                 );
    vectorize<Jet, float >(event.jet_vector, [](const emjet::Jet& obj ){return obj.chf                 ;}, otree->jet_chf                 );
    vectorize<Jet, float >(event.jet_vector, [](const emjet::Jet& obj ){return obj.nhf                 ;}, otree->jet_nhf                 );
    vectorize<Jet, float >(event.jet_vector, [](const emjet::Jet& obj ){return obj.phf                 ;}, otree->jet_phf                 );
    vectorize<Jet, int   >(event.jet_vector, [](const emjet::Jet& obj ){return obj.nPromptTracks       ;}, otree->jet_nPromptTracks       );
    vectorize<Jet, int   >(event.jet_vector, [](const emjet::Jet& obj ){return obj.nDispTracks         ;}, otree->jet_nDispTracks         );
    vectorize<Jet, int   >(event.jet_vector, [](const emjet::Jet& obj ){return obj.nSV                 ;}, otree->jet_nSV                 );
    vectorize<Jet, float >(event.jet_vector, [](const emjet::Jet& obj ){return obj.medianLogIpSig      ;}, otree->jet_medianLogIpSig      );
    vectorize<Jet, int   >(event.jet_vector, [](const emjet::Jet& obj ){return obj.missHits            ;}, otree->jet_missHits            );
    vectorize<Jet, int   >(event.jet_vector, [](const emjet::Jet& obj ){return obj.muonHits            ;}, otree->jet_muonHits            );
    vectorize<Jet, float >(event.jet_vector, [](const emjet::Jet& obj ){return obj.alphaMax            ;}, otree->jet_alphaMax            );
    vectorize<Jet, int   >(event.jet_vector, [](const emjet::Jet& obj ){return obj.nDarkPions          ;}, otree->jet_nDarkPions          );
    vectorize<Jet, float >(event.jet_vector, [](const emjet::Jet& obj ){return obj.minDRDarkPion       ;}, otree->jet_minDRDarkPion       );
    //[[[end]]]
  }
  // Jet-Track-level variables
  {
    for (const auto jet : event.jet_vector) {
      //[[[cog
      //import vars_EmJetAnalyzer as m
      //template_string = "auto $name = vectorize_new<Track, $cpptype>(jet.track_vector, [](const Track& obj ){return obj.$name;}); otree->track_$name.push_back($name);"
      //for vardict in m.jet_track_vardicts: m.replaceSingleLine(template_string, vardict)
      //]]]
      auto index                = vectorize_new<Track, int   >(jet.track_vector, [](const Track& obj ){return obj.index               ;}); otree->track_index               .push_back(index               );
      auto source               = vectorize_new<Track, int   >(jet.track_vector, [](const Track& obj ){return obj.source              ;}); otree->track_source              .push_back(source              );
      auto pt                   = vectorize_new<Track, float >(jet.track_vector, [](const Track& obj ){return obj.pt                  ;}); otree->track_pt                  .push_back(pt                  );
      auto eta                  = vectorize_new<Track, float >(jet.track_vector, [](const Track& obj ){return obj.eta                 ;}); otree->track_eta                 .push_back(eta                 );
      auto phi                  = vectorize_new<Track, float >(jet.track_vector, [](const Track& obj ){return obj.phi                 ;}); otree->track_phi                 .push_back(phi                 );
      auto pca_r                = vectorize_new<Track, float >(jet.track_vector, [](const Track& obj ){return obj.pca_r               ;}); otree->track_pca_r               .push_back(pca_r               );
      auto pca_eta              = vectorize_new<Track, float >(jet.track_vector, [](const Track& obj ){return obj.pca_eta             ;}); otree->track_pca_eta             .push_back(pca_eta             );
      auto pca_phi              = vectorize_new<Track, float >(jet.track_vector, [](const Track& obj ){return obj.pca_phi             ;}); otree->track_pca_phi             .push_back(pca_phi             );
      auto algo                 = vectorize_new<Track, int   >(jet.track_vector, [](const Track& obj ){return obj.algo                ;}); otree->track_algo                .push_back(algo                );
      auto originalAlgo         = vectorize_new<Track, int   >(jet.track_vector, [](const Track& obj ){return obj.originalAlgo        ;}); otree->track_originalAlgo        .push_back(originalAlgo        );
      auto nHits                = vectorize_new<Track, int   >(jet.track_vector, [](const Track& obj ){return obj.nHits               ;}); otree->track_nHits               .push_back(nHits               );
      auto nMissInnerHits       = vectorize_new<Track, int   >(jet.track_vector, [](const Track& obj ){return obj.nMissInnerHits      ;}); otree->track_nMissInnerHits      .push_back(nMissInnerHits      );
      auto nTrkLayers           = vectorize_new<Track, int   >(jet.track_vector, [](const Track& obj ){return obj.nTrkLayers          ;}); otree->track_nTrkLayers          .push_back(nTrkLayers          );
      auto nMissInnerTrkLayers  = vectorize_new<Track, int   >(jet.track_vector, [](const Track& obj ){return obj.nMissInnerTrkLayers ;}); otree->track_nMissInnerTrkLayers .push_back(nMissInnerTrkLayers );
      auto nMissOuterTrkLayers  = vectorize_new<Track, int   >(jet.track_vector, [](const Track& obj ){return obj.nMissOuterTrkLayers ;}); otree->track_nMissOuterTrkLayers .push_back(nMissOuterTrkLayers );
      auto nMissTrkLayers       = vectorize_new<Track, int   >(jet.track_vector, [](const Track& obj ){return obj.nMissTrkLayers      ;}); otree->track_nMissTrkLayers      .push_back(nMissTrkLayers      );
      auto nPxlLayers           = vectorize_new<Track, int   >(jet.track_vector, [](const Track& obj ){return obj.nPxlLayers          ;}); otree->track_nPxlLayers          .push_back(nPxlLayers          );
      auto nMissInnerPxlLayers  = vectorize_new<Track, int   >(jet.track_vector, [](const Track& obj ){return obj.nMissInnerPxlLayers ;}); otree->track_nMissInnerPxlLayers .push_back(nMissInnerPxlLayers );
      auto nMissOuterPxlLayers  = vectorize_new<Track, int   >(jet.track_vector, [](const Track& obj ){return obj.nMissOuterPxlLayers ;}); otree->track_nMissOuterPxlLayers .push_back(nMissOuterPxlLayers );
      auto nMissPxlLayers       = vectorize_new<Track, int   >(jet.track_vector, [](const Track& obj ){return obj.nMissPxlLayers      ;}); otree->track_nMissPxlLayers      .push_back(nMissPxlLayers      );
      auto ipXY                 = vectorize_new<Track, float >(jet.track_vector, [](const Track& obj ){return obj.ipXY                ;}); otree->track_ipXY                .push_back(ipXY                );
      auto ipZ                  = vectorize_new<Track, float >(jet.track_vector, [](const Track& obj ){return obj.ipZ                 ;}); otree->track_ipZ                 .push_back(ipZ                 );
      auto ipXYSig              = vectorize_new<Track, float >(jet.track_vector, [](const Track& obj ){return obj.ipXYSig             ;}); otree->track_ipXYSig             .push_back(ipXYSig             );
      auto dRToJetAxis          = vectorize_new<Track, float >(jet.track_vector, [](const Track& obj ){return obj.dRToJetAxis         ;}); otree->track_dRToJetAxis         .push_back(dRToJetAxis         );
      auto distanceToJet        = vectorize_new<Track, float >(jet.track_vector, [](const Track& obj ){return obj.distanceToJet       ;}); otree->track_distanceToJet       .push_back(distanceToJet       );
      //[[[end]]]
    }
  }
  // Jet-Vertex-level variables
  {
    for (const auto jet : event.jet_vector) {
      //[[[cog
      //import vars_EmJetAnalyzer as m
      //template_string = "auto $name = vectorize_new<Vertex, $cpptype>(jet.vertex_vector, [](const Vertex& obj ){return obj.$name;}); otree->vertex_$name.push_back($name);"
      //for vardict in m.jet_vertex_vardicts: m.replaceSingleLine(template_string, vardict)
      //]]]
      auto index                = vectorize_new<Vertex, int   >(jet.vertex_vector, [](const Vertex& obj ){return obj.index               ;}); otree->vertex_index               .push_back(index               );
      auto source               = vectorize_new<Vertex, int   >(jet.vertex_vector, [](const Vertex& obj ){return obj.source              ;}); otree->vertex_source              .push_back(source              );
      auto x                    = vectorize_new<Vertex, float >(jet.vertex_vector, [](const Vertex& obj ){return obj.x                   ;}); otree->vertex_x                   .push_back(x                   );
      auto y                    = vectorize_new<Vertex, float >(jet.vertex_vector, [](const Vertex& obj ){return obj.y                   ;}); otree->vertex_y                   .push_back(y                   );
      auto z                    = vectorize_new<Vertex, float >(jet.vertex_vector, [](const Vertex& obj ){return obj.z                   ;}); otree->vertex_z                   .push_back(z                   );
      auto xError               = vectorize_new<Vertex, float >(jet.vertex_vector, [](const Vertex& obj ){return obj.xError              ;}); otree->vertex_xError              .push_back(xError              );
      auto yError               = vectorize_new<Vertex, float >(jet.vertex_vector, [](const Vertex& obj ){return obj.yError              ;}); otree->vertex_yError              .push_back(yError              );
      auto zError               = vectorize_new<Vertex, float >(jet.vertex_vector, [](const Vertex& obj ){return obj.zError              ;}); otree->vertex_zError              .push_back(zError              );
      auto deltaR               = vectorize_new<Vertex, float >(jet.vertex_vector, [](const Vertex& obj ){return obj.deltaR              ;}); otree->vertex_deltaR              .push_back(deltaR              );
      auto Lxy                  = vectorize_new<Vertex, float >(jet.vertex_vector, [](const Vertex& obj ){return obj.Lxy                 ;}); otree->vertex_Lxy                 .push_back(Lxy                 );
      auto mass                 = vectorize_new<Vertex, float >(jet.vertex_vector, [](const Vertex& obj ){return obj.mass                ;}); otree->vertex_mass                .push_back(mass                );
      auto chi2                 = vectorize_new<Vertex, float >(jet.vertex_vector, [](const Vertex& obj ){return obj.chi2                ;}); otree->vertex_chi2                .push_back(chi2                );
      auto ndof                 = vectorize_new<Vertex, float >(jet.vertex_vector, [](const Vertex& obj ){return obj.ndof                ;}); otree->vertex_ndof                .push_back(ndof                );
      auto pt2sum               = vectorize_new<Vertex, float >(jet.vertex_vector, [](const Vertex& obj ){return obj.pt2sum              ;}); otree->vertex_pt2sum              .push_back(pt2sum              );
      //[[[end]]]
    }
  }
}

// write_jet_to_event(emjet::Jet jet, emjet::OutputTree* otree)
// {
//   otree.jet_pt.push_back(jet.pt);
//   track_pt = vectorize(jet.track_vector, auto[](auto track){return track.pt;} )
// }

#endif
