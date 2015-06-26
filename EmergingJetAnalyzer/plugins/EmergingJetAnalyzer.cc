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
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

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

#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "TLorentzVector.h"
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

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::Service<TFileService> fs;

			edm::EDGetTokenT< reco::PFJetCollection > jetCollectionToken_;

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
      TH1F * h_jet_cef_trackFrac05;
      TH1F * h_jet_nef_trackFrac05;
      TH1F * h_jet_chf_trackFrac05;
      TH1F * h_jet_nhf_trackFrac05;

      TH1F * h_rVtx;
      TH1F * h_chi2Vtx;
      TH1F * h_nTracksVtx;
      TH1F * h_massVtx;
      TH1F * h_ptVtx;

      TH1F * h_jet_phf_trackFrac05;

      TH1F * h_jetIpTracks;
      TH1F * h_trackDrJet;
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


      TH1F * h_medianLogIpSig;
      TH1F * h_medianIp;
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
EmergingJetAnalyzer::EmergingJetAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

//    h_ptWtPixHits_vs_chargeFraction = fs->make<TH2F>( "h_ptWtPixHits_vs_chargeFraction"  , ";p_{T} wt. PIX hits;charge fraction", 10,  0., 4., 10, 0., 1. );
//    h_ptWtPixHits_vs_emFraction = fs->make<TH2F>( "h_ptWtPixHits_vs_emFraction"  , ";p_{T} wt. PIX hits;EM fraction", 10,  0., 4., 10, 0., 1. );

//    h_ptWtD0 = fs->make<TH1F>( "h_ptWtD0",";p_{T} weighted d0 [cm];" , 25, -3., 3.);
//    h_ptWtD0_trackFrac05 = fs->make<TH1F>( "h_ptWtD0_trackFrac05",";p_{T} weighted d0 [cm];" , 25, -3., 3.);
//    h_jetEtaPhi = fs->make<TH2F>( "h_jetEtaPhi" , ";jet #eta;jet #phi;" , 25, -2.5, 2.5, 25, -TMath::Pi(), TMath::Pi() );

    h_nVtxInEvent = fs->make<TH1F>( "h_nVtxInEvent",  ";nVtx;",101,-0.5,100.5);
    h_rVtx      = fs->make<TH1F>("h_rVtx",";vertex r [cm];",101,-0.5,100.5);
    h_chi2Vtx   = fs->make<TH1F>("h_chi2Vtx",";vertex #chi^{2};",16,-0.5,16);
    h_nTracksVtx= fs->make<TH1F>("h_nTracksVtx",";N tracks/vertex;",10,-0.5,9.5);
    h_massVtx   = fs->make<TH1F>("h_massVtx",";mass [GeV];",20,0.,10.);
    h_ptVtx     = fs->make<TH1F>("h_ptVtx",";pT [GeV];",50,0.,100.);

//    h_rVtxGenTracks = fs->make<TH1F>( "h_rVtxGenTracks", ";vertex r [cm];",300.,0.,150.);

    h_jetPtVsIndex = fs->make<TH2F>( "h_jetPtVsIndex", ";p_{T} [GeV];index", 100,0.,1000.,4,-0.5,3.5);
    h_jetIpTracks  = fs->make<TH1F>("h_jetIpTracks",";N disp trks;",20,-0.5,19.5);
    h_jetMatchedVertices  = fs->make<TH1F>("h_jetMatchedVertices",";N disp vtx;",20,-0.5,19.5);
    h_jetDispMuons  = fs->make<TH1F>("h_jetDispMuons",";N disp STA #mu;",20,-0.5,19.5);
    h_trackDrJet   = fs->make<TH1F>("h_trackDrJet",";#DeltaR(track, jet);",20.,0.,1.5);
    h_vertexDrJet   = fs->make<TH1F>("h_vertexDrJet",";#DeltaR(track, jet);",20.,0.,1.5);
    h_muonDrJet   = fs->make<TH1F>("h_muonDrJet",";#DeltaR(track, jet);",20.,0.,1.5);
    h_eventIpTracks  = fs->make<TH1F>("h_eventIpTracks",";N disp trks;",50,-0.5,99.5);
    h_eventIpTracks0_1  = fs->make<TH1F>("h_eventIpTracks0_1",";N disp trks;",50,-0.5,99.5);
    h_eventIpTracks1_0  = fs->make<TH1F>("h_eventIpTracks1_0",";N disp trks;",50,-0.5,99.5);
    h_eventIpTracks10_  = fs->make<TH1F>("h_eventIpTracks10_",";N disp trks;",50,-0.5,99.5);
    h_eventVertices  = fs->make<TH1F>("h_eventVertices",";N disp vtx;",20,-0.5,19.5);
    h_eventDispMuons  = fs->make<TH1F>("h_eventDispMuons",";N disp STA #mu;",20,-0.5,19.5);
    

    h_eventSumLogIp = fs->make<TH1F>( "h_eventSumLogIp", ";#Sigma(ln[ipSig]);",100,0.,200.);
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
    h_jetNeutralFraction =  fs->make<TH1F>("h_jetNeutralFraction",";neutral fraction;",20,0.,1.);

    h_sumEt      = fs->make<TH1F>("h_sumEt",";#Sigma(E_{T}) [GeV];",130,0.,13000.);

    h_logTrackDxy   = fs->make<TH1F>("h_logTrackDxy",";log(dxy);",100.,-3.,3.);
    h_logTrackDz    = fs->make<TH1F>("h_logTrackDz", ";log(dz);", 100.,-3.,3.);

    h_medianLogIpSig = fs->make<TH1F>("h_medianLogIpSig",";median log(ipSig);",100.,-3.,3.);
    h_medianIp       = fs->make<TH1F>("h_medianIp",      ";median IP;",        100.,-3.,3.);

//    h_jet_cef_trackFrac05     = fs->make<TH1F>("h_jet_cef_trackFrac05",";charged EM fraction;",20,0.,1.);
//    h_jet_nef_trackFrac05     = fs->make<TH1F>("h_jet_nef_trackFrac05",";neutral EM fraction;",20,0.,1.);
//    h_jet_chf_trackFrac05     = fs->make<TH1F>("h_jet_chf_trackFrac05",";charged had fraction;",20,0.,1.);
//    h_jet_nhf_trackFrac05     = fs->make<TH1F>("h_jet_nhf_trackFrac05",";neutral had fraction;",20,0.,1.);
//    h_jet_phf_trackFrac05     = fs->make<TH1F>("h_jet_phf_trackFrac05",";photon fraction;",20,0.,1.);

	jetCollectionToken_ = consumes< reco::PFJetCollection > (iConfig.getParameter<edm::InputTag>("srcJets"));
    
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

   edm::ESHandle<TransientTrackBuilder> theB;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

   edm::Handle<reco::PFJetCollection> pfjetH;
   iEvent.getByToken(jetCollectionToken_, pfjetH);

//   Handle<reco::BeamSpot> theBeamSpotHandle;
//   iEvent.getByLabel("offlineBeamSpot", theBeamSpotHandle);
//   const reco::BeamSpot* theBeamSpot = theBeamSpotHandle.product();

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

   edm::Handle<reco::PFMET> pfmet;
   iEvent.getByLabel("pfMet",pfmet);

   edm::Handle<reco::VertexCollection> primary_vertices;
   iEvent.getByLabel("offlinePrimaryVerticesWithBS",primary_vertices);
   const reco::Vertex& primary_vertex = primary_vertices->at(0);

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

   std::vector<reco::TransientTrack> ipTracks;
   std::vector<float> logIpSig;
   std::vector<float> ipVector;
   float ipSig = 3.;
   int ip0_1 = 0.;
   int ip1_0 = 0.;
   int ip10_ = 0.;
   float logSum = 0.;
   for (std::vector<reco::TransientTrack>::iterator itk = generalTracks.begin(); itk < generalTracks.end(); ++itk) {
       if (itk->track().pt() < 1.) continue;
       auto dxy_ipv = IPTools::absoluteTransverseImpactParameter(*itk, primary_vertex);
       h_logTrackDz->Fill(TMath::Log(fabs(itk->track().dz() - primary_vertex.position().z())));
       h_logTrackDxy->Fill(TMath::Log(dxy_ipv.second.value()));
//       std::cout << "Track with value, significance " << dxy_ipv.second.value() << "\t" << dxy_ipv.second.significance() << std::endl;
       if (dxy_ipv.second.value() < 0.05) continue;
       if (dxy_ipv.second.significance() < ipSig) continue;
       logIpSig.push_back(TMath::Log(dxy_ipv.second.significance()));
       ipVector.push_back(dxy_ipv.second.value());
       logSum += TMath::Log(dxy_ipv.second.significance());
//       if (itk->track().dxy(*theBeamSpot) < 0.05 && itk->track().dxy(*theBeamSpot) > -0.05) continue;
//       if (itk->track().dxy(*theBeamSpot) / itk->track().dxyError() < ipSig) continue;
//       std::cout << "\tq dxy " << itk->track().charge() << "\t" << itk->track().dxy(*theBeamSpot) << std::endl;
//       std::cout << "\tImpact parameter significance " << itk->track().dxy(*theBeamSpot) / itk->track().dxyError() << std::endl;

       ipTracks.push_back(*itk);
       if (dxy_ipv.second.value() > 0.1) ip0_1++; 
       if (dxy_ipv.second.value() > 1.0) ip1_0++; 
       if (dxy_ipv.second.value() > 10.) ip10_++; 
   }

   std::sort(logIpSig.begin(),logIpSig.end());
   std::sort(ipVector.begin(),   ipVector.end());
   if (logIpSig.size() != 0) {
        h_medianLogIpSig->Fill(logIpSig[logIpSig.size()/2]);
        h_medianIp->Fill(ipVector[ipVector.size()/2]);
   }

   h_eventIpTracks0_1->Fill(ip0_1);
   h_eventIpTracks1_0->Fill(ip1_0);
   h_eventIpTracks10_->Fill(ip10_);

   h_eventSumLogIp->Fill(logSum);

   // try my own vertex reco


   AdaptiveVertexReconstructor avr (2.0, 6.0, 0.5, true );
   std::vector<TransientVertex> theVertices = avr.vertices(ipTracks);
   std::vector<TLorentzVector>  vertexVectors;
   for (size_t ivtx = 0; ivtx < theVertices.size(); ++ivtx) {
        if (theVertices[ivtx].normalisedChiSquared() > 15.) continue;
        h_rVtx->Fill(theVertices[ivtx].position().perp());
        h_chi2Vtx->Fill(theVertices[ivtx].normalisedChiSquared());
        h_nTracksVtx->Fill(theVertices[ivtx].refittedTracks().size());
        TLorentzVector cand;
        // loop over the tracks
        for (std::vector<reco::TransientTrack>::const_iterator itk = theVertices[ivtx].refittedTracks().begin(); itk != theVertices[ivtx].refittedTracks().end(); ++itk) {
            TrajectoryStateClosestToPoint trajectory = itk->trajectoryStateClosestToPoint(theVertices[ivtx].position());
//            int charge = itk->track().charge();
            GlobalVector momentum = trajectory.momentum();
            TLorentzVector trackVector;  trackVector.SetPtEtaPhiM(momentum.perp(), momentum.eta(), momentum.phi(), 0.13957);
            cand += trackVector;
        }
        vertexVectors.push_back(cand);
        h_massVtx->Fill(cand.M());
        h_ptVtx->Fill(cand.Pt());
   }
   h_nVtxInEvent->Fill(theVertices.size());

//   return;
   
   int ijet = 0;
   for ( reco::PFJetCollection::const_iterator jet = pfjetH->begin(); jet != pfjetH->end(); ++jet ) {

       if (jet->pt() < 50.) continue;
       if (fabs(jet->eta()) > 2.5) continue;
       h_jetPtVsIndex->Fill(jet->pt(), ijet);
       TLorentzVector jetVector; jetVector.SetPtEtaPhiM(jet->pt(),jet->eta(),jet->phi(),0.);
       ++ijet;
//       std::cout << "Jet pt,eta,phi " << jet->pt() << "\t" << jet->eta() << "\t" << jet->phi() << std::endl;

       //std::cout << "PFJet pt,eta,phi" << "\t" << jet->pt() << "\t" << jet->eta() << "\t" << jet->phi() << std::endl;
       //std::cout << "\t neutralEM:  " << jet->neutralEmEnergyFraction() << std::endl;
       //std::cout << "\tneutralHad:  " << jet->neutralHadronEnergyFraction() << std::endl;
       float trackFrac = 0.;
       for (reco::TrackRefVector::iterator ijt = jet->getTrackRefs().begin(); ijt != jet->getTrackRefs().end(); ++ijt) {
           reco::TransientTrack itk = theB->build(*ijt);
           if (itk.track().pt() < 1.) continue;
           auto dxy_ipv = IPTools::absoluteTransverseImpactParameter(itk, primary_vertex);
           //       std::cout << "Track with value, significance " << dxy_ipv.second.value() << "\t" << dxy_ipv.second.significance() << std::endl;
           if (dxy_ipv.second.value() > 0.05) continue;
           if (dxy_ipv.second.significance() > ipSig) continue;

           trackFrac += itk.track().pt();
       }
       trackFrac /= jet->pt();
       h_trackFrac->Fill(trackFrac);

       // now try to match the displaced tracks to this jet
       int matchedTracks = 0;
//       for (std::vector<reco::TransientTrack>::iterator itk = ipTracks.begin(); itk < ipTracks.end(); ++itk) {

//           if (!itk->innermostMeasurementState().isValid()) continue;
//           TLorentzVector trackVector; trackVector.SetPtEtaPhiM(itk->innermostMeasurementState().globalPosition().perp(), itk->innermostMeasurementState().globalPosition().eta(), itk->innermostMeasurementState().globalPosition().phi(), 0.13957);
//           h_trackDrJet->Fill(trackVector.DeltaR(jetVector));
//           if (trackVector.DeltaR(jetVector) < 0.5) ++matchedTracks;

//       }
       h_jetIpTracks->Fill(matchedTracks);

       // same with the vertices
       int matchedVertices = 0;
       for (std::vector<TransientVertex>::iterator ivtx = theVertices.begin(); ivtx != theVertices.end(); ++ivtx) {
            TLorentzVector vertexPosition;  vertexPosition.SetPtEtaPhiM(ivtx->position().perp(),ivtx->position().eta(),ivtx->position().phi(),0.);
            h_vertexDrJet->Fill(vertexPosition.DeltaR(jetVector));
            if (vertexPosition.DeltaR(jetVector) < 0.5) ++matchedVertices;
       }
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
       h_jet_nhf->Fill(jet->neutralHadronEnergyFraction());
       h_jet_phf->Fill(jet->photonEnergyFraction());

       h_jetNeutralFraction->Fill(jet->neutralHadronEnergyFraction() + jet->neutralEmEnergyFraction());

   }

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


   h_eventIpTracks->Fill(ipTracks.size());
   h_eventVertices->Fill(theVertices.size());
   h_eventDispMuons->Fill(dispMuons);
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
