import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")


process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_FULL', '')

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/group/phys_exotica/EmergingJets/EmergingJets_ModelD_TuneCUETP8M1_13TeV_pythia8Mod/RECO/150525_131119/0000/step2_1.root',
        '/store/group/phys_exotica/EmergingJets/EmergingJets_ModelD_TuneCUETP8M1_13TeV_pythia8Mod/RECO/150525_131119/0000/step2_10.root',
        '/store/group/phys_exotica/EmergingJets/EmergingJets_ModelD_TuneCUETP8M1_13TeV_pythia8Mod/RECO/150525_131119/0000/step2_2.root',
        '/store/group/phys_exotica/EmergingJets/EmergingJets_ModelD_TuneCUETP8M1_13TeV_pythia8Mod/RECO/150525_131119/0000/step2_3.root',
        '/store/group/phys_exotica/EmergingJets/EmergingJets_ModelD_TuneCUETP8M1_13TeV_pythia8Mod/RECO/150525_131119/0000/step2_4.root',
        '/store/group/phys_exotica/EmergingJets/EmergingJets_ModelD_TuneCUETP8M1_13TeV_pythia8Mod/RECO/150525_131119/0000/step2_5.root',
        '/store/group/phys_exotica/EmergingJets/EmergingJets_ModelD_TuneCUETP8M1_13TeV_pythia8Mod/RECO/150525_131119/0000/step2_6.root',
        '/store/group/phys_exotica/EmergingJets/EmergingJets_ModelD_TuneCUETP8M1_13TeV_pythia8Mod/RECO/150525_131119/0000/step2_7.root',
        '/store/group/phys_exotica/EmergingJets/EmergingJets_ModelD_TuneCUETP8M1_13TeV_pythia8Mod/RECO/150525_131119/0000/step2_8.root',
        '/store/group/phys_exotica/EmergingJets/EmergingJets_ModelD_TuneCUETP8M1_13TeV_pythia8Mod/RECO/150525_131119/0000/step2_9.root'
    )
)

process.demo = cms.EDAnalyzer('EmergingJetAnalyzer'
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("histo_D.root") )

process.p = cms.Path(process.demo)
