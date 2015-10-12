############################################################
# Important things to change before running
# - Default value for inputFiles
# - Global Tag (between releases)
# - datasets.py import location
############################################################

import FWCore.ParameterSet.Config as cms
# Command line argument parsing
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process('SKIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')


## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

## Options and Output Report
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False),
        # SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
options.maxEvents = -1 # -1 means all events
# options.inputFiles= 'root://xrootd-cms.infn.it//store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00CC714A-F86B-E411-B99A-0025904B5FB8.root'
# options.inputFiles= 'file:/afs/cern.ch/work/y/yoshin/public/RunIISpring15DR74/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/00CC714A-F86B-E411-B99A-0025904B5FB8.root'
# options.inputFiles= '/store/group/phys_exotica/EmergingJets/EmergingJets_ModelA_TuneCUETP8M1_13TeV_pythia8Mod/RECO/150715_195547/0000/aodsim_10.root'
options.inputFiles= '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/96EF1A5F-8115-E511-AF17-02163E0125CE.root'
options.outputFile = ''
options.register ('crab',
                  0, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Set to 1 to run on CRAB.")
options.crab = 0
options.register ('data',
                  0, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Set to 1 for data.")
options.data = 0


# get and parse the command line arguments
options.parseArguments()


# Unscheduled mode
# process.options.allowUnscheduled = cms.untracked.bool(False)
process.options.allowUnscheduled = cms.untracked.bool(True)



### =========== Global configuration ==========================

## MET flavor configuration -> may be time consuming
## NoPu and MVA met need to be generated with some refeernce objects (e;g. leptons)
## look for the corresponding area in the config file to set your own definition


### ===========================================================
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V9', '')

process.source = cms.Source("PoolSource",
		# fileNames = cms.untracked.vstring("file:input.root"),
    skipEvents = cms.untracked.uint32(0)
)
# if options.inputFiles :
process.source.fileNames = cms.untracked.vstring(options.inputFiles)

if options.data and not options.crab:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = 'goodList.json').getVLuminosityBlockRange()

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
    )

##
## To remove the "begin job processing line " from printing
##
#process.MessageLogger = cms.Service("MessageLogger",
#   destinations = cms.untracked.vstring('cout',),
#   cout         = cms.untracked.PSet(threshold = cms.untracked.string('ERROR'))
#)
# process.MessageLogger.cout.INFO.limit = 0
# process.MessageLogger.cerr.threshold = "DEBUG"
# process.MessageLogger.categories.append('JetFilter')

# process.MessageLogger.cerr.FwkReport.reportEvery = 1
# process.MessageLogger.cerr.default.limit = 100
#process.options.wantSummary = False

############################################################
# PAT Default Configuration
############################################################
# process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
# process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
# process.load("PhysicsTools.PatAlgos.cleaningLayer1.cleanPatCandidates_cff")
# process.load("PhysicsTools.PatAlgos.selectionLayer1.countPatCandidates_cff")

# from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning

############################################################
# Recreate miniAOD
############################################################
if options.data: 
    # customisation of the process.
    process.load('Configuration.StandardSequences.PAT_cff')
    # Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
    from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllData 
    #call to customisation function miniAOD_customizeAllData imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
    process = miniAOD_customizeAllData(process)
else: 
    # customisation of the process.
    process.load('Configuration.StandardSequences.PATMC_cff')
    # Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
    from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 
    #call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
    process = miniAOD_customizeAllMC(process)
############################################################

############################################################
# Private modules
############################################################
process.zJetFilter = cms.EDFilter("ZJetFilter",
    srcMuons = cms.InputTag("slimmedMuons"),
    srcElectrons = cms.InputTag("slimmedElectrons"),
    srcJets = cms.InputTag("slimmedJets"),
    minPtMuon = cms.double(20.0),
    minPtElectron = cms.double(30.0),
    minZMass = cms.double(80.0),
    maxZMass = cms.double(100.0),
    maxDeltaPhi = cms.double(0.4),
    minPtSelectedJet = cms.double(50.0),
    maxPtAdditionalJets = cms.double(50.0),
)
# process.jetFilter = cms.EDFilter("JetFilter",
#     srcJets = cms.InputTag("ak4PFJetsCHS"),
#     # srcJets = cms.InputTag("patJets"),
#     # additionalCut = cms.string(""),
#     additionalCut = cms.string("abs(eta) < 2.5 && pt > 50.0"),
#     jetCuts = cms.VPSet(
#         cms.PSet(
#             minPt = cms.double(400.0),
#             maxEta = cms.double(2.5),
#             stringCut = cms.string(""),
#         ),
#         cms.PSet(
#             minPt = cms.double(200.0),
#             maxEta = cms.double(2.5),
#             stringCut = cms.string(""),
#         ),
#         cms.PSet(
#             minPt = cms.double(125.0),
#             maxEta = cms.double(2.5),
#             stringCut = cms.string(""),
#         ),
#         cms.PSet(
#             minPt = cms.double(50.0),
#             maxEta = cms.double(2.5),
#             stringCut = cms.string(""),
#         ),
#     )
# )

# process.demo = cms.EDAnalyzer('EmergingJetAnalyzer'
# )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("histo.root")
)

# process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
# process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
#                                    src = cms.InputTag("genParticles"),                                                                 
#                                    printP4 = cms.untracked.bool(False),
#                                    printPtEtaPhi = cms.untracked.bool(True),
#                                    printVertex = cms.untracked.bool(False),
#                                    printStatus = cms.untracked.bool(True),
#                                    printIndex = cms.untracked.bool(False),
#                                    # status = cms.untracked.vint32( 3 )
#                                    )
# process.printDecay = cms.EDAnalyzer("ParticleDecayDrawer",
#     src = cms.InputTag("genParticles"),
#     printP4 = cms.untracked.bool(True),
#     printPtEtaPhi = cms.untracked.bool(True),
#     printVertex = cms.untracked.bool(False)
#   )
#
process.p = cms.Path( process.zJetFilter )



############################################################
# PATH definition
############################################################
# process.muonPath = cms.Path( process.countPatMuons )
# process.electronPath = cms.Path( process.countPatElectrons )

# process.jetPath = cms.Path( process.jetFilter )




process.out = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('AODSIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fileName = cms.untracked.string('file:aodsim.root'),
    outputCommands = AODSIMEventContent.outputCommands,
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    )
)
if options.outputFile=='.root' or options.outputFile.find('_numEvent')==0:
    print """
############################################################
Warning: outputFile unspecified. Writing to aodsim.root
############################################################
"""
else:
    process.out.fileName = cms.untracked.string(options.outputFile)
    print """
############################################################
Writing to """ + options.outputFile + """
############################################################
"""


process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('cout'),
    categories = cms.untracked.vstring('ZJetFilter'),
    debugModules = cms.untracked.vstring('zJetFilter'),
    cout         = cms.untracked.PSet(
        threshold = cms.untracked.string('DEBUG'),
        # DEBUG = cms.untracked.int32(0),
        WARNING = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
        ),
    )
)
# storage
process.outpath = cms.EndPath(process.out) #dummy


