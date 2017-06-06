import FWCore.ParameterSet.Config as cms

########################################
# Command line argument parsing
########################################
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')
options.register ('crab',
                  0, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Set to 1 to run on CRAB.")
options.register ('data',
                  0, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Set to 1 for data.")
sample_options = ['signal', 'background', 'wjet', 'gjet'] # Valid options for sample
options.register ('sample',
                  'signal', # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "Specify type of sample. Valid values: %s" % sample_options)
steps_options = ['skim', 'analyze'] # Valid options for steps
options.register ('steps',
                  [],
                  VarParsing.VarParsing.multiplicity.list, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "Steps to execute. Possible values: skim, analyze.")
options.steps = ['skim', 'analyze'] # default value
options.register ('doHLT',
                  0, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Set to 1 to turn on HLT selection for skim step.")
options.register ('doJetFilter',
                  1, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Set to 1 to turn on JetFilter for skim step.")
# Get and parse the command line arguments
options.parseArguments()
print ''
print 'Printing options:'
print options
print 'Only the following options are used: crab, data, sample, steps, doHLT, doJetFilter'
print ''

# Check validity of command line arguments
if options.sample not in sample_options:
    print 'Invalid sample type. Setting to sample type to signal.'
    options.sample = 'signal'
for step in options.steps:
    if step not in steps_options:
        print "Skipping invalid steps: %s" % step
        options.steps.remove(step)

print options.steps

from EmergingJetAnalysis.Configuration.emjetTools import *

process = cms.Process('TEST')
if 'skim' in options.steps and len(options.steps)==1:
    # If only running skim, change process name
    process.setName_('SKIM')

########################################
# Stable configuration
########################################
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

# Unscheduled execution
# process.options.allowUnscheduled = cms.untracked.bool(False)
process.options.allowUnscheduled = cms.untracked.bool(True)

########################################
# Skim
########################################
import os
cmssw_version = os.environ['CMSSW_VERSION']
skimStep = cms.Sequence()
if 'skim' in options.steps:
    print ''
    print '####################'
    print 'Adding Skim step'
    print '####################'
    print ''
    if options.sample=='wjet':
        skimStep = addWJetSkim(process, options.data)
        if 'CMSSW_7_4_12' in cmssw_version:
            process.wJetFilter.electronID = cms.string('cutBasedElectronID-Spring15-25ns-V1-standalone-medium')
        elif 'CMSSW_7_4_1_patch4' in cmssw_version:
            process.wJetFilter.electronID = cms.string('cutBasedElectronID-CSA14-50ns-V1-standalone-medium')
    elif options.sample=='gjet':
       skimStep = addGJetSkim(process, options.data)
    else:
        skimStep = addSkim(process, options.data, doJetFilter=options.doJetFilter, doHLT=options.doHLT)

########################################
# Analyze
########################################
analyzeStep = cms.Sequence()
if 'analyze' in options.steps:
    print ''
    print '####################'
    print 'Adding Analyze step'
    print '####################'
    print ''
    analyzeStep = addAnalyze(process, options.data, options.sample)

########################################
# Testing step
########################################
testing = 0
testingStep = cms.Sequence()
if testing:
    testingStep = addTesting(process, options.data, options.sample)

process.p = cms.Path( skimStep * testingStep * analyzeStep )

########################################
# Configure EDM Output
########################################
if 'skim' in options.steps and len(options.steps)==1:
    # If only running skim, add AOD/AODSIM and jetFilter/wJetFilter to output
    print ''
    print '####################'
    print 'Adding EDM output'
    print '####################'
    print ''
    addEdmOutput(process, options.data, options.sample)
else:
    # Otherwise only save EDM output of jetFilter and wJetFilter
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('output.root'),
        outputCommands = cms.untracked.vstring('drop *'),
    )
    if options.sample=='wjet' : process.out.outputCommands.extend(cms.untracked.vstring('keep *_wJetFilter_*_*',))
    elif options.sample=='gjet': process.out.outputCommands.extend(cms.untracked.vstring('keep *_gJetFilter_*_*',))
    else                      : process.out.outputCommands.extend(cms.untracked.vstring('keep *_jetFilter_*_*',))

########################################
# Generic configuration
########################################
if 'CMSSW_7_4_12' in cmssw_version:
    globalTags=['74X_mcRun2_design_v2','74X_dataRun2_Prompt_v3']
elif 'CMSSW_7_4_1_patch4' in cmssw_version:
    globalTags=['MCRUN2_74_V9','74X_dataRun2_Prompt_v0']
elif 'CMSSW_7_6_3' in cmssw_version:
    globalTags=['76X_mcRun2_asymptotic_RunIIFall15DR76_v1','76X_dataRun2_16Dec2015_v0']
elif 'CMSSW_8_0_26_patch1' in cmssw_version:
    # globalTags=['80X_mcRun2_asymptotic_2016_miniAODv2_v1','80X_dataRun2_2016SeptRepro_v7']
    globalTags=['80X_mcRun2_asymptotic_2016_TrancheIV_v8','80X_dataRun2_2016SeptRepro_v7']
else: print 'No global tag specified for CMSSW_VERSION: %s' % cmssw_version
print 'CMSSW_VERSION is %s' % cmssw_version
print 'Using the following global tags [MC, DATA]:'
print globalTags
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTags[options.data], '')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.MessageLogger.cerr.FwkReport.limit = 20
process.MessageLogger.cerr.default.limit = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(300) )

process.source = cms.Source("PoolSource",
    # eventsToProcess = cms.untracked.VEventRange("1:36:3523-1:36:3523"),
    fileNames = cms.untracked.vstring(
        # File with single dark pions
        # 'file:/afs/cern.ch/user/y/yoshin/work/public/temp/step2_dark.root'
        # Model A
        # '/store/group/phys_exotica/EmergingJets/EmergingJets_ModelA_TuneCUETP8M1_13TeV_pythia8Mod/AODSIM-v1/160201_201550/0000/aodsim_1.root',
        'file:/store/user/yoshin/EmJetMC/AODSIM/EmergingJets_mass_X_d_1000_mass_pi_d_2_tau_pi_d_5_TuneCUETP8M1_13TeV_pythia8Mod/AODSIM-2017-04-30-test/170430_233938/0000/aodsim_8.root',
        # 80X WJet MC
        # '/store/mc/RunIISummer16DR80Premix/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/023963AD-CBBE-E611-AC12-D4AE526A048B.root',
        # 80X Data
        # '/store/data/Run2016B/JetHT/AOD/23Sep2016-v1/90000/0024BC8F-AC82-E611-9019-001E675817A4.root'
        # 80X QCD
        # '/store/mc/RunIISummer16DR80Premix/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/007BFA96-C3B0-E611-90E0-047D7BD6DD64.root'
    ),
)


if options.data==0:
    # Produce PDF weights (maximum is 3)
    process.pdfWeights = cms.EDProducer("PdfWeightProducer",
        # Fix POWHEG if buggy (this PDF set will also appear on output,
        # so only two more PDF sets can be added in PdfSetNames if not "")
        #FixPOWHEG = cms.untracked.string("cteq66.LHgrid"),
        #GenTag = cms.untracked.InputTag("genParticles"),
        PdfInfoTag = cms.untracked.InputTag("generator"),
        PdfSetNames = cms.untracked.vstring(
                "CT14nlo.LHgrid",
                "NNPDF30_nlo_as_0118.LHgrid",
                "NNPDF23_lo_as_0130_qed.LHgrid",
                # , "MRST2006nnlo.LHgrid"
                # , "NNPDF10_100.LHgrid"
        )
    )


process.TFileService = cms.Service("TFileService", fileName = cms.string("ntuple.root") )
# process.TFileService = cms.Service("TFileService", fileName = cms.string("ntuple-74X-DRtest-20170425.root") )
# process.TFileService = cms.Service("TFileService", fileName = cms.string("ntuple-74X-QCD-20170425.root") )

# storage
process.outpath = cms.EndPath(process.out)

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
# # Needed for GetTrackTrajInfo
# process.load("RecoTracker.Configuration.RecoTracker_cff")

# process.load('Configuration.StandardSequences.Reconstruction_cff') #new for navigation
# process.load('Configuration.StandardSequences.GeometryExtended_cff') #new for navigation
# # process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff') #new for navigation
# # process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
# # process.load('JetMETCorrections.Configuration.CorrectedJetProducers_cff')
# # # #get the jet energy corrections from the db file
# # process.load("CondCore.CondDB.CondDB_cfi") 
