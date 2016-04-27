import FWCore.ParameterSet.Config as cms

#set up a process , for e.g. POSTSKIM in this case
processName = "POSTSKIM"
process = cms.Process(processName)

# this inputs the input files
process.source = cms.Source ("PoolSource",
                             fileNames=cms.untracked.vstring(
                                 'file:/afs/cern.ch/user/y/yoshin/www/emjet-skim/myoutput101_numEvent100.root',
                             )
)



# loads your analyzer
process.pfjetToTrack = cms.EDProducer('PFJetTrackAnalyzer',
                                      srcJets = cms.InputTag("ak4PFJetsCHS"),
)

# talk to output module
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("output.root"),
                               outputCommands = cms.untracked.vstring(
                                   'keep *',
                               )
)

# Defines which modules and sequences to run
process.mypath = cms.Path(process.pfjetToTrack)

# A list of analyzers or output modules to be run after all paths have been run.
process.outpath = cms.EndPath(process.out)
