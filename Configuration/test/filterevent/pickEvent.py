import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')
# add a list of strings for events to process
options.register ('eventsToProcess',
                                  '',
                                  VarParsing.multiplicity.list,
                                  VarParsing.varType.string,
                                  "Events to process")
options.parseArguments()

process = cms.Process("PickEvent")
process.source = cms.Source ("PoolSource",
          fileNames = cms.untracked.vstring (
              #file name, can be either local or remote. If local, add file: in the front.
              'file:PVFailure/14F69527-D485-E611-9B5D-02163E014789.root',
          ),
          eventsToProcess = cms.untracked.VEventRange (
              # options.eventsToProcess
              #event coordinate: run:lumi:event
              '281641:60:99516202',
          )                               
)

process.Out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string ('PVFailure_data.root')#output file name
)

process.end = cms.EndPath(process.Out)
