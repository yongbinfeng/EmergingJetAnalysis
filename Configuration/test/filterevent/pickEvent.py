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
              'file:5E77AE43-D6B0-E611-BC6B-FA163ED28262.root'
          ),
          eventsToProcess = cms.untracked.VEventRange (
              # options.eventsToProcess
              #event coordinate: run:lumi:event
              '1:39938:115533454',
          )                               
)

process.Out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string ('alphamax_newworse.root')#output file name
)

process.end = cms.EndPath(process.Out)
