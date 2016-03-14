#!/bin/bash

STARTINGDIR=${PWD}
echo $STARTINGDIR
# cd cogFiles
export PYTHONPATH=${PWD}/cogFiles
echo "WARNING: Replacing EmergingJetAnalyzer/interface/OutputTree.h"
cog.py -r ${STARTINGDIR}/EmergingJetAnalyzer/interface/OutputTree.h
echo "WARNING: Replacing EmergingJetAnalyzer/plugins/EmergingJetAnalyzer.cc"
cog.py -r ${STARTINGDIR}/EmergingJetAnalyzer/plugins/EmergingJetAnalyzer.cc
USER_CXXFLAGS="-Wno-error=unused-variable -Wno-error=unused-but-set-variable -DEDM_ML_DEBUG -g" scram b -v -j8
cd $STARTINGDIR
