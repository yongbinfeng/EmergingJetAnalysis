"""MultiCRAB script template for EmergingJet analysis"""
MC=0
DATA=1

# this is for another set of Gamma +Jet MC samples. For trackref check.
#   fake rate check

jobname = 'Analysis-20170426'        # Jobname
psetname = 'Configuration/test/test_cfg.py'      # Path to pset
postfix = 'v1'                # Postfix to job, increment for each major version
dryrun = 0

from datasets import dataset
from wrappers import submit, submit_newthread
datasets = [
    #dataset( "WJetsToLNuInclusive" , "/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM" , MC, 1, 2000, splitting='FileBased', priority=30, label='wjet' ),
    # dataset( "WJetSkimMuon"          , "/SingleMuon/yoshin-WJetSkim-ede3f21fae18a825b193df32c86b780e/USER"                                                             , DATA , 10000 , 10000000 , splitting='EventAwareLumiBased' , priority=30 , label='wjet'     , inputDBS='phys03' ) ,
    dataset( "GJets_HT-40To100_RunIISummer16DR80Premix" , "/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM" , MC, 1, 1000, splitting='FileBased', priority=30, label='gjet' ),
    dataset( "GJets_HT-100To200_RunIISummer16DR80Premix" , "/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM" , MC, 1, 1000, splitting='FileBased', priority=30, label='gjet' ),
    dataset( "GJets_HT-200To400_RunIISummer16DR80Premix" , "/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM" , MC, 1, 1000, splitting='FileBased', priority=30, label='gjet' ),
    dataset( "GJets_HT-400To600_RunIISummer16DR80Premix" , "/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM" , MC, 1, 1000, splitting='FileBased', priority=30, label='gjet' ),
    dataset( "GJets_HT-600ToInf_RunIISummer16DR80Premix" , "/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM" , MC, 1, 1000, splitting='FileBased', priority=30, label='gjet' ),
]

if __name__ == '__main__':

    ############################################################
    ## Common settings
    ############################################################
    from WMCore.Configuration import Configuration
    import time
    config = Configuration()

    config.section_("General")
    config.General.workArea = 'crabTasks/' + 'crab_' + jobname + time.strftime("-%Y-%m%d") + '-' + postfix
    tasklistFileName = config.General.workArea + '.txt'
    if not dryrun: tasklistFile = open(tasklistFileName, 'a+')

    config.section_("JobType")
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = psetname
    # config.JobType.scriptExe = 'crab_script.sh'                                # Script to run instead of cmsRun
    # config.JobType.inputFiles = ['emjet-basicJetAnalyzer.py','crab_script.py'] # Additional input files
    # config.JobType.outputFiles = ['histo.root']                                # Collect non-EDM outputs if any

    config.section_("Data")
    config.Data.splitting = 'FileBased'
    config.Data.publication = True
    config.Data.outputDatasetTag = jobname
    config.Data.ignoreLocality = True

    config.section_("Site")
    config.Site.storageSite = "T3_US_UMD"
    config.Site.blacklist = ['T2_US_Purdue']
    # config.Site.whitelist = ['T2_CH_CERN']
    # config.Site.ignoreGlobalBlacklist = True

    ############################################################
    ## Dataset specific settings
    ############################################################
    for dataset in datasets:
        alias = dataset.alias
        pyCfgParams = ['crab=1'] # Temporary list, must be written to config.JobType.pyCfgParams before submitting
        config.General.requestName   = ('%s-%s' + time.strftime("-%Y-%m%d-%H%M%S")) % (jobname, alias)
        config.Data.outLFNDirBase = '/store/user/yofeng/80X/ntuple/%s-%s/%s/' % (jobname, postfix, alias)
        config.Data.inputDataset = dataset.fullpath
        config.Data.unitsPerJob  = dataset.unitsPerJob
        config.Data.totalUnits   = dataset.totalUnits
        isData                   = dataset.isData
        config.Data.splitting    = dataset.splitting
        config.JobType.priority  = dataset.priority
        config.Data.inputDBS     = dataset.inputDBS
        label                    = dataset.label
        doHLT                    = dataset.doHLT
        doJetFilter              = dataset.doJetFilter
        # MC specific settings:
        if not isData:
            pyCfgParams.append('data=0')
        # Data specific settings:
        if isData:
            pyCfgParams.append('data=1')
            config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
        pyCfgParams.append('steps=skim,analyze')
        # Label specific settings
        pyCfgParams.append('sample='+label)
        pyCfgParams.append('doHLT=%d' % doHLT)
        pyCfgParams.append('doJetFilter=%d' % doJetFilter)
        config.JobType.pyCfgParams = pyCfgParams

        if not dryrun:
            res = submit_newthread(config, dryrun=dryrun)
            taskId = res['uniquerequestname'].split(':')[0]
            filePath = '%s%s/%s/%s/' % ( config.Data.outLFNDirBase, config.Data.inputDataset.split('/')[1], config.Data.outputDatasetTag, taskId )
            print 'filePath:'
            print filePath
            tasklistFile.write(filePath)
            tasklistFile.write('\n')

    if not dryrun: tasklistFile.close()

