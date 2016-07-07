from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'RunIISpring15MiniAODv2_WJetsToLNu_madgraphMLM_Asympt25ns'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'miniAODv2_7415_25ns_MC.py'
config.JobType.outputFiles = ['OutTree.root']

config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/gurrola/'
config.Data.publication = True
config.Data.publishDataName = 'WJetsToLNu_madgraphMLM_Asympt25ns_TNTMaker_MiniAODv2'

config.Site.storageSite = 'T2_US_Vanderbilt'
