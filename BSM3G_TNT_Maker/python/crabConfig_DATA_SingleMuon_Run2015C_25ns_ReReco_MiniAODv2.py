from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'SingleMuon_Run2015C_ReReco_25ns_MiniAODv2_December2015'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'miniAODv2_7415_25ns_DATA_ReReco2012C.py'
config.JobType.outputFiles = ['OutTree.root']

config.Data.inputDataset = '/SingleMuon/Run2015C_25ns-05Oct2015-v1/MINIAOD'
config.Data.lumiMask = 'Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/gurrola/'
config.Data.publication = True
config.Data.publishDataName = 'SingleMuon_Run2015C_ReReco_25ns_MiniAODv2_December2015'

config.Site.storageSite = 'T2_US_Vanderbilt'
