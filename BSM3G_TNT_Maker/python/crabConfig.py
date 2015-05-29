from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'NTuples_Data2015_v3'
config.section_('JobType')
config.JobType.psetName = 'miniAOD.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['OutTree.root']
config.section_('Data')
config.Data.inputDataset = '/DoubleMuParked/CMSSW_7_4_0_pre9_ROOT6-GR_R_74_V8_1Apr_RelVal_dm2012D-v2/MINIAOD'
config.Data.unitsPerJob = 50
config.Data.splitting = 'LumiBased'
config.Data.outLFN = '/store/user/florez/Ntuples_Data2015_Test'
config.Data.publication = False
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'
