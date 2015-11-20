import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')

process.load('Configuration.Geometry.GeometryRecoDB_cff')
#process.load('Configuration.Geometry.GeometryIdeal_cff')

process.load("Configuration.StandardSequences.MagneticField_38T_PostLS1_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")

#needs to be commented out or else it complains about: 'two EventSetup Producers want to deliver type="CaloSubdetectorGeometry" label="ZDC"'
#process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# NOTE: the pick the right global tag!
#    for Spring15 50ns MC: global tag is 'auto:run2_mc_50ns'
#    for Spring15 25ns MC: global tag is 'auto:run2_mc'
#    for Run 2 data: global tag is 'auto:run2_data'
#  as a rule, find the "auto" global tag in $CMSSW_RELEASE_BASE/src/Configuration/AlCa/python/autoCond.py
#  This auto global tag will look up the "proper" global tag
#  that is typically found in the DAS under the Configs for given dataset
#  (although it can be "overridden" by requirements of a given release)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    '/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/00759690-D16E-E511-B29E-00261894382D.root'
  )
)

# START ELECTRON ID SECTION
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string("OutTree.root")
)

process.TNT = cms.EDAnalyzer("BSM3G_TNT_Maker",
    # boolean variables
    debug_ = cms.bool(False),
    filltriggerinfo     = cms.bool(True),
    fillmuoninfo        = cms.bool(True),
    fillelectronpatinfo = cms.bool(True),
    filltauinfo         = cms.bool(True),
    fillgeninfo         = cms.bool(True),
    fillgenweightinfo   = cms.bool(True),
    fillruninfo         = cms.bool(True),
    fillPVinfo          = cms.bool(True),
    filljetinfo         = cms.bool(True),
    fillMETinfo         = cms.bool(True),
    fillphotoninfo      = cms.bool(True),   

    # make a super tiny ntuple, only with a few branches?
    super_TNT  = cms.bool(False),
    # is data or MC?
    is_data     = cms.bool(False),
 
    # input tags 
    triggerResults      = cms.InputTag( 'TriggerResults', '', 'HLT' ),
    triggerPrescales    = cms.InputTag("patTrigger"),
    vertices            = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot            = cms.InputTag("offlineBeamSpot"),
    pileupInfo          = cms.InputTag("slimmedAddPileupInfo"),
    muons               = cms.InputTag("slimmedMuons"),
    patElectrons        = cms.InputTag("slimmedElectrons"),
    taus                = cms.InputTag("slimmedTaus"),
    photons             = cms.InputTag("slimmedPhotons"),
    jets                = cms.InputTag("slimmedJets"),
    jetsPUPPI           = cms.InputTag("slimmedJetsPuppi"),
    mets                = cms.InputTag("slimmedMETs"),
    metsPUPPI           = cms.InputTag("slimmedMETsPuppi"),
    metsNoHF            = cms.InputTag("slimmedMETsNoHF"),
    packedPFCandidates  = cms.InputTag("packedPFCandidates"),
    electronVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    electronLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    eleHEEPIdMap        = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),

    # muon cuts
    Muon_pt_min              = cms.double(8.0),
    Muon_eta_max             = cms.double(2.4),
    Muon_vtx_ndof_min        = cms.int32(4),
    Muon_vtx_rho_max         = cms.int32(2),
    Muon_vtx_position_z_max  = cms.double(24.),

    # electron cuts
    patElectron_pt_min       = cms.double(10.0),
    patElectron_eta_max      = cms.double(2.4),
    patElectron_vtx_ndof_min        = cms.int32(4),
    patElectron_vtx_rho_max         = cms.int32(2),
    patElectron_vtx_position_z_max  = cms.double(24.),

    # tau cuts
    Tau_pt_min  = cms.double(20.),
    Tau_eta_max = cms.double(2.3),
    Tau_vtx_ndof_min        = cms.int32(4),
    Tau_vtx_rho_max         = cms.int32(2),
    Tau_vtx_position_z_max  = cms.double(24.),

    # jet cuts
    Jet_pt_min   = cms.double(30.),

    # photon cuts 
    Photon_pt_min   = cms.double(10.0),
    Photon_eta_max  = cms.double(5.0),    

    # primary vertex cuts
    Pvtx_ndof_min   = cms.double(4.),
    Pvtx_vtx_max  = cms.double(24.),
    Pvtx_vtxdxy_max = cms.double(24.),

)

# filter out anomalous MET from detector noise, cosmic rays, and beam halo particles
process.load("RecoMET.METFilters.metFilters_cff")

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24),
                                           maxd0 = cms.double(2)
)

##___________________________HCAL_Noise_Filter________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
   reverseDecision = cms.bool(False)
)

process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
   reverseDecision = cms.bool(False)
)

# met filters using triggerResults
process.load("NtupleMaker.BSM3G_TNT_Maker.triggerReq_cfi")

# re-calculate MET significance and the covariance matrix
process.load("RecoMET/METProducers.METSignificance_cfi")
process.load("RecoMET/METProducers.METSignificanceParams_cfi")

process.p = cms.Path(process.primaryVertexFilter * 
                     process.CSCTightHaloFilterPAT *
                     process.eeBadScFilterPAT *
                     process.HBHENoiseFilterResultProducer * 
                     process.ApplyBaselineHBHENoiseFilter * 
                     process.egmGsfElectronIDSequence * 
                     process.METSignificance * 
                     process.TNT
)
#process.p = cms.Path(process.TNT)
