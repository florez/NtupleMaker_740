import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# NOTE: the pick the right global tag!
# for PHYS14 scenario PU4bx50 : global tag is ???
# for PHYS14 scenario PU20bx25: global tag is PHYS14_25_V1
# as a rule, find the global tag in the DAS under the Configs for given dataset
process.GlobalTag.globaltag = 'PHYS14_25_V1::All'
#process.GlobalTag.globaltag = 'GR_R_74_V8'
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #'/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/20000/002E1374-5F84-E211-83C4-20CF305616D0.root',
      '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0432E62A-7A6C-E411-87BB-002590DB92A8.root',
    )
)

#
# START ELECTRON ID SECTION
#
# Set up everything that is needed to compute electron IDs and
# add the ValueMaps with ID decisions into the event data stream
#
# Load tools and function definitions
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
# overwrite a default parameter: for miniAOD, the collection name is a slimmed one
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')
from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)
# Define which IDs we want to produce
# Each of these two example IDs contains all four standard
# cut-based ID working points (only two WP of the PU20bx25 are actually used here).
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff']
#Add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
# Do not forget to add the egmGsfElectronIDSequence to the path,
# as in the example below!
#
# END ELECTRON ID SECTION
#

process.TFileService = cms.Service("TFileService",
fileName = cms.string("OutTree.root")
)

process.TNT = cms.EDAnalyzer("BSM3G_TNT_Maker",
    # boolean variables
    debug_ = cms.bool(False),
    fillmuoninfo        = cms.bool(False),
    fillelectronpatinfo = cms.bool(True),
    filltauinfo         = cms.bool(True),
    fillgeninfo         = cms.bool(False),
    fillPVinfo          = cms.bool(False),
    filljetinfo         = cms.bool(False),
    fillMETinfo         = cms.bool(False),
    fillphotoninfo      = cms.bool(True),   

    # make a super tiny ntuple, only with a few branches?
    super_TNT  = cms.bool(False),
    # is data or MC?
    is_data     = cms.bool(False),
 
    # input tags 
    vertices            = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons               = cms.InputTag("slimmedMuons"),
    patElectrons        = cms.InputTag("slimmedElectrons"),
    taus                = cms.InputTag("slimmedTaus"),
    photons             = cms.InputTag("slimmedPhotons"),
    jets                = cms.InputTag("slimmedJets"),
    fatjets             = cms.InputTag("slimmedJetsAK8"),
    mets                = cms.InputTag("slimmedMETs"),
    bits                = cms.InputTag("TriggerResults","","HLT"),
    prescales           = cms.InputTag("patTrigger"),
    objects             = cms.InputTag("selectedPatTrigger"),  
    electronVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
    electronLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
    electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
    electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
    eleHEEPIdMap        = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51"),

    # muon cuts
    Muon_pt_min              = cms.double(10.0),
    Muon_eta_max             = cms.double(2.4),
    Muon_vtx_ndof_min        = cms.int32(4),
    Muon_vtx_rho_max         = cms.int32(2),
    Muon_vtx_position_z_max  = cms.double(24.),

    # electron cuts
    patElectron_pt_min       = cms.double(10.0),
    patElectron_eta_max      = cms.double(2.4),

    # tau cuts
    Tau_pt_min  = cms.double(20.),
    Tau_eta_max = cms.double(2.3),

    # jet cuts
    Jet_pt_min   = cms.double(25.),

    # photon cuts 
    Photon_pt_min   = cms.double(5.0),
    Photon_eta_max  = cms.double(5.0),    

    # primary vertex cuts
    Pvtx_ndof_min   = cms.double(4.),
    Pvtx_vtx_max  = cms.double(24.),
    Pvtx_vtxdxy_max = cms.double(24.),
)

process.p = cms.Path(process.egmGsfElectronIDSequence * process.TNT)
#process.p = cms.Path(process.TNT)
