import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# NOTE: the pick the right global tag!
# for PHYS14 scenario PU4bx50 : global tag is ???
# for PHYS14 scenario PU20bx25: global tag is PHYS14_25_V1
# as a rule, find the global tag in the DAS under the Configs for given dataset
#process.GlobalTag.globaltag = 'PHYS14_25_V1'
process.GlobalTag.globaltag = 'GR_R_74_V8'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'DATAFILE',
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
from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
# Define which IDs we want to produce
# Each of these two example IDs contains all four standard
# cut-based ID working points (only two WP of the PU20bx25 are actually used here).
#Add them to the VID producer
# Do not forget to add the egmGsfElectronIDSequence to the path,
# as in the example below!
#
# END ELECTRON ID SECTION
#

process.TFileService = cms.Service("TFileService",
fileName = cms.string("/eos/uscms/store/user/florez/NTuples_Data2015/OUTPUTFILE.root")
)

process.TNT = cms.EDAnalyzer("BSM3G_TNT_Maker",
    # boolean variables
    debug_ = cms.bool(False),
    fillmuoninfo        = cms.bool(True),
    filltauinfo         = cms.bool(True),
    fillgeninfo         = cms.bool(False),
    fillPVinfo          = cms.bool(True),
    filljetinfo         = cms.bool(True),
    fillMETinfo         = cms.bool(True),

    # make a super tiny ntuple, only with a few branches?
    super_TNT  = cms.bool(False),

    # is data or MC?
    is_data     = cms.bool(True), 
    
    # input tags 
    vertices    = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons       = cms.InputTag("slimmedMuons"),
    taus        = cms.InputTag("slimmedTaus"),
    photons     = cms.InputTag("slimmedPhotons"),
    jets        = cms.InputTag("slimmedJets"),
    fatjets     = cms.InputTag("slimmedJetsAK8"),
    mets        = cms.InputTag("slimmedMETs"),
    bits        = cms.InputTag("TriggerResults","","HLT"),
    prescales   = cms.InputTag("patTrigger"),
    objects     = cms.InputTag("selectedPatTrigger"),  

    # muon cuts
    Muon_pt_min              = cms.double(10.0),
    Muon_eta_max             = cms.double(2.4),
    Muon_vtx_ndof_min        = cms.int32(4),
    Muon_vtx_rho_max         = cms.int32(2),
    Muon_vtx_position_z_max  = cms.double(24.),

    Tau_pt_min    = cms.double(20.),
    Tau_eta_max   = cms.double(2.3),

    # jet cuts
    Jet_pt_min    = cms.double(25.),

    # primary vertex cuts
    Pvtx_ndof_min   = cms.double(4.),
    Pvtx_vtx_max    = cms.double(24.),
    Pvtx_vtxdxy_max = cms.double(24.),
)

process.p = cms.Path(process.TNT)
