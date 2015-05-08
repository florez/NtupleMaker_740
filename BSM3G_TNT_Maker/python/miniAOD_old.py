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
     #'root://xrootd.unl.edu//store/mc/Phys14DR/DYToMuMu_M-50_Tune4C_13TeV-pythia8/MINIAODSIM/PU40bx25_tsg_castor_PHYS14_25_V1-v2/00000/622CAFBA-BD9A-E411-BE11-002481E14FFC.root'
      #'file:reco.root'
      #'file:/uscms_data/d2/florez/BSM3G_TNT_Maker/CMSSW_7_2_0/src/NtupleMaker/BSM3G_TNT_Maker/python/0432E62A-7A6C-E411-87BB-002590DB92A8.root'
      #'/store/mc/Phys14DR/DYToMuMu_M-50_Tune4C_13TeV-pythia8/MINIAODSIM/PU20bx25_castor_PHYS14_25_V1-v3/00000/08CFF3B5-B99D-E411-9FBC-002590200AE8.root',
      #'/store/mc/Phys14DR/QCD_Pt-15to30_Tune4C_13TeV_pythia8/MINIAODSIM/AVE30BX50_tsg_castor_PHYS14_ST_V1-v2/00000/06180282-6B93-E411-8B12-002590200B40.root',
      '/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/20000/002E1374-5F84-E211-83C4-20CF305616D0.root',
      #'/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/08967E6D-6CDC-E411-92DF-0025905B859E.root',
      #'/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/180224F9-F6DB-E411-BDE9-0025905A6056.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/1A1A1EA4-6EDC-E411-8275-0025905A48BC.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/1AF47D60-8CDC-E411-86C2-0025905B85D8.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/22DC35AE-6EDC-E411-BE3F-0025905A48EC.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/2888A658-74DC-E411-9435-0025905A6080.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/28A77600-F7DB-E411-9F2A-002354EF3BE2.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/4AC2D46C-6CDC-E411-8F82-0025905A48F2.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/56E8CF6C-71DC-E411-A14E-0025905A6080.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/7E2AADA0-6EDC-E411-9B5F-0025905A60D6.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/84ACA646-7EDC-E411-A9D5-0025905B858E.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/88C2E163-70DC-E411-BC6F-0025905B859E.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/8A07B311-6EDC-E411-9472-0025905B861C.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/8A087DAF-6EDC-E411-B12D-0025905B858E.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/9C8F946D-8CDC-E411-85AF-0025905B858C.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/AA5CD26A-6CDC-E411-9512-0025905B861C.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/B41E8F5A-71DC-E411-89DA-0025905A6060.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/BE37FF6C-6CDC-E411-9927-0025905B85B2.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/DC4190B1-6EDC-E411-840F-0025905A60D2.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/E4E66694-1CDC-E411-9DD6-0025905B85D8.root',
      # '/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleMuParked/MINIAOD/GR_R_74_V8_1Apr_RelVal_dm2012D-v2/00000/FA4AEEAE-6EDC-E411-82F2-0025905A6060.root'
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
fileName = cms.string("OutTree.root")
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
