import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

HBHENoiseFilterPAT = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::PAT'),
    triggerConditions = (
      'Flag_HBHENoiseFilter',
    ),
    l1tResults = '',
    throw = False
)

HBHENoiseIsoFilterPAT = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::PAT'),
    triggerConditions = (
      'Flag_HBHENoiseIsoFilter',
    ),
    l1tResults = '',
    throw = False
)

EcalDeadCellTriggerPrimitiveFilterPAT = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::PAT'),
    triggerConditions = (
      'Flag_EcalDeadCellTriggerPrimitiveFilter',
    ),
    l1tResults = '',
    throw = False
)

eeBadScFilterPAT = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::PAT'),
    triggerConditions = (
      'Flag_eeBadScFilter',
    ),
    l1tResults = '',
    throw = False
)

CSCTightHaloFilterPAT = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::PAT'),
    triggerConditions = (
      'Flag_globalTightHalo2016Filter',
    ),
    l1tResults = '',
    throw = False
)

goodVerticesFilterPAT = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::PAT'),
    triggerConditions = (
      'Flag_goodVertices',
    ),
    l1tResults = '',
    throw = False
)





HBHENoiseFilterRECO = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::RECO'),
    triggerConditions = (
      'Flag_HBHENoiseFilter',
    ),
    l1tResults = '',
    throw = False
)

HBHENoiseIsoFilterRECO = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::RECO'),
    triggerConditions = (
      'Flag_HBHENoiseIsoFilter',
    ),
    l1tResults = '',
    throw = False
)

EcalDeadCellTriggerPrimitiveFilterRECO = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::RECO'),
    triggerConditions = (
      'Flag_EcalDeadCellTriggerPrimitiveFilter',
    ),
    l1tResults = '',
    throw = False
)

eeBadScFilterRECO = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::RECO'),
    triggerConditions = (
      'Flag_eeBadScFilter',
    ),
    l1tResults = '',
    throw = False
)

CSCTightHaloFilterRECO = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::RECO'),
    triggerConditions = (
      'Flag_globalTightHalo2016Filter',
    ),
    l1tResults = '',
    throw = False
)

goodVerticesFilterRECO = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults::RECO'),
    triggerConditions = (
      'Flag_goodVertices',
    ),
    l1tResults = '',
    throw = False
)




HBHENoiseFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults'),
    triggerConditions = (
      'Flag_HBHENoiseFilter',
    ),
    l1tResults = '',
    throw = False
)

HBHENoiseIsoFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults'),
    triggerConditions = (
      'Flag_HBHENoiseIsoFilter',
    ),
    l1tResults = '',
    throw = False
)

EcalDeadCellTriggerPrimitiveFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults'),
    triggerConditions = (
      'Flag_EcalDeadCellTriggerPrimitiveFilter',
    ),
    l1tResults = '',
    throw = False
)

eeBadScFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults'),
    triggerConditions = (
      'Flag_eeBadScFilter',
    ),
    l1tResults = '',
    throw = False
)

CSCTightHaloFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults'),
    triggerConditions = (
      'Flag_globalTightHalo2016Filter',
    ),
    l1tResults = '',
    throw = False
)

goodVerticesFilter = hlt.triggerResultsFilter.clone(
    hltResults = cms.InputTag('TriggerResults'),
    triggerConditions = (
      'Flag_goodVertices',
    ),
    l1tResults = '',
    throw = False
)


