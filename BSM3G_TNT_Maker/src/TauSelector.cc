#include "NtupleMaker/BSM3G_TNT_Maker/interface/TauSelector.h"

TauSelector::TauSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):
 baseTree(name,tree,debug){
  tauToken_       = iConfig.getParameter<edm::InputTag>("taus");
  _vertexInputTag = iConfig.getParameter<edm::InputTag>("vertices");
  _Tau_pt_min     = iConfig.getParameter<double>("Tau_pt_min");
  _Tau_eta_max    = iConfig.getParameter<double>("Tau_eta_max");
  _super_TNT      = iConfig.getParameter<bool>("super_TNT");
  SetBranches();
}

TauSelector::~TauSelector(){
  delete tree_;
}

void TauSelector::Fill(const edm::Event& iEvent){
  Clear();
  
  int taucoun = 0;
  edm::Handle<pat::TauCollection> taus;
  iEvent.getByLabel(tauToken_, taus);
  for (const pat::Tau &tau : *taus) {
    if (tau.pt() < _Tau_pt_min) continue;
    if(fabs(tau.eta()) > _Tau_eta_max) continue;
    
    Tau_eta.push_back( tau.eta());
    Tau_phi.push_back( tau.phi());
    Tau_pt.push_back( tau.pt());
    Tau_energy.push_back( tau.energy());
    Tau_charge.push_back( tau.charge());
    
    // tau id. discriminators
    Tau_decayModeFinding.push_back( tau.tauID("decayModeFinding"));
    Tau_decayModeFindingNewDMs.push_back( tau.tauID("decayModeFindingNewDMs"));
    Tau_chargedIsoPtSum.push_back(tau.tauID("chargedIsoPtSum"));
    Tau_neutralIsoPtSum.push_back(tau.tauID("neutralIsoPtSum"));
    Tau_againstMuonTight3.push_back( tau.tauID("againstMuonTight3")); 
    Tau_againstElectronMVATightMVA5.push_back( tau.tauID("againstElectronTightMVA5"));
    Tau_nProngs.push_back(tau.signalChargedHadrCands().size());

    if(!_super_TNT){
      Tau_puCorrPtSum.push_back(tau.tauID("puCorrPtSum"));
      Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
      Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
      Tau_byTightCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));
      Tau_byLooseIsolationMVA3newDMwLT.push_back(tau.tauID("byLooseIsolationMVA3newDMwLT"));
      Tau_byLooseIsolationMVA3newDMwoLT.push_back(tau.tauID("byLooseIsolationMVA3newDMwoLT"));
      Tau_byLooseIsolationMva3oldDMwLT.push_back(tau.tauID("byLooseIsolationMVA3oldDMwLT"));
      Tau_byLooseIsolationMVA3oldDMwoLT.push_back(tau.tauID("byLooseIsolationMVA3oldDMwoLT"));
      Tau_byMediumIsolationMVA3newDMwLT.push_back(tau.tauID("byMediumIsolationMVA3newDMwLT"));
      Tau_byMediumIsolationMVA3newDMwoLT.push_back(tau.tauID("byMediumIsolationMVA3newDMwoLT"));
      Tau_byMediumIsolationMva3oldDMwLT.push_back(tau.tauID("byMediumIsolationMVA3oldDMwLT"));
      Tau_byMediumIsolationMVA3oldDMwoLT.push_back(tau.tauID("byMediumIsolationMVA3oldDMwoLT"));
      Tau_byTightIsolationMVA3newDMwLT.push_back(tau.tauID("byTightIsolationMVA3newDMwLT"));
      Tau_byTightIsolationMVA3newDMwoLT.push_back(tau.tauID("byTightIsolationMVA3newDMwoLT"));
      Tau_byTightIsolationMva3oldDMwLT.push_back(tau.tauID("byTightIsolationMVA3oldDMwLT"));
      Tau_byTightIsolationMVA3oldDMwoLT.push_back(tau.tauID("byTightIsolationMVA3oldDMwoLT"));
    
      // discriminators against electrons/muons
    
      Tau_againstMuonLoose2.push_back( tau.tauID("againstMuonLoose2"));
      Tau_againstMuonLoose3.push_back( tau.tauID("againstMuonLoose3"));
      Tau_againstMuonTight2.push_back( tau.tauID("againstMuonTight2"));
      Tau_againstMuonTight3.push_back( tau.tauID("againstMuonTight3"));
    
      Tau_againstElectronMVALooseMVA5.push_back( tau.tauID("againstElectronLooseMVA5"));
      Tau_againstElectronMVAMediumMVA5.push_back( tau.tauID("againstElectronMediumMVA5"));
    
      const reco::CandidatePtr hadTauLeadChargedCand = tau.leadChargedHadrCand();                                                                   
      Tau_leadChargedCandPt.push_back(hadTauLeadChargedCand.isNonnull() ? hadTauLeadChargedCand->pt() : -1.);      
      Tau_leadChargedCandCharge.push_back(hadTauLeadChargedCand.isNonnull() ? hadTauLeadChargedCand->charge() : -2);   
      Tau_leadChargedCandEta.push_back(hadTauLeadChargedCand.isNonnull() ? hadTauLeadChargedCand->eta() : -999);
      Tau_leadChargedCandPhi.push_back(hadTauLeadChargedCand.isNonnull() ? hadTauLeadChargedCand->phi() : -999);
    
    }
    taucoun++;
  }
  
}

void TauSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  
  
  AddBranch(&Tau_eta,                         "Tau_eta");
  AddBranch(&Tau_phi,                         "Tau_phi");
  AddBranch(&Tau_pt,                          "Tau_pt");
  AddBranch(&Tau_energy,                      "Tau_energy");
  AddBranch(&Tau_charge,                      "Tau_charge");
  AddBranch(&Tau_decayModeFinding,            "Tau_decayModeFinding");
  AddBranch(&Tau_decayModeFindingNewDMs,      "Tau_decayModeFindingNewDMs");
  AddBranch(&Tau_chargedIsoPtSum,             "Tau_chargedIsoPtSum");
  AddBranch(&Tau_neutralIsoPtSum,             "Tau_neutralIsoPtSum");
  AddBranch(&Tau_againstMuonTight3,           "Tau_againstMuonTight3"); 
  AddBranch(&Tau_againstElectronMVATightMVA5, "Tau_againstElectronMVATightMVA5");
  AddBranch(&Tau_nProngs,                     "Tau_nProngs");

  if(!_super_TNT){ 
    AddBranch(&Tau_puCorrPtSum,             "Tau_puCorrPtSum");
    AddBranch(&Tau_byLooseCombinedIsolationDeltaBetaCorr,         "Tau_byLooseCombinedIsolationDeltaBetaCorr");
    AddBranch(&Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits,    "Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits");
    AddBranch(&Tau_byMediumCombinedIsolationDeltaBetaCorr,        "Tau_byMediumCombinedIsolationDeltaBetaCorr");
    AddBranch(&Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits,   "Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits");
    AddBranch(&Tau_byTightCombinedIsolationDeltaBetaCorr,         "Tau_byTightCombinedIsolationDeltaBetaCorr");
    AddBranch(&Tau_byTightCombinedIsolationDeltaBetaCorr3Hits,    "Tau_byTightCombinedIsolationDeltaBetaCorr3Hits");
    AddBranch(&Tau_byLooseIsolationMVA3newDMwLT,                  "Tau_byLooseIsolationMVA3newDMwLT");
    AddBranch(&Tau_byLooseIsolationMVA3newDMwoLT,                 "Tau_byLooseIsolationMVA3newDMwoLT");
    AddBranch(&Tau_byLooseIsolationMva3oldDMwLT,                  "Tau_byLooseIsolationMva3oldDMwLT");
    AddBranch(&Tau_byLooseIsolationMVA3oldDMwoLT,                 "Tau_byLooseIsolationMVA3oldDMwoLT");
    AddBranch(&Tau_byMediumIsolationMVA3newDMwLT,                 "Tau_byMediumIsolationMVA3newDMwLT");
    AddBranch(&Tau_byMediumIsolationMVA3newDMwoLT,                "Tau_byMediumIsolationMVA3newDMwoLT");
    AddBranch(&Tau_byMediumIsolationMva3oldDMwLT,                 "Tau_byMediumIsolationMva3oldDMwLT");
    AddBranch(&Tau_byMediumIsolationMVA3oldDMwoLT,                "Tau_byMediumIsolationMVA3oldDMwoLT");
    AddBranch(&Tau_byTightIsolationMVA3newDMwLT,                  "Tau_byTightIsolationMVA3newDMwLT");
    AddBranch(&Tau_byTightIsolationMVA3newDMwoLT,                 "Tau_byTightIsolationMVA3newDMwoLT");
    AddBranch(&Tau_byTightIsolationMva3oldDMwLT,                  "Tau_byTightIsolationMva3oldDMwLT");
    AddBranch(&Tau_byTightIsolationMVA3oldDMwoLT,                 "Tau_byTightIsolationMVA3oldDMwoLT");
    AddBranch(&Tau_againstMuonLoose2,                             "Tau_againstMuonLoose2");
    AddBranch(&Tau_againstMuonLoose3,                             "Tau_againstMuonLoose3");
    AddBranch(&Tau_againstMuonTight2,                             "Tau_againstMuonTight2");
    AddBranch(&Tau_againstElectronMVALooseMVA5,                   "Tau_againstElectronMVALooseMVA5");
    AddBranch(&Tau_againstElectronMVAMediumMVA5,                  "Tau_againstElectronMVAMediumMVA5");
    AddBranch(&Tau_byVLooseCombinedIsolationDeltaBetaCorr,        "Tau_byVLooseCombinedIsolationDeltaBetaCorr");
    AddBranch(&Tau_leadChargedCandPt,                              "Tau_leadChargedCandPt");
    AddBranch(&Tau_leadChargedCandCharge,                          "Tau_leadChargedCandCharge");
    AddBranch(&Tau_leadChargedCandEta,                             "Tau_leadChargedCandEta");
    AddBranch(&Tau_leadChargedCandPhi,                             "Tau_leadChargedCandPhi");
  }
}

void TauSelector::Clear(){
  
  Tau_eta.clear(); 
  Tau_phi.clear(); 
  Tau_pt.clear(); 
  Tau_energy.clear(); 
  Tau_charge.clear(); 
  Tau_decayModeFinding.clear(); 
  Tau_decayModeFindingNewDMs.clear(); 
  Tau_chargedIsoPtSum.clear(); 
  Tau_neutralIsoPtSum.clear(); 
  Tau_puCorrPtSum.clear(); 
  Tau_byLooseCombinedIsolationDeltaBetaCorr.clear(); 
  Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.clear(); 
  Tau_byMediumCombinedIsolationDeltaBetaCorr.clear(); 
  Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.clear(); 
  Tau_byTightCombinedIsolationDeltaBetaCorr.clear(); 
  Tau_byTightCombinedIsolationDeltaBetaCorr3Hits.clear(); 
  Tau_byLooseIsolationMVA3newDMwLT.clear();
  Tau_byLooseIsolationMVA3newDMwoLT.clear();
  Tau_byLooseIsolationMva3oldDMwLT.clear();
  Tau_byLooseIsolationMVA3oldDMwoLT.clear();
  Tau_byMediumIsolationMVA3newDMwLT.clear();
  Tau_byMediumIsolationMVA3newDMwoLT.clear();
  Tau_byMediumIsolationMva3oldDMwLT.clear();
  Tau_byMediumIsolationMVA3oldDMwoLT.clear();
  Tau_byTightIsolationMVA3newDMwLT.clear();
  Tau_byTightIsolationMVA3newDMwoLT.clear();
  Tau_byTightIsolationMva3oldDMwLT.clear();
  Tau_byTightIsolationMVA3oldDMwoLT.clear();
  Tau_againstMuonLoose2.clear();
  Tau_againstMuonLoose3.clear();
  Tau_againstMuonTight2.clear();
  Tau_againstMuonTight3.clear();
  Tau_againstElectronMVALooseMVA5.clear();
  Tau_againstElectronMVAMediumMVA5.clear();
  Tau_againstElectronMVATightMVA5.clear();
  Tau_byVLooseCombinedIsolationDeltaBetaCorr.clear();
  Tau_leadChargedCandPt.clear();
  Tau_leadChargedCandCharge.clear();
  Tau_leadChargedCandEta.clear();
  Tau_leadChargedCandPhi.clear();
  Tau_nProngs.clear();
}
