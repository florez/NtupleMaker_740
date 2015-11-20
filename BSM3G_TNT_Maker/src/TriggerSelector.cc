#include "NtupleMaker/BSM3G_TNT_Maker/interface/TriggerSelector.h"

TriggerSelector::TriggerSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout << "BSM3G TNT Maker: In the TriggerSelector Constructor --> getting parameters & calling SetBranches()." << std::endl;
  triggerResultsTag_  = iConfig.getParameter<edm::InputTag>("triggerResults");
  triggerPrescales_   = iConfig.getParameter<edm::InputTag>("triggerPrescales");
  SetBranches();
}

TriggerSelector::~TriggerSelector(){
  delete tree_;
}

void TriggerSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();
  
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByLabel(triggerResultsTag_, triggerBits);
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByLabel(triggerPrescales_, triggerPrescales);

  if(debug_) std::cout << "     TriggerSelector: Cleared the vectors, grabbed the trigger collection handles, and looping over triggers." << std::endl;

  if(triggerBits.isValid() == true ){
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    //std::cout << "\n === TRIGGER PATHS === " << std::endl;
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
      //std::cout << "Trigger " << names.triggerName(i) << ", prescale " << triggerPrescales->getPrescaleForIndex(i) << ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") << std::endl;
      if(names.triggerName(i).find("HLT_") != string::npos) {
      Trigger_names.push_back(names.triggerName(i));
        if(triggerBits->accept(i)) { Trigger_decision.push_back(1); }
        else { Trigger_decision.push_back(0); }
        Trigger_prescale.push_back(triggerPrescales->getPrescaleForIndex(i));
      }
    }
  }

}

void TriggerSelector::SetBranches(){
  if(debug_) std::cout << "     TriggerSelector: Setting branches by calling AddBranch of baseTree." << std::endl;
  
  AddBranch(&Trigger_decision,   "Trigger_decision");
  AddBranch(&Trigger_names,   "Trigger_names");
  AddBranch(&Trigger_prescale,   "Trigger_prescale");

  if(debug_) std::cout << "     TriggerSelector: Finished setting branches." << std::endl;
}

void TriggerSelector::Clear(){
  Trigger_decision.clear();
  Trigger_names.clear();  
  Trigger_prescale.clear();
}

void TriggerSelector::startTrigger(edm::EventSetup const& iSetup, edm::Run const & iRun){
  bool changed(true);
  hltConfig_.init(iRun,iSetup,"HLT",changed);
}
