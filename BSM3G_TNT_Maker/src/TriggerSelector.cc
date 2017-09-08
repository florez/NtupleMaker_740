#include "NtupleMaker/BSM3G_TNT_Maker/interface/TriggerSelector.h"

#include <algorithm>


TriggerSelector::TriggerSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector&& iCC, std::vector<string>& _triggernames):baseTree(name,tree,debug){
  if(debug) std::cout << "BSM3G TNT Maker: In the TriggerSelector Constructor --> getting parameters & calling SetBranches()." << std::endl;
  triggerResultsTag_  = iCC.consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"));
  triggerPrescales_   = iCC.consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescales"));
  triggernames= &_triggernames;
  SetBranches();
}

TriggerSelector::~TriggerSelector(){
  delete tree_;
}

void TriggerSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerResultsTag_, triggerBits);
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  if(debug_) std::cout << "     TriggerSelector: Cleared the vectors, grabbed the trigger collection handles, and looping over triggers." << std::endl;

  if(triggerBits.isValid() == true ){
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    //std::cout << "\n === TRIGGER PATHS === " << std::endl;
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
      //std::cout << "Trigger " << names.triggerName(i) << ", prescale " << triggerPrescales->getPrescaleForIndex(i) << ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") << std::endl;
      if(names.triggerName(i).find("HLT_") != string::npos) {
        vector<string>::iterator nameit= std::find(triggernames->begin(), triggernames->end(), names.triggerName(i));
        if(nameit == triggernames->end()){
          triggernames->push_back(names.triggerName(i));
          //this will be the last element
          if(triggerBits->accept(i)) {
            Trigger_decision.push_back(triggernames->size()-1);
            Trigger_prescale.push_back(triggerPrescales->getPrescaleForIndex(i));
          }
        }else{
          ptrdiff_t pos=nameit-triggernames->begin();
          if(triggerBits->accept(i)) {
            Trigger_decision.push_back(int(pos));
            Trigger_prescale.push_back(triggerPrescales->getPrescaleForIndex(i));
          }
        }
      }
    }
  }


}

void TriggerSelector::SetBranches(){
  if(debug_) std::cout << "     TriggerSelector: Setting branches by calling AddBranch of baseTree." << std::endl;

  AddBranch(&Trigger_decision,   "Trigger_decision");
  AddBranch(&Trigger_prescale,   "Trigger_prescale");

  if(debug_) std::cout << "     TriggerSelector: Finished setting branches." << std::endl;
}

void TriggerSelector::Clear(){
  Trigger_decision.clear();
  Trigger_prescale.clear();
}

void TriggerSelector::startTrigger(edm::EventSetup const& iSetup, edm::Run const & iRun){
}
