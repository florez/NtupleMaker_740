#include "NtupleMaker/BSM3G_TNT_Maker/interface/TriggerSelector.h"

TriggerSelector::TriggerSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  triggerResultsTag_  = iConfig.getParameter<edm::InputTag>("triggerResults");
  SetBranches();
}

TriggerSelector::~TriggerSelector(){
  delete tree_;
}

void TriggerSelector::Fill(const edm::Event& iEvent){
  
  if(debug_)    std::cout<<"getting met info"<<std::endl;
   
  Clear();
  
  edm::Handle<edm::TriggerResults>      triggerResultsHandle_;
  iEvent.getByLabel(triggerResultsTag_, triggerResultsHandle_);

  if(triggerResultsHandle_.isValid() == true ){
     const edm::TriggerNames vTrigNames = iEvent.triggerNames(*triggerResultsHandle_);
     const unsigned int ntrigs = triggerResultsHandle_->size();

     for(unsigned int u = 0; u < hltConfig_.triggerNames().size(); u++){
       const unsigned int triggerIndex(hltConfig_.triggerIndex(hltConfig_.triggerNames().at(u)));
       if( (ntrigs > triggerIndex) ){
           Trigger_names.push_back(hltConfig_.triggerNames().at(u));
           if(triggerResultsHandle_->accept(triggerIndex) ){
             Trigger_decision.push_back(1);
           }else {
             Trigger_decision.push_back(0);
           }
       }
     }
  }
  if(debug_)    std::cout<<"got Trigger info"<<std::endl;
}

void TriggerSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  
  AddBranch(&Trigger_decision,   "Trigger_decision");
  AddBranch(&Trigger_names,   "Trigger_names");
  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void TriggerSelector::Clear(){
  Trigger_decision.clear();
  Trigger_names.clear();  
}

void TriggerSelector::startTrigger(edm::EventSetup const& iSetup, edm::Run const & iRun){
  bool changed(true);
  hltConfig_.init(iRun,iSetup,"HLT",changed);
}
