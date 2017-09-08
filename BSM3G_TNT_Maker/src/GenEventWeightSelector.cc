// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#include "NtupleMaker/BSM3G_TNT_Maker/interface/GenEventWeightSelector.h"

GenEventWeightSelector::GenEventWeightSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector&& iCC):baseTree(name,tree,debug) {
  if(debug) std::cout << "BSM3G TNT Maker: In the GenEventWeightSelector Constructor --> calling SetBranches()." << std::endl;
  generatorToken_               = iCC.consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genweights"));
  _is_data                      = iConfig.getParameter<bool>("is_data");
  if(!_is_data){
    SetBranches();
  }
}

GenEventWeightSelector::~GenEventWeightSelector(){
  delete tree_;
}

void GenEventWeightSelector::Fill(const edm::Event& iEvent){
  Clear(); 

  if(!_is_data) {
    if(debug_) std::cout << "     GenEventWeightSelector: Re-initialized the event weight and grabbing the gen event info collection handle." << std::endl;

    // grabbing the gen event weight info
    edm::Handle<GenEventInfoProduct> genEvt;
    iEvent.getByToken(generatorToken_,genEvt);
    weightevt = genEvt->weight();
  
    if(debug_) std::cout << "     GenEventWeightSelector: Extracted the event weight." << std::endl;
  }
}

void GenEventWeightSelector::SetBranches(){
  if(debug_) std::cout << "     GenEventWeightSelector: Setting branches by calling AddBranch of baseTree." << std::endl;

  AddBranch(&weightevt            ,"weightevt");

  if(debug_) std::cout << "     GenEventWeightSelector: Finished setting branches." << std::endl;
}

void GenEventWeightSelector::Clear(){
  weightevt = 1;
}
