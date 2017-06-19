// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#include "NtupleMaker/BSM3G_TNT_Maker/interface/RunInfoSelector.h"

RunInfoSelector::RunInfoSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector&& iCC):baseTree(name,tree,debug) {
  if(debug) std::cout << "BSM3G TNT Maker: In the RunInfoSelector Constructor --> calling SetBranches()." << std::endl;

  SetBranches();
}

RunInfoSelector::~RunInfoSelector(){
  delete tree_;
}

void RunInfoSelector::Fill(const edm::Event& iEvent){
  Clear();

  if(debug_) std::cout << "     RunInfoSelector: Re-initialized the variables and extracting the run number, event number, and lumi block." << std::endl;

  runNumber = iEvent.id().run();
  eventNumber = iEvent.id().event();
  if( iEvent.isRealData() ) {
    lumiBlock = iEvent.luminosityBlock();
  } else {
    lumiBlock = -1;
  }


}

void RunInfoSelector::SetBranches(){
  if(debug_) std::cout << "     RunInfoSelector: Setting branches by calling AddBranch of baseTree." << std::endl;

  AddBranch(&runNumber            ,"runNumber");
  AddBranch(&eventNumber          ,"eventNumber");
  AddBranch(&lumiBlock            ,"lumiBlock");


  if(debug_) std::cout << "     RunInfoSelector: Finished setting branches." << std::endl;
}

void RunInfoSelector::Clear(){
  runNumber = -1;
  eventNumber = -1;
  lumiBlock = -1;

}
