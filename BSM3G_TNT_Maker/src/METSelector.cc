// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#include "NtupleMaker/BSM3G_TNT_Maker/interface/METSelector.h"

METSelector::METSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in METSelector constructor"<<std::endl;
  metToken_ = iConfig.getParameter<edm::InputTag>("mets");
  if(debug) std::cout<<"in pileup constructor: calling SetBrances()"<<std::endl;
  _is_data   = iConfig.getParameter<bool>("is_data");
  _super_TNT = iConfig.getParameter<bool>("super_TNT");
  SetBranches();
}

METSelector::~METSelector(){
  delete tree_;
}

void METSelector::Fill(const edm::Event& iEvent){
  
  if(debug_)    std::cout<<"getting met info"<<std::endl;
  
  Clear();
  
  //edm::Handle<edm::View<pat::MET> > mets;
  edm::Handle<pat::METCollection> mets;
  iEvent.getByLabel(metToken_, mets);
  const pat::MET &met = mets->front();
 
  if(debug_)    std::cout<<"Filling met branches"<<std::endl;
 
  Met_pt = met.pt();
  Met_phi = met.phi();
  Met_sumEt = met.sumEt();  
  if(!_super_TNT){
    Met_px = met.px();
    Met_py = met.py();
    Met_pz = met.pz();
    if (!_is_data) Gen_Met = met.genMET()->pt();
    Met_shiftedPtUp = met.shiftedPt(pat::MET::JetEnUp);
    Met_shiftedPtDown  = met.shiftedPt(pat::MET::JetEnDown);
  }
  if(debug_)    std::cout<<"got MET info"<<std::endl;
}

void METSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  
  AddBranch(&Met_pt,            "Met_pt");
  AddBranch(&Met_sumEt,         "Met_sumEt");
  AddBranch(&Met_phi,           "Met_phi");

  if(!_super_TNT){
    AddBranch(&Met_px,            "Met_px");
    AddBranch(&Met_py,            "Met_py");
    AddBranch(&Met_pz,            "Met_pz");
    AddBranch(&Gen_Met,           "Gen_Met");
    AddBranch(&Met_shiftedPtUp,   "Met_shiftedPtUp");
    AddBranch(&Met_shiftedPtDown, "Met_shiftedPtDown");
  }
  if(debug_)    std::cout<<"set branches"<<std::endl;
}


void METSelector::Clear(){
  Met_pt = -9999999999;                                                                                                                                    
  Met_px = -9999999999;
  Met_py = -9999999999;
  Met_pz = -9999999999;
  Met_phi  = -9999999999;
  Met_sumEt  = -9999999999 ;
  Gen_Met  = -9999999999;
  Met_shiftedPtUp  = -9999999999;
  Met_shiftedPtDown   = -9999999999;
  
}
