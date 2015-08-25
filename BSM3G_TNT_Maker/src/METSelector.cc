// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#include "NtupleMaker/BSM3G_TNT_Maker/interface/METSelector.h"

METSelector::METSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout<<"in METSelector constructor"<<std::endl;
  metToken_ = iConfig.getParameter<edm::InputTag>("mets");
  puppi_metToken_ = iConfig.getParameter<edm::InputTag>("metsPUPPI");
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
  
  // grab handle to the met collections (as of 08/21/2015, the current POG recommendation is type1 PF Met)
  edm::Handle<pat::METCollection> mets;
  iEvent.getByLabel(metToken_, mets);
  const pat::MET &met = mets->front();
  edm::Handle<pat::METCollection> puppimets;
  iEvent.getByLabel(puppi_metToken_, puppimets);
  const pat::MET &puppimet = puppimets->front();
 
  if(debug_)    std::cout<<"Filling met branches"<<std::endl;
 
  // type-1 PF Met:  propagation of the jet energy corrections to MET. 
  // The jet collection is AK4 CHS jets (non-CHS ak4PFJets) with Pt>10 GeV/c using the L1L2L3-L1 scheme.
  Met_type1PF_pt = met.pt();
  Met_type1PF_px = met.px();
  Met_type1PF_py = met.py();
  Met_type1PF_pz = met.pz();
  Met_type1PF_phi = met.phi();
  Met_type1PF_sumEt = met.sumEt();  

  // PUPPI Met:  Met reconstructed with the PUPPI pileup mitigation algorithm.
  Met_puppi_pt = puppimet.pt();
  Met_puppi_px = puppimet.px();
  Met_puppi_py = puppimet.py();
  Met_puppi_pz = puppimet.pz();
  Met_puppi_phi = puppimet.phi();
  Met_puppi_sumEt = puppimet.sumEt();

  // store additional information such as generator level
  if(!_super_TNT){
    if (!_is_data) Gen_Met = met.genMET()->pt();
    Met_type1PF_shiftedPtUp = met.shiftedPt(pat::MET::JetEnUp);
    Met_type1PF_shiftedPtDown  = met.shiftedPt(pat::MET::JetEnDown);
//    Met_puppi_shiftedPtUp = puppimet.shiftedPt(pat::MET::JetEnUp);
//    Met_puppi_shiftedPtDown  = puppimet.shiftedPt(pat::MET::JetEnDown);
  }
  if(debug_)    std::cout<<"got MET info"<<std::endl;
}

void METSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  
  AddBranch(&Met_type1PF_pt,            "Met_type1PF_pt");
  AddBranch(&Met_type1PF_px,            "Met_type1PF_px");
  AddBranch(&Met_type1PF_py,            "Met_type1PF_py");
  AddBranch(&Met_type1PF_pz,            "Met_type1PF_pz");
  AddBranch(&Met_type1PF_sumEt,         "Met_type1PF_sumEt");
  AddBranch(&Met_type1PF_phi,           "Met_type1PF_phi");
  AddBranch(&Met_puppi_pt,              "Met_puppi_pt");
  AddBranch(&Met_puppi_px,              "Met_puppi_px");
  AddBranch(&Met_puppi_py,              "Met_puppi_py");
  AddBranch(&Met_puppi_pz,              "Met_puppi_pz");
  AddBranch(&Met_puppi_sumEt,           "Met_puppi_sumEt");
  AddBranch(&Met_puppi_phi,             "Met_puppi_phi");
  if(!_super_TNT){
    AddBranch(&Gen_Met,           "Gen_Met");
    AddBranch(&Met_type1PF_shiftedPtUp,   "Met_type1PF_shiftedPtUp");
    AddBranch(&Met_type1PF_shiftedPtDown, "Met_type1PF_shiftedPtDown");
//    AddBranch(&Met_puppi_shiftedPtUp,     "Met_puppi_shiftedPtUp");
//    AddBranch(&Met_puppi_shiftedPtDown,   "Met_puppi_shiftedPtDown");
  }
  if(debug_)    std::cout<<"set branches"<<std::endl;
}


void METSelector::Clear(){
  Met_type1PF_pt = -9999999999;                                                                                                                                    
  Met_type1PF_px = -9999999999;
  Met_type1PF_py = -9999999999;
  Met_type1PF_pz = -9999999999;
  Met_type1PF_phi  = -9999999999;
  Met_type1PF_sumEt  = -9999999999 ;
  Met_type1PF_shiftedPtUp  = -9999999999;
  Met_type1PF_shiftedPtDown   = -9999999999;
  Met_puppi_pt = -9999999999;                                                                                                                                    
  Met_puppi_px = -9999999999;
  Met_puppi_py = -9999999999;
  Met_puppi_pz = -9999999999;
  Met_puppi_phi  = -9999999999;
  Met_puppi_sumEt  = -9999999999 ;
//  Met_puppi_shiftedPtUp  = -9999999999;
//  Met_puppi_shiftedPtDown   = -9999999999;
  Gen_Met  = -9999999999;
}
