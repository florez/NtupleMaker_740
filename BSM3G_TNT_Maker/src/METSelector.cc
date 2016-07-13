// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#include "NtupleMaker/BSM3G_TNT_Maker/interface/METSelector.h"

METSelector::METSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector&& iCC):baseTree(name,tree,debug){
  if(debug) std::cout << "BSM3G TNT Maker: In the METSelector Constructor --> getting parameters & calling SetBranches()." << std::endl;
//  metToken_ = iConfig.getParameter<edm::InputTag>("mets");
  metToken_               = iCC.consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
//  puppi_metToken_ = iConfig.getParameter<edm::InputTag>("metsPUPPI");
  puppi_metToken_               = iCC.consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metsPUPPI"));
//  metNoHFToken_ = iConfig.getParameter<edm::InputTag>("metsNoHF");
  metNoHFToken_               = iCC.consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metsNoHF"));
  metCovToken_               = iCC.consumes<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> > >(iConfig.getParameter<edm::InputTag>("metsCov"));
  _is_data   = iConfig.getParameter<bool>("is_data");
  _super_TNT = iConfig.getParameter<bool>("super_TNT");
  SetBranches();
}

METSelector::~METSelector(){
  delete tree_;
}

void METSelector::Fill(const edm::Event& iEvent){
  Clear();
  
  // grab handle to the met collections
  edm::Handle<pat::METCollection> mets;
//  iEvent.getByLabel(metToken_, mets);
  iEvent.getByToken(metToken_, mets);
  const pat::MET &met = mets->front();
  edm::Handle<pat::METCollection> puppimets;
//  iEvent.getByLabel(puppi_metToken_, puppimets);
  iEvent.getByToken(puppi_metToken_, puppimets);
  const pat::MET &puppimet = puppimets->front();
  edm::Handle<pat::METCollection> metsNoHF;
//  iEvent.getByLabel(metNoHFToken_, metsNoHF);
  iEvent.getByToken(metNoHFToken_, metsNoHF);
  const pat::MET &metNoHF = metsNoHF->front();
  edm::Handle<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> > > metCov;
//  iEvent.getByLabel("METSignificance","METCovariance", metCov);
  iEvent.getByToken(metCovToken_, metCov);
  AlgebraicSymMatrix22 metCovMatrix = *(metCov.product());

  if(debug_) std::cout << "     METSelector: Cleared the vectors, grabbed the met collection handles, and extracting met values." << std::endl;
 
  // type-1 PF Met:  propagation of the jet energy corrections to MET. 
  // The jet collection is AK4 CHS jets (non-CHS ak4PFJets) with Pt>10 GeV/c using the L1L2L3-L1 scheme.
  Met_type1PF_pt = met.pt();
  Met_type1PF_px = met.px();
  Met_type1PF_py = met.py();
  Met_type1PF_pz = met.pz();
  Met_type1PF_phi = met.phi();
  Met_type1PF_sumEt = met.sumEt();  
  //Met_type1PF_cov00 = met.getSignificanceMatrix()(0,0);  
  //Met_type1PF_cov01 = met.getSignificanceMatrix()(0,1);  
  //Met_type1PF_cov10 = met.getSignificanceMatrix()(1,0);  
  //Met_type1PF_cov11 = met.getSignificanceMatrix()(1,1);  
  Met_type1PF_cov00 = metCovMatrix[0][0];  
  Met_type1PF_cov01 = metCovMatrix[0][1];  
  Met_type1PF_cov10 = metCovMatrix[1][0];  
  Met_type1PF_cov11 = metCovMatrix[1][1];  

  // PUPPI Met:  Met reconstructed with the PUPPI pileup mitigation algorithm.
  Met_puppi_pt = puppimet.pt();
  Met_puppi_px = puppimet.px();
  Met_puppi_py = puppimet.py();
  Met_puppi_pz = puppimet.pz();
  Met_puppi_phi = puppimet.phi();
  Met_puppi_sumEt = puppimet.sumEt();

  // PUPPI Met:  Met reconstructed with the PUPPI pileup mitigation algorithm.
  Met_NoHF_pt = metNoHF.pt();
  Met_NoHF_px = metNoHF.px();
  Met_NoHF_py = metNoHF.py();
  Met_NoHF_pz = metNoHF.pz();
  Met_NoHF_phi = metNoHF.phi();
  Met_NoHF_sumEt = metNoHF.sumEt();

  // store additional information such as generator level
  if(!_super_TNT){
    if (!_is_data) Gen_Met = met.genMET()->pt();
    Met_type1PF_shiftedPtUp = met.shiftedPt(pat::MET::JetEnUp);
    Met_type1PF_shiftedPtDown  = met.shiftedPt(pat::MET::JetEnDown);
//    Met_puppi_shiftedPtUp = puppimet.shiftedPt(pat::MET::JetEnUp);
//    Met_puppi_shiftedPtDown  = puppimet.shiftedPt(pat::MET::JetEnDown);
  }

}

void METSelector::SetBranches(){
  if(debug_) std::cout << "     METSelector: Setting branches by calling AddBranch of baseTree." << std::endl;
  
  AddBranch(&Met_type1PF_pt,            "Met_type1PF_pt");
  AddBranch(&Met_type1PF_px,            "Met_type1PF_px");
  AddBranch(&Met_type1PF_py,            "Met_type1PF_py");
  AddBranch(&Met_type1PF_pz,            "Met_type1PF_pz");
  AddBranch(&Met_type1PF_sumEt,         "Met_type1PF_sumEt");
  AddBranch(&Met_type1PF_phi,           "Met_type1PF_phi");
  AddBranch(&Met_type1PF_cov00,         "Met_type1PF_cov00");
  AddBranch(&Met_type1PF_cov01,         "Met_type1PF_cov01");
  AddBranch(&Met_type1PF_cov10,         "Met_type1PF_cov10");
  AddBranch(&Met_type1PF_cov11,         "Met_type1PF_cov11");
  AddBranch(&Met_puppi_pt,              "Met_puppi_pt");
  AddBranch(&Met_puppi_px,              "Met_puppi_px");
  AddBranch(&Met_puppi_py,              "Met_puppi_py");
  AddBranch(&Met_puppi_pz,              "Met_puppi_pz");
  AddBranch(&Met_puppi_sumEt,           "Met_puppi_sumEt");
  AddBranch(&Met_puppi_phi,             "Met_puppi_phi");
  AddBranch(&Met_NoHF_pt,               "Met_NoHF_pt");
  AddBranch(&Met_NoHF_px,               "Met_NoHF_px");
  AddBranch(&Met_NoHF_py,               "Met_NoHF_py");
  AddBranch(&Met_NoHF_pz,               "Met_NoHF_pz");
  AddBranch(&Met_NoHF_sumEt,            "Met_NoHF_sumEt");
  AddBranch(&Met_NoHF_phi,              "Met_NoHF_phi");
  if(!_super_TNT){
    AddBranch(&Gen_Met,           "Gen_Met");
    AddBranch(&Met_type1PF_shiftedPtUp,   "Met_type1PF_shiftedPtUp");
    AddBranch(&Met_type1PF_shiftedPtDown, "Met_type1PF_shiftedPtDown");
//    AddBranch(&Met_puppi_shiftedPtUp,     "Met_puppi_shiftedPtUp");
//    AddBranch(&Met_puppi_shiftedPtDown,   "Met_puppi_shiftedPtDown");
  }

  if(debug_) std::cout << "     METSelector: Finished setting branches." << std::endl;
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
  Met_type1PF_cov00  = -9999999999;
  Met_type1PF_cov01  = -9999999999;
  Met_type1PF_cov10  = -9999999999;
  Met_type1PF_cov11  = -9999999999;
  Met_puppi_pt = -9999999999;                                                                                                                                    
  Met_puppi_px = -9999999999;
  Met_puppi_py = -9999999999;
  Met_puppi_pz = -9999999999;
  Met_puppi_phi  = -9999999999;
  Met_puppi_sumEt  = -9999999999 ;
  Met_NoHF_pt = -9999999999;                                                                                                                                    
  Met_NoHF_px = -9999999999;
  Met_NoHF_py = -9999999999;
  Met_NoHF_pz = -9999999999;
  Met_NoHF_phi  = -9999999999;
  Met_NoHF_sumEt  = -9999999999 ;
//  Met_puppi_shiftedPtUp  = -9999999999;
//  Met_puppi_shiftedPtDown   = -9999999999;
  Gen_Met  = -9999999999;
}
