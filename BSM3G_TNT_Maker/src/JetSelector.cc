// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#include "NtupleMaker/BSM3G_TNT_Maker/interface/JetSelector.h"

JetSelector::JetSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  if(debug) std::cout << "BSM3G TNT Maker: In the JetSelector Constructor --> getting parameters & calling SetBranches()." << std::endl;
  jetToken_       = iConfig.getParameter<edm::InputTag>("jets");
  puppi_jetToken_ = iConfig.getParameter<edm::InputTag>("jetsPUPPI");
  _vertexInputTag = iConfig.getParameter<edm::InputTag>("vertices");
  _Jet_pt_min     = iConfig.getParameter<double>("Jet_pt_min");
  _super_TNT      = iConfig.getParameter<bool>("super_TNT");
  SetBranches();
}

JetSelector::~JetSelector(){
  delete tree_;
}

void JetSelector::Fill(const edm::Event& iEvent){
  Clear();
  
  // grab handle to the jet collections
  edm::Handle<pat::JetCollection> jets;                                       
  iEvent.getByLabel(jetToken_, jets);                                         
  edm::Handle<pat::JetCollection> puppijets;                                       
  iEvent.getByLabel(puppi_jetToken_, puppijets);                                         

  if(debug_) std::cout << "     JetSelector: Cleared the vectors, grabbed the jet collection handle, and looping over jets." << std::endl;
  
  // loop over "standard" ak4 CHS PF jets
  for (const pat::Jet &j : *jets) { 
    if (j.pt() < _Jet_pt_min) continue; // only keep jets passing some user defined pt cut (defined in miniAOD.py)

    // fill root tree with "necessary" information:  kinematics, ID discriminators
    Jet_pt.push_back(j.pt());         
    Jet_eta.push_back(j.eta());       
    Jet_phi.push_back(j.phi());       
    Jet_energy.push_back(j.energy());
    Jet_bDiscriminator.push_back(j.bDiscriminator("combinedSecondaryVertexBJetTags"));
    Jet_bDiscriminator_CISVV2.push_back(j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"));
    Jet_bDiscriminator_pfCISVV2.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    Jet_vtxMass.push_back(j.userFloat("vtxMass"));
    Jet_decayLength.push_back(j.userFloat("vtx3DVal"));
    Jet_decayLengthSignificance.push_back(j.userFloat("vtx3DSig"));
    Jet_pileupId.push_back(j.userFloat("pileupJetId:fullDiscriminant"));
    Jet_partonFlavour.push_back(j.partonFlavour());
    Jet_neutralHadEnergyFraction.push_back(j.neutralHadronEnergyFraction());                               
    Jet_neutralEmEmEnergyFraction.push_back(j.neutralEmEnergyFraction());                                   
    Jet_chargedHadronEnergyFraction.push_back(j.chargedHadronEnergyFraction());                               
    Jet_chargedEmEnergyFraction.push_back(j.chargedEmEnergyFraction());                              
    Jet_muonEnergyFraction.push_back(j.muonEnergyFraction());                                  
    Jet_numberOfConstituents.push_back(j.chargedMultiplicity() + j.neutralMultiplicity());                                  
    Jet_chargedMultiplicity.push_back(j.chargedMultiplicity());

    // store additional information such as the "raw" uncorrected pt
    if(!_super_TNT){
      Jet_mass.push_back(j.mass()); 
      Jet_electronEnergy.push_back(j.electronEnergy());                               
      Jet_photonEnergy.push_back(j.photonEnergy());                              
      UncorrJet_pt.push_back(j.correctedJet("Uncorrected").pt());                
    }
  } 
  
  // loop over ak4 PF jets reconstructed with the PUPPI pileup mitigation algorithms
  for (const pat::Jet &j : *puppijets) { 
    if (j.pt() < _Jet_pt_min) continue; // only keep jets passing some user defined pt cut (defined in miniAOD.py)

    // fill root tree with "necessary" information:  kinematics, ID discriminators
    Jet_puppi_pt.push_back(j.pt());         
    Jet_puppi_eta.push_back(j.eta());       
    Jet_puppi_phi.push_back(j.phi());       
    Jet_puppi_energy.push_back(j.energy());
    Jet_puppi_bDiscriminator.push_back(j.bDiscriminator("combinedSecondaryVertexBJetTags"));
    Jet_puppi_bDiscriminator_CISVV2.push_back(j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"));
    Jet_puppi_bDiscriminator_pfCISVV2.push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    Jet_puppi_vtxMass.push_back(j.userFloat("vtxMass"));
    Jet_puppi_decayLength.push_back(j.userFloat("vtx3DVal"));
    Jet_puppi_decayLengthSignificance.push_back(j.userFloat("vtx3DSig"));
    Jet_puppi_pileupId.push_back(j.userFloat("pileupJetId:fullDiscriminant"));
    Jet_puppi_partonFlavour.push_back(j.partonFlavour());
    Jet_puppi_neutralHadEnergyFraction.push_back(j.neutralHadronEnergyFraction());                               
    Jet_puppi_neutralEmEmEnergyFraction.push_back(j.neutralEmEnergyFraction());                                   
    Jet_puppi_chargedHadronEnergyFraction.push_back(j.chargedHadronEnergyFraction());                               
    Jet_puppi_chargedEmEnergyFraction.push_back(j.chargedEmEnergyFraction());                              
    Jet_puppi_muonEnergyFraction.push_back(j.muonEnergyFraction());                                  
    Jet_puppi_numberOfConstituents.push_back(j.chargedMultiplicity() + j.neutralMultiplicity());                                  
    Jet_puppi_chargedMultiplicity.push_back(j.chargedMultiplicity());

    // store additional information such as the "raw" uncorrected pt
    if(!_super_TNT){
      Jet_puppi_mass.push_back(j.mass()); 
      Jet_puppi_electronEnergy.push_back(j.electronEnergy());                               
      Jet_puppi_photonEnergy.push_back(j.photonEnergy());                              
      UncorrJet_puppi_pt.push_back(j.correctedJet("Uncorrected").pt());                
    }
  } 
}

void JetSelector::SetBranches(){
  if(debug_) std::cout << "     JetSelector: Setting branches by calling AddBranch of baseTree." << std::endl;

  AddBranch(&Jet_pt,                  		"Jet_pt");
  AddBranch(&Jet_eta,                 		"Jet_eta");
  AddBranch(&Jet_phi,                 		"Jet_phi");
  AddBranch(&Jet_energy,             		"Jet_energy");
  AddBranch(&Jet_bDiscriminator,      		"Jet_bDiscriminator");
  AddBranch(&Jet_bDiscriminator_CISVV2,      	"Jet_bDiscriminator_CISVV2");
  AddBranch(&Jet_bDiscriminator_pfCISVV2,      	"Jet_bDiscriminator_pfCISVV2");
  AddBranch(&Jet_vtxMass,            		"Jet_vtxMass");
  AddBranch(&Jet_decayLength,            	"Jet_decayLength");
  AddBranch(&Jet_decayLengthSignificance,	"Jet_decayLengthSignificance");
  AddBranch(&Jet_pileupId,            		"Jet_pileupId");
  AddBranch(&Jet_partonFlavour,       		"Jet_partonFlavour");
  AddBranch(&Jet_neutralHadEnergyFraction,    	"Jet_neutralHadEnergyFraction");
  AddBranch(&Jet_neutralEmEmEnergyFraction,   	"Jet_neutralEmEmEnergyFraction");
  AddBranch(&Jet_chargedHadronEnergyFraction, 	"Jet_chargedHadronEnergyFraction");
  AddBranch(&Jet_chargedEmEnergyFraction,     	"Jet_chargedEmEnergyFraction");
  AddBranch(&Jet_muonEnergyFraction,          	"Jet_muonEnergyFraction");
  AddBranch(&Jet_numberOfConstituents,		"Jet_numberOfConstituents");
  AddBranch(&Jet_chargedMultiplicity, 		"Jet_chargedMultiplicity");
  AddBranch(&Jet_puppi_pt,             		"Jet_puppi_pt");
  AddBranch(&Jet_puppi_eta,              	"Jet_puppi_eta");
  AddBranch(&Jet_puppi_phi,                 	"Jet_puppi_phi");
  AddBranch(&Jet_puppi_energy,             	"Jet_puppi_energy");
  AddBranch(&Jet_puppi_bDiscriminator,      	"Jet_puppi_bDiscriminator");
  AddBranch(&Jet_puppi_bDiscriminator_CISVV2,   "Jet_puppi_bDiscriminator_CISVV2");
  AddBranch(&Jet_puppi_bDiscriminator_pfCISVV2, "Jet_puppi_bDiscriminator_pfCISVV2");
  AddBranch(&Jet_puppi_vtxMass,       		"Jet_puppi_vtxMass");
  AddBranch(&Jet_puppi_decayLength,            	"Jet_puppi_decayLength");
  AddBranch(&Jet_puppi_decayLengthSignificance,	"Jet_puppi_decayLengthSignificance");
  AddBranch(&Jet_puppi_pileupId,            	"Jet_puppi_pileupId");
  AddBranch(&Jet_puppi_partonFlavour,       	"Jet_puppi_partonFlavour");
  AddBranch(&Jet_puppi_neutralHadEnergyFraction,"Jet_puppi_neutralHadEnergyFraction");
  AddBranch(&Jet_puppi_neutralEmEmEnergyFraction,"Jet_puppi_neutralEmEmEnergyFraction");
  AddBranch(&Jet_puppi_chargedHadronEnergyFraction,"Jet_puppi_chargedHadronEnergyFraction");
  AddBranch(&Jet_puppi_chargedEmEnergyFraction, "Jet_puppi_chargedEmEnergyFraction");
  AddBranch(&Jet_puppi_muonEnergyFraction,      "Jet_puppi_muonEnergyFraction");
  AddBranch(&Jet_puppi_numberOfConstituents,	"Jet_puppi_numberOfConstituents");
  AddBranch(&Jet_puppi_chargedMultiplicity, 	"Jet_puppi_chargedMultiplicity");
  if(!_super_TNT){
    AddBranch(&Jet_mass,               		"Jet_mass");
    AddBranch(&Jet_electronEnergy,      	"Jet_electronEnergy");
    AddBranch(&Jet_photonEnergy,        	"Jet_photonEnergy");
    AddBranch(&UncorrJet_pt,            	"UncorrJet_pt");
    AddBranch(&Jet_puppi_mass,         		"Jet_puppi_mass");
    AddBranch(&Jet_puppi_electronEnergy,      	"Jet_puppi_electronEnergy");
    AddBranch(&Jet_puppi_photonEnergy,        	"Jet_puppi_photonEnergy");
    AddBranch(&UncorrJet_puppi_pt,            	"UncorrJet_puppi_pt");
  }

  if(debug_) std::cout << "     JetSelector: Finished setting branches." << std::endl;
}

void JetSelector::Clear(){
  
  Jet_pt.clear();
  Jet_eta.clear();
  Jet_phi.clear();
  Jet_energy.clear();
  Jet_bDiscriminator.clear();
  Jet_bDiscriminator_CISVV2.clear();
  Jet_bDiscriminator_pfCISVV2.clear();
  Jet_vtxMass.clear();
  Jet_decayLength.clear();
  Jet_decayLengthSignificance.clear();
  Jet_pileupId.clear();
  Jet_partonFlavour.clear();
  Jet_mass.clear();
  Jet_neutralHadEnergyFraction.clear();
  Jet_neutralEmEmEnergyFraction.clear();
  Jet_chargedHadronEnergyFraction.clear();
  Jet_chargedEmEnergyFraction.clear();
  Jet_muonEnergyFraction.clear();
  Jet_numberOfConstituents.clear();
  Jet_chargedMultiplicity.clear();
  Jet_electronEnergy.clear();
  Jet_photonEnergy.clear();
  UncorrJet_pt.clear();
  Jet_puppi_pt.clear();
  Jet_puppi_eta.clear();
  Jet_puppi_phi.clear();
  Jet_puppi_energy.clear();
  Jet_puppi_bDiscriminator.clear();
  Jet_puppi_bDiscriminator_CISVV2.clear();
  Jet_puppi_bDiscriminator_pfCISVV2.clear();
  Jet_puppi_vtxMass.clear();
  Jet_puppi_decayLength.clear();
  Jet_puppi_decayLengthSignificance.clear();
  Jet_puppi_pileupId.clear();
  Jet_puppi_partonFlavour.clear();
  Jet_puppi_mass.clear();
  Jet_puppi_neutralHadEnergyFraction.clear();
  Jet_puppi_neutralEmEmEnergyFraction.clear();
  Jet_puppi_chargedHadronEnergyFraction.clear();
  Jet_puppi_chargedEmEnergyFraction.clear();
  Jet_puppi_muonEnergyFraction.clear();
  Jet_puppi_numberOfConstituents.clear();
  Jet_puppi_chargedMultiplicity.clear();
  Jet_puppi_electronEnergy.clear();
  Jet_puppi_photonEnergy.clear();
  UncorrJet_puppi_pt.clear();
  
}
