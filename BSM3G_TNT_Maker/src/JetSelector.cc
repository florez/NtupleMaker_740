#include "NtupleMaker/BSM3G_TNT_Maker/interface/JetSelector.h"

JetSelector::JetSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  
  jetToken_       = iConfig.getParameter<edm::InputTag>("jets");
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
  
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByLabel(_vertexInputTag, vtx_h);
  
  reco::VertexCollection::const_iterator firstGoodVertex = vtx_h->end();
  // require a good vertex 
  
  edm::Handle<pat::JetCollection> jets;                                       
  iEvent.getByLabel(jetToken_, jets);                                         
  
  for (const pat::Jet &j : *jets) { 
    if (j.pt() < _Jet_pt_min) continue;

    Jet_pt.push_back(j.pt());         
    Jet_eta.push_back(j.eta());       
    Jet_phi.push_back(j.phi());       
    Jet_energy.push_back(j.energy());
    Jet_bDiscriminator.push_back(j.bDiscriminator("combinedSecondaryVertexBJetTags"));
    Jet_pileupId.push_back(j.userFloat("pileupJetId:fullDiscriminant"));

    if(!_super_TNT){
//      Jet_neutralHadEnergy.push_back(j.neutralHadronEnergy());                               
      Jet_neutralHadEnergyFraction.push_back(j.neutralHadronEnergyFraction());                               
//      Jet_neutralEmEmEnergy.push_back(j.neutralEmEnergy());                                   
      Jet_neutralEmEmEnergyFraction.push_back(j.neutralEmEnergyFraction());                                   
//      Jet_chargedHadronEnergy.push_back(j.chargedHadronEnergy());                               
      Jet_chargedHadronEnergyFraction.push_back(j.chargedHadronEnergyFraction());                               
//      Jet_chargedEmEnergy.push_back(j.chargedEmEnergy());                              
      Jet_chargedEmEnergyFraction.push_back(j.chargedEmEnergyFraction());                              
//      Jet_muonEnergy.push_back(j.muonEnergy());                                  
      Jet_muonEnergyFraction.push_back(j.muonEnergyFraction());                                  

      Jet_numberOfConstituents.push_back(j.chargedMultiplicity() + j.neutralMultiplicity());                                  
      Jet_chargedMultiplicity.push_back(j.chargedMultiplicity());

      Jet_mass.push_back(j.mass()); 
      Jet_electronEnergy.push_back(j.electronEnergy());                               
      Jet_photonEnergy.push_back(j.photonEnergy());                              
      UncorrJet_pt.push_back(j.correctedJet("Uncorrected").pt());                
    }

  } 
  
}

void JetSelector::SetBranches(){
  
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  AddBranch(&Jet_pt,                  "Jet_pt");
  AddBranch(&Jet_eta,                 "Jet_eta");
  AddBranch(&Jet_phi,                 "Jet_phi");
  AddBranch(&Jet_energy,              "Jet_energy");
  AddBranch(&Jet_bDiscriminator,      "Jet_bDiscriminator");
  AddBranch(&Jet_pileupId,            "Jet_pileupId");

  if(!_super_TNT){
    AddBranch(&Jet_mass,                "Jet_mass");
    AddBranch(&Jet_neutralHadEnergyFraction,    "Jet_neutralHadEnergyFraction");
    AddBranch(&Jet_neutralEmEmEnergyFraction,   "Jet_neutralEmEmEnergyFraction");
    AddBranch(&Jet_chargedHadronEnergyFraction, "Jet_chargedHadronEnergyFraction");
    AddBranch(&Jet_chargedEmEnergyFraction,     "Jet_chargedEmEnergyFraction");
    AddBranch(&Jet_muonEnergyFraction,          "Jet_muonEnergyFraction");
    AddBranch(&Jet_numberOfConstituents,"Jet_numberOfConstituents");
    AddBranch(&Jet_chargedMultiplicity, "Jet_chargedMultiplicity");
    AddBranch(&Jet_electronEnergy,      "Jet_electronEnergy");
    AddBranch(&Jet_photonEnergy,        "Jet_photonEnergy");
    AddBranch(&UncorrJet_pt,            "UncorrJet_pt");
  }
  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void JetSelector::Clear(){
  
  Jet_pt.clear();
  Jet_eta.clear();
  Jet_phi.clear();
  Jet_energy.clear();
  Jet_bDiscriminator.clear();
  Jet_pileupId.clear();
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
  
}
