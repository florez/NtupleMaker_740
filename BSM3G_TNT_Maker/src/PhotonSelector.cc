#include "NtupleMaker/BSM3G_TNT_Maker/interface/PhotonSelector.h"

PhotonSelector::PhotonSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector&& iCC): baseTree(name,tree,debug){
  if(debug) std::cout << "BSM3G TNT Maker: In the PhotonSelector Constructor --> getting parameters & calling SetBranches()." << std::endl;
//  _PhotonToken      = iConfig.getParameter<edm::InputTag>("photons");
  _PhotonToken      = iCC.consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photons"));
  _Photon_pt_min    = iConfig.getParameter<double>("Photon_pt_min");
  _Photon_eta_max   = iConfig.getParameter<double>("Photon_eta_max");
  SetBranches();
}

PhotonSelector::~PhotonSelector(){
  delete tree_;
}

void PhotonSelector::Fill(const edm::Event& iEvent){
  Clear();
  
  edm::Handle<edm::View<pat::Photon> > photon_handle;
//  iEvent.getByLabel(_PhotonToken, photon_handle);
  iEvent.getByToken(_PhotonToken, photon_handle);

  if(debug_) std::cout << "     PhotonSelector: Cleared the vectors, grabbed the photon collection handle, and looping over photons." << std::endl;
 
  if (photon_handle.isValid()) { 
    for(edm::View<pat::Photon>::const_iterator ph = photon_handle->begin(); ph != photon_handle->end(); ph++){
      if (ph->pt() < _Photon_pt_min) continue;
      if (fabs(ph->eta()) > _Photon_eta_max) continue;  
      Photon_pt.push_back(ph->pt());
      Photon_eta.push_back(ph->eta());
      Photon_phi.push_back(ph->phi());
      Photon_energy.push_back(ph->energy());
      Photon_et.push_back(ph->et());
      Photon_HoverE.push_back(ph->hadTowOverEm());
      Photon_phoR9.push_back(ph->r9());
      Photon_SigmaIEtaIEta.push_back(ph->see());
      Photon_SigmaIPhiIPhi.push_back(ph->spp());
      Photon_PFChIso.push_back(ph->chargedHadronIso());
      Photon_PFPhoIso.push_back(ph->photonIso());
      Photon_PFNeuIso.push_back(ph->neutralHadronIso());
      Photon_EleVeto.push_back((int)ph->passElectronVeto());
      Photon_hasPixelSeed.push_back((int)ph->hasPixelSeed());
    }
  }
}

void PhotonSelector::SetBranches(){
  if(debug_) std::cout << "     PhotonSelector: Setting branches by calling AddBranch of baseTree." << std::endl;
  
  AddBranch(&Photon_pt             ,"Photon_pt");
  AddBranch(&Photon_eta            ,"Photon_eta");
  AddBranch(&Photon_phi            ,"Photon_phi");
  AddBranch(&Photon_energy         ,"Photon_energy");
  AddBranch(&Photon_et             ,"Photon_et");
  AddBranch(&Photon_HoverE         ,"Photon_HoverE");
  AddBranch(&Photon_phoR9          ,"Photon_phoR9");
  AddBranch(&Photon_SigmaIEtaIEta  ,"Photon_SigmaIEtaIEta");
  AddBranch(&Photon_SigmaIPhiIPhi  ,"Photon_SigmaIPhiIPhi");
  AddBranch(&Photon_PFChIso        ,"Photon_PFChIso");
  AddBranch(&Photon_PFPhoIso       ,"Photon_PFPhoIso");
  AddBranch(&Photon_PFNeuIso       ,"Photon_PFNeuIso");
  AddBranch(&Photon_EleVeto        ,"Photon_EleVeto");
  AddBranch(&Photon_hasPixelSeed   ,"Photon_hasPixelSeed");

  if(debug_) std::cout << "     PhotonSelector: Finished setting branches." << std::endl;
}

void PhotonSelector::Clear(){
  
  Photon_pt.clear();
  Photon_eta.clear();
  Photon_phi.clear();
  Photon_energy.clear();
  Photon_et.clear();
  Photon_HoverE.clear();
  Photon_phoR9.clear();
  Photon_SigmaIEtaIEta.clear();
  Photon_SigmaIPhiIPhi.clear();
  Photon_PFChIso.clear();
  Photon_PFPhoIso.clear();
  Photon_PFNeuIso.clear();
  Photon_EleVeto.clear();
  Photon_hasPixelSeed.clear();
}
