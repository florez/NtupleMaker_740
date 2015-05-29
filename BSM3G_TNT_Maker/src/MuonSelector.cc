#include "NtupleMaker/BSM3G_TNT_Maker/interface/MuonSelector.h"

MuonSelector::MuonSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):baseTree(name,tree,debug){
  
  _muonToken               = iConfig.getParameter<edm::InputTag>("muons");
  _vertexInputTag          = iConfig.getParameter<edm::InputTag>("vertices");
  _Muon_pt_min             = iConfig.getParameter<double>("Muon_pt_min");
  _Muon_eta_max            = iConfig.getParameter<double>("Muon_eta_max");
  _Muon_vtx_ndof_min       = iConfig.getParameter<int>("Muon_vtx_ndof_min");
  _Muon_vtx_rho_max        = iConfig.getParameter<int>("Muon_vtx_rho_max");
  _Muon_vtx_position_z_max = iConfig.getParameter<double>("Muon_vtx_position_z_max");
  _super_TNT               = iConfig.getParameter<bool>("super_TNT");
  SetBranches();
}

MuonSelector::~MuonSelector(){
  delete tree_;
}

void MuonSelector::Fill(const edm::Event& iEvent){
  Clear();
  
  int mucoun = 0;
  
  edm::Handle<edm::View<pat::Muon> > muon_h;
  iEvent.getByLabel(_muonToken, muon_h);
  
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByLabel(_vertexInputTag, vtx_h);
  
  reco::VertexCollection::const_iterator firstGoodVertex = vtx_h->end();
  for (reco::VertexCollection::const_iterator it = vtx_h->begin(); it != firstGoodVertex; it++)
    {
      isGoodVertex(*it);
      firstGoodVertex = it;
      break;
    }
  // require a good vertex 
  if (firstGoodVertex == vtx_h->end()) return;
  
  for(edm::View<pat::Muon>::const_iterator mu = muon_h->begin(); mu != muon_h->end(); mu++){
    if (mu->pt() < _Muon_pt_min) continue;
    if (fabs(mu->eta()) > _Muon_eta_max) continue;  
    
    reco::TrackRef gtk = mu->globalTrack();
    
    Muon_pt.push_back(mu->pt());
    Muon_eta.push_back(mu->eta());
    Muon_phi.push_back(mu->phi());
    Muon_p.push_back(mu->p());
    Muon_charge.push_back(mu->charge());

    if (!_super_TNT){
      if(gtk.isNonnull()){
        Muon_chi2.push_back(gtk->normalizedChi2());
        Muon_validHits.push_back(gtk->hitPattern().numberOfValidMuonHits()); 
      } else {                
        Muon_chi2.push_back(-9999);
        Muon_validHits.push_back(-9999);  
      }
      Muon_matchedStat.push_back(mu->numberOfMatchedStations());
      Muon_dxy.push_back(mu->muonBestTrack()->dxy(firstGoodVertex->position())); 
    
      if(mu->innerTrack().isNonnull()){
        Muon_validHitsInner.push_back(mu->innerTrack()->hitPattern().numberOfValidPixelHits());
        Muon_TLayers.push_back(mu->innerTrack()->hitPattern().trackerLayersWithMeasurement());
      } else {
        Muon_validHitsInner.push_back(-999);
        Muon_TLayers.push_back(-999);
      }
    
      Muon_dz.push_back(mu->muonBestTrack()->dz(firstGoodVertex->position()));
      Muon_isoNeutralHadron.push_back((mu->pfIsolationR04().sumNeutralHadronEt));
      Muon_isoPhoton.push_back((mu->pfIsolationR04().sumPhotonEt));
      Muon_isoPU.push_back((mu->pfIsolationR04().sumPUPt));
    }
    Muon_loose.push_back(mu->isLooseMuon());
    Muon_tight.push_back(mu->isTightMuon(*firstGoodVertex));
    Muon_soft.push_back(mu->isSoftMuon(*firstGoodVertex));
    Muon_isHighPt.push_back(mu->isHighPtMuon(*firstGoodVertex));
    Muon_pf.push_back(mu->isPFMuon());   
    Muon_energy.push_back(mu->energy());
    Muon_isoSum.push_back((mu->trackIso() + mu->ecalIso() + mu->hcalIso()));
    Muon_isoCharParPt.push_back((mu->pfIsolationR04().sumChargedParticlePt));
    Muon_isoCharged.push_back((mu->pfIsolationR04().sumChargedHadronPt));
    
    mucoun++;
  }
  
  
}

void MuonSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;

  AddBranch(&Muon_pt                ,"Muon_pt");
  AddBranch(&Muon_eta               ,"Muon_eta");
  AddBranch(&Muon_phi               ,"Muon_phi");
  AddBranch(&Muon_p                 ,"Muon_p");
  AddBranch(&Muon_energy            ,"Muon_energy");
  AddBranch(&Muon_charge            ,"Muon_charge");
  AddBranch(&Muon_tight             ,"Muon_tight");
  AddBranch(&Muon_loose             ,"Muon_loose");
  AddBranch(&Muon_soft              ,"Muon_soft");
  AddBranch(&Muon_isHighPt          ,"Muon_isHighPt");
  AddBranch(&Muon_pf                ,"Muon_pf");
  AddBranch(&Muon_isoCharged        ,"Muon_isoCharged");
  AddBranch(&Muon_isoSum            ,"Muon_isoSum");
  AddBranch(&Muon_isoCharParPt      ,"Muon_isoCharParPt");
  if (!_super_TNT){  
    AddBranch(&Muon_chi2              ,"Muon_chi2");
    AddBranch(&Muon_validHits         ,"Muon_validHits");
    AddBranch(&Muon_validHitsInner    ,"Muon_validHitsInner");
    AddBranch(&Muon_matchedStat       ,"Muon_matchedStat");
    AddBranch(&Muon_dxy               ,"Muon_dxy");
    AddBranch(&Muon_TLayers           ,"Muon_TLayers");
    AddBranch(&Muon_dz                ,"Muon_dz");
    AddBranch(&Muon_isoNeutralHadron  ,"Muon_isoNeutralHadron");
    AddBranch(&Muon_isoPhoton         ,"Muon_isoPhoton");
    AddBranch(&Muon_isoPU             ,"Muon_isoPU");
  }
  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void MuonSelector::Clear(){
  
  Muon_pt.clear();
  Muon_eta.clear();
  Muon_phi.clear();
  Muon_energy.clear();
  Muon_dz.clear();
  Muon_tight.clear();
  Muon_loose.clear();
  Muon_soft.clear();
  Muon_pf.clear();   
  Muon_isHighPt.clear();
  Muon_isoCharged.clear();
  Muon_isoNeutralHadron.clear();
  Muon_isoPhoton.clear();
  Muon_isoPU.clear();
  Muon_p.clear(); 
  Muon_charge.clear(); 
  Muon_chi2.clear(); 
  Muon_validHits.clear();
  Muon_validHitsInner.clear(); 
  Muon_matchedStat.clear(); 
  Muon_dxy.clear(); 
  Muon_TLayers.clear(); 
  
}

bool MuonSelector::isGoodVertex(const reco::Vertex& vtx)
{
  if (vtx.isFake()) return false;
  if (vtx.ndof() < _Muon_vtx_ndof_min) return false;
  if (vtx.position().Rho() > _Muon_vtx_rho_max) return false;
  if (fabs(vtx.position().Z()) > _Muon_vtx_position_z_max) return false;
  return true;
}
