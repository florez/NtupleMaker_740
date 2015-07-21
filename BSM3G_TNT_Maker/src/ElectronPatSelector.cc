#include "NtupleMaker/BSM3G_TNT_Maker/interface/ElectronPatSelector.h"

ElectronPatSelector::ElectronPatSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic): 
  baseTree(name,tree,debug),
  electronVetoIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"))),
  electronLooseIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronLooseIdMap"))),
  electronMediumIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMediumIdMap"))),
  electronTightIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap"))),
  eleHEEPIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap")))
{
  
  SetBranches();
  _vertexInputTag        = iConfig.getParameter<edm::InputTag>("vertices");
  _beamSpot              = iConfig.getParameter<edm::InputTag>("beamSpot");
  _patElectronToken      = iConfig.getParameter<edm::InputTag>("patElectrons");
  _patElectron_pt_min    = iConfig.getParameter<double>("patElectron_pt_min");
  _patElectron_eta_max   = iConfig.getParameter<double>("patElectron_eta_max");
}

ElectronPatSelector::~ElectronPatSelector(){
  delete tree_;
}

void ElectronPatSelector::Fill(const edm::Event& iEvent){
  Clear();
  
  int elcoun = 0;
  
  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByLabel(_vertexInputTag, vtx);
  if (vtx->empty()) return;
  
  edm::Handle<edm::View<pat::Electron> > electron_pat;
  iEvent.getByLabel(_patElectronToken, electron_pat);
  
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
  
  iEvent.getByToken(electronVetoIdMapToken_,veto_id_decisions);
  iEvent.getByToken(electronLooseIdMapToken_,loose_id_decisions);
  iEvent.getByToken(electronMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(electronTightIdMapToken_,tight_id_decisions);  
  iEvent.getByToken(eleHEEPIdMapToken_, heep_id_decisions);

  // Get BeamSpot Information
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(_beamSpot, beamSpotHandle);

  if ( beamSpotHandle.isValid() )
   {
     beamSpot = *beamSpotHandle;

   } else {
    edm::LogInfo("MyAnalyzer")
      << "No beam spot available from EventSetup \n";
   }

  math::XYZPoint point(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());


  reco::VertexCollection::const_iterator firstGoodVertex = vtx->end();
  for(edm::View<pat::Electron>::const_iterator el = electron_pat->begin(); el != electron_pat->end(); el++){
    if (el->pt() < _patElectron_pt_min) continue;
    if (fabs(el->eta()) > _patElectron_eta_max) continue;  
    
    patElectron_pt.push_back(el->pt());
    patElectron_eta.push_back(el->eta());
    patElectron_phi.push_back(el->phi());
    patElectron_charge.push_back(el->charge());
    patElectron_energy.push_back(el->energy());
    
    // Look up the ID decision for this electron in
    // the ValueMap object and store it. We need a Ptr object as the key.
    const Ptr<pat::Electron> elPtr(electron_pat, el - electron_pat->begin() );
    bool isPassVeto   = (*veto_id_decisions)[ elPtr ];
    bool isPassLoose  = (*loose_id_decisions)[ elPtr ];
    bool isPassMedium = (*medium_id_decisions)[ elPtr ];
    bool isPassTight  = (*tight_id_decisions)[ elPtr ];
    bool isHEEPId     = (*heep_id_decisions)[ elPtr ];
 
    passVetoId_.push_back( isPassVeto );
    passLooseId_.push_back( isPassLoose );
    passMediumId_.push_back( isPassMedium );
    passTightId_.push_back( isPassTight );
    passHEEPId_.push_back( isHEEPId );   

    // Variables for life-time (impact parameter) studies  
    patElectron_gsfTrack_normChi2.push_back(el->gsfTrack()->normalizedChi2());
    patElectron_gsfTrack_ndof.push_back(el->gsfTrack()->ndof());
    patElectron_gsfTrack_vtx.push_back(el->gsfTrack()->vx()); 
    patElectron_gsfTrack_vty.push_back(el->gsfTrack()->vy());
    patElectron_gsfTrack_vtz.push_back(el->gsfTrack()->vz());
    patElectron_dxy.push_back( (-1) * el->gsfTrack()->dxy(firstGoodVertex->position() ) );
    patElectron_dxyError.push_back(el->gsfTrack()->d0Error());
    patElectron_dz.push_back( el->gsfTrack()->dz( firstGoodVertex->position() ) );
    patElectron_dxy_bs.push_back ( el->gsfTrack()->dxy(point) );
    expectedMissingInnerHits.push_back(el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) );
    passConversionVeto_.push_back( el->passConversionVeto() );
    
    // Isolation
    reco::GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();
    // Compute isolation with delta beta correction for PU
    isoChargedHadrons_.push_back( pfIso.sumChargedHadronPt );
    isoNeutralHadrons_.push_back( pfIso.sumNeutralHadronEt );
    isoPhotons_.push_back( pfIso.sumPhotonEt );
    isoPU_.push_back( pfIso.sumPUPt );
    
    elcoun++;
  }
  
}

void ElectronPatSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;
  
  AddBranch(&patElectron_pt                ,"patElectron_pt");
  AddBranch(&patElectron_eta               ,"patElectron_eta");
  AddBranch(&patElectron_phi               ,"patElectron_phi");
  AddBranch(&patElectron_energy            ,"patElectron_energy");
  AddBranch(&patElectron_charge            ,"patElectron_charge");
  AddBranch(&passVetoId_                   ,"patElectron_isPassVeto");          
  AddBranch(&passLooseId_                  ,"patElectron_isPassLoose");
  AddBranch(&passMediumId_                 ,"patElectron_isPassMedium");
  AddBranch(&passTightId_                  ,"patElectron_isPassTight");
  AddBranch(&passHEEPId_                   ,"patElectron_isPassHEEPId");
  AddBranch(&patElectron_dxy               ,"patElectron_dxy");
  AddBranch(&patElectron_dxyError          ,"patElectron_dxyError");
  AddBranch(&patElectron_dz                ,"patElectron_dz");
  AddBranch(&patElectron_dxy_bs            ,"patElectron_dxy_bs");
  AddBranch(&patElectron_gsfTrack_normChi2 ,"patElectron_gsfTrack_normChi2");
  AddBranch(&patElectron_gsfTrack_ndof     ,"patElectron_gsfTrack_ndof");
  AddBranch(&patElectron_gsfTrack_vtx      ,"patElectron_gsfTrack_vtx");
  AddBranch(&patElectron_gsfTrack_vty      ,"patElectron_gsfTrack_vty");
  AddBranch(&patElectron_gsfTrack_vtz      ,"patElectron_gsfTrack_vtz");
  AddBranch(&expectedMissingInnerHits      ,"patElectron_expectedMissingInnerHits");
  AddBranch(&passConversionVeto_           ,"patElectron_passConversionVeto"); 
  AddBranch(&isoChargedHadrons_            ,"patElectron_isoChargedHadrons");
  AddBranch(&isoNeutralHadrons_            ,"patElectron_isoNeutralHadrons");
  AddBranch(&isoPhotons_                   ,"patElectron_isoPhotons");
  AddBranch(&isoPU_                        ,"patElectron_isoPU");
  
  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void ElectronPatSelector::Clear(){
  
  patElectron_pt.clear();
  patElectron_eta.clear();
  patElectron_phi.clear();
  patElectron_energy.clear();
  patElectron_charge.clear(); 
  passVetoId_.clear();
  passLooseId_.clear();
  passMediumId_.clear();
  passTightId_.clear();  
  passHEEPId_.clear();
  patElectron_dxy.clear();
  patElectron_dxyError.clear();
  patElectron_gsfTrack_normChi2.clear();
  patElectron_gsfTrack_ndof.clear();
  patElectron_gsfTrack_vtx.clear(); 
  patElectron_gsfTrack_vty.clear();
  patElectron_gsfTrack_vtz.clear();
  patElectron_dz.clear();
  patElectron_dxy_bs.clear();
  expectedMissingInnerHits.clear();
  passConversionVeto_.clear();
  isoChargedHadrons_.clear();
  isoNeutralHadrons_.clear();
  isoPhotons_.clear();
  isoPU_.clear();
}
