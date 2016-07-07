// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#include "NtupleMaker/BSM3G_TNT_Maker/interface/ElectronPatSelector.h"

using namespace edm;
using namespace std;
using namespace reco;

ElectronPatSelector::ElectronPatSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic): 
  baseTree(name,tree,debug),
  electronVetoIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"))),
  electronLooseIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronLooseIdMap"))),
  electronMediumIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMediumIdMap"))),
  electronTightIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap"))),
  eleHEEPIdMapToken_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap"))),
  electronMVAwp1Token_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMVA_wp1"))),
  electronMVAwp2Token_(ic.consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMVA_wp2")))
{
  if(debug) std::cout << "BSM3G TNT Maker: In the ElectronPatSelector Constructor --> getting parameters & calling SetBranches()." << std::endl;
  _patElectronToken      	  = iConfig.getParameter<edm::InputTag>("patElectrons");
  _vertexInputTag       	  = iConfig.getParameter<edm::InputTag>("vertices");
  _beamSpot              	  = iConfig.getParameter<edm::InputTag>("beamSpot");
  _patElectron_pt_min    	  = iConfig.getParameter<double>("patElectron_pt_min");
  _patElectron_eta_max   	  = iConfig.getParameter<double>("patElectron_eta_max");
  _patElectron_vtx_ndof_min       = iConfig.getParameter<int>("patElectron_vtx_ndof_min");
  _patElectron_vtx_rho_max        = iConfig.getParameter<int>("patElectron_vtx_rho_max");
  _patElectron_vtx_position_z_max = iConfig.getParameter<double>("patElectron_vtx_position_z_max");
  _super_TNT               	  = iConfig.getParameter<bool>("super_TNT");
  SetBranches();

}

ElectronPatSelector::~ElectronPatSelector(){
  delete tree_;
}

void ElectronPatSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();

  if(debug_) std::cout << "     ElectronPatSelector: Cleared the vectors, grabbing the electron/PV/BeamSpot collection handles." << std::endl;
  
  int elcoun = 0;
  
  // grab handle to the electron collection
  edm::Handle<edm::View<pat::Electron> > electron_pat;
  iEvent.getByLabel(_patElectronToken, electron_pat);

  if(debug_) std::cout << "     ElectronPatSelector: Looping over PVs to extract the position of the best PV." << std::endl;

  // grab handle to the vertex collection
  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByLabel(_vertexInputTag, vtx);
  reco::VertexCollection::const_iterator firstGoodVertex = vtx->end();
  for (reco::VertexCollection::const_iterator it = vtx->begin(); it != vtx->end(); it++) {
    if (isGoodVertex(*it)) {
      firstGoodVertex = it;
      break;
    }
  }
  if (firstGoodVertex == vtx->end()) return; // skip event if there are no good PVs
  GlobalPoint thepv(firstGoodVertex->position().x(),firstGoodVertex->position().y(),firstGoodVertex->position().z());

  if(debug_) std::cout << "     ElectronPatSelector: Extracting the beamspot position." << std::endl;

  // Get BeamSpot Information
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(_beamSpot, beamSpotHandle);
  if ( beamSpotHandle.isValid() ) { beamSpot = *beamSpotHandle; }
  else { edm::LogInfo("MyAnalyzer") << "No beam spot available from EventSetup \n"; }
  math::XYZPoint point(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());
  GlobalPoint thebs(beamSpot.x0(),beamSpot.y0(),beamSpot.z0());

  if(debug_) std::cout << "     ElectronPatSelector: Grabbing the electron ID value maps." << std::endl;

  // Get electron ID value maps
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
  edm::Handle<edm::ValueMap<bool> > mva_wp1_id_decisions;
  edm::Handle<edm::ValueMap<bool> > mva_wp2_id_decisions;
  iEvent.getByToken(electronVetoIdMapToken_,veto_id_decisions);
  iEvent.getByToken(electronLooseIdMapToken_,loose_id_decisions);
  iEvent.getByToken(electronMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(electronTightIdMapToken_,tight_id_decisions);  
  iEvent.getByToken(eleHEEPIdMapToken_, heep_id_decisions);
  iEvent.getByToken(electronMVAwp1Token_, mva_wp1_id_decisions);
  iEvent.getByToken(electronMVAwp2Token_, mva_wp2_id_decisions);

  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  if(debug_) std::cout << "     ElectronPatSelector: Looping over electrons." << std::endl;

  // loop over electrons
  int nElec=0;
  for(edm::View<pat::Electron>::const_iterator el = electron_pat->begin(); el != electron_pat->end(); el++) {
    nElec++;

    if(debug_) std::cout << "          Electron info: Elec #" << nElec << ": pt = " << el->pt() << ", eta = " << el->eta() << ", gsfTrackIsNonnull = " << el->gsfTrack().isNonnull() << ", trackIsNonnull = " << el->closestCtfTrackRef().isNonnull() << std::endl;

    if (el->pt() < _patElectron_pt_min) continue;
    if (fabs(el->eta()) > _patElectron_eta_max) continue;  
    if (!(el->gsfTrack().isNonnull())) continue;
    if (!(el->closestCtfTrackRef().isNonnull())) continue;

    if(debug_) std::cout << "          Electron info: Elec #" << nElec << " passed kinematic cuts & isNonNull cuts" << std::endl;

    // fill root tree with "necessary" information:  kinematics, ID, isolation
    patElectron_pt.push_back(el->pt());
    patElectron_eta.push_back(el->eta());
    patElectron_phi.push_back(el->phi());
    patElectron_energy.push_back(el->energy());
    patElectron_charge.push_back(el->charge());
    
    // Look up the ID decision for this electron in
    // the ValueMap object and store it. We need a Ptr object as the key.
    const Ptr<pat::Electron> elPtr(electron_pat, el - electron_pat->begin() );
    bool isPassVeto   = (*veto_id_decisions)[ elPtr ];
    bool isPassLoose  = (*loose_id_decisions)[ elPtr ];
    bool isPassMedium = (*medium_id_decisions)[ elPtr ];
    bool isPassTight  = (*tight_id_decisions)[ elPtr ];
    bool isHEEPId     = (*heep_id_decisions)[ elPtr ];
    bool isMVAwp1     = (*mva_wp1_id_decisions)[ elPtr ];
    bool isMVAwp2     = (*mva_wp2_id_decisions)[ elPtr ]; 

    passVetoId_.push_back( isPassVeto );
    passLooseId_.push_back( isPassLoose );
    passMediumId_.push_back( isPassMedium );
    passTightId_.push_back( isPassTight );
    passHEEPId_.push_back( isHEEPId );   
    passMVAwp1Id_.push_back(isMVAwp1);
    passMVAwp2Id_.push_back(isMVAwp2);
    // Isolation
    reco::GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();
    // Compute isolation with delta beta correction for PU
    isoChargedHadrons_.push_back( pfIso.sumChargedHadronPt );
    isoNeutralHadrons_.push_back( pfIso.sumNeutralHadronEt );
    isoPhotons_.push_back( pfIso.sumPhotonEt );
    isoPU_.push_back( pfIso.sumPUPt );

    // conversion veto/info
    passConversionVeto_.push_back( el->passConversionVeto() );
    expectedMissingInnerHits.push_back(el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) );

    // store additional information such as lifetime variables (needed for tau analyses)
    if (!_super_TNT){

      // Variables for life-time (impact parameter) studies  
      patElectron_gsfTrack_normChi2.push_back(el->gsfTrack()->normalizedChi2());
      patElectron_gsfTrack_ndof.push_back(el->gsfTrack()->ndof());
      patElectron_gsfTrack_vtx.push_back(el->gsfTrack()->vx()); 
      patElectron_gsfTrack_vty.push_back(el->gsfTrack()->vy());
      patElectron_gsfTrack_vtz.push_back(el->gsfTrack()->vz());
      patElectron_gsfTrack_dxy_pv.push_back(el->gsfTrack()->dxy(firstGoodVertex->position()));
      patElectron_gsfTrack_dxy_bs.push_back(el->gsfTrack()->dxy(point));
      patElectron_gsfTrack_dz_pv.push_back(el->gsfTrack()->dz(firstGoodVertex->position()));
      patElectron_gsfTrack_dz_bs.push_back(el->gsfTrack()->dz(point));
      patElectron_dxyError.push_back(el->gsfTrack()->d0Error());
    
      // point of closest approach (PCA) to the beamspot and primary vertex
      TransientTrack elecTransTkPtr = theB->build(*(el->closestCtfTrackRef()));
//      TransientTrack elecTransTkPtr = theB->build(*(el->gsfTrack()));
      GlobalPoint patElectron_pca_bs = elecTransTkPtr.trajectoryStateClosestToPoint(thebs).position();
//      GlobalPoint patElectron_pca_bs(el->TrackPositionAtVtx().X(), el->TrackPositionAtVtx().Y(),el->TrackPositionAtVtx().Z());
      GlobalPoint patElectron_pca_pv = elecTransTkPtr.trajectoryStateClosestToPoint(thepv).position();
      patElectron_gsfTrack_PCAx_bs.push_back(patElectron_pca_bs.x());
      patElectron_gsfTrack_PCAy_bs.push_back(patElectron_pca_bs.y());
      patElectron_gsfTrack_PCAz_bs.push_back(patElectron_pca_bs.z());
      patElectron_gsfTrack_PCAx_pv.push_back(patElectron_pca_pv.x());
      patElectron_gsfTrack_PCAy_pv.push_back(patElectron_pca_pv.y());
      patElectron_gsfTrack_PCAz_pv.push_back(patElectron_pca_pv.z());

      // extract track fit errors
      const float elecMass = 0.000510998928;
      float elecSigma = elecMass*1e-6;
      float chi2 = 0.0;
      float ndf = 0.0;
      KinematicParticleFactoryFromTransientTrack pFactory;
      RefCountedKinematicParticle elecParticle = pFactory.particle(elecTransTkPtr, elecMass, chi2, ndf, elecSigma);
      patElectron_gsfTrackFitErrorMatrix_00.push_back(elecParticle->stateAtPoint(patElectron_pca_bs).kinematicParametersError().matrix()(0,0));
      patElectron_gsfTrackFitErrorMatrix_01.push_back(elecParticle->stateAtPoint(patElectron_pca_bs).kinematicParametersError().matrix()(0,1));
      patElectron_gsfTrackFitErrorMatrix_02.push_back(elecParticle->stateAtPoint(patElectron_pca_bs).kinematicParametersError().matrix()(0,2));
      patElectron_gsfTrackFitErrorMatrix_11.push_back(elecParticle->stateAtPoint(patElectron_pca_bs).kinematicParametersError().matrix()(1,1));
      patElectron_gsfTrackFitErrorMatrix_12.push_back(elecParticle->stateAtPoint(patElectron_pca_bs).kinematicParametersError().matrix()(1,2));
      patElectron_gsfTrackFitErrorMatrix_22.push_back(elecParticle->stateAtPoint(patElectron_pca_bs).kinematicParametersError().matrix()(2,2));

    }

    elcoun++;
  }
  
}

void ElectronPatSelector::SetBranches(){
  if(debug_) std::cout << "     ElectronPatSelector: Setting branches by calling AddBranch of baseTree." << std::endl;
  
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
  AddBranch(&passMVAwp1Id_                 ,"patElectron_passMV1wp1Id");
  AddBranch(&passMVAwp2Id_                 ,"patElectron_passMV2wp1Id");
  AddBranch(&isoChargedHadrons_            ,"patElectron_isoChargedHadrons");
  AddBranch(&isoNeutralHadrons_            ,"patElectron_isoNeutralHadrons");
  AddBranch(&isoPhotons_                   ,"patElectron_isoPhotons");
  AddBranch(&isoPU_                        ,"patElectron_isoPU");
  AddBranch(&expectedMissingInnerHits      ,"patElectron_expectedMissingInnerHits");
  AddBranch(&passConversionVeto_           ,"patElectron_passConversionVeto"); 

  if (!_super_TNT){
    AddBranch(&patElectron_gsfTrack_dxy_pv   ,"patElectron_gsfTrack_dxy_pv");
    AddBranch(&patElectron_gsfTrack_dxy_bs   ,"patElectron_gsfTrack_dxy_bs");
    AddBranch(&patElectron_dxyError          ,"patElectron_dxyError");
    AddBranch(&patElectron_gsfTrack_dz_pv    ,"patElectron_gsfTrack_dz_pv");
    AddBranch(&patElectron_gsfTrack_dz_bs    ,"patElectron_gsfTrack_dz_bs");
    AddBranch(&patElectron_gsfTrack_normChi2 ,"patElectron_gsfTrack_normChi2");
    AddBranch(&patElectron_gsfTrack_ndof     ,"patElectron_gsfTrack_ndof");
    AddBranch(&patElectron_gsfTrack_vtx      ,"patElectron_gsfTrack_vtx");
    AddBranch(&patElectron_gsfTrack_vty      ,"patElectron_gsfTrack_vty");
    AddBranch(&patElectron_gsfTrack_vtz      ,"patElectron_gsfTrack_vtz");
    AddBranch(&patElectron_gsfTrack_PCAx_bs  ,"patElectron_gsfTrack_PCAx_bs");
    AddBranch(&patElectron_gsfTrack_PCAy_bs  ,"patElectron_gsfTrack_PCAy_bs");
    AddBranch(&patElectron_gsfTrack_PCAz_bs  ,"patElectron_gsfTrack_PCAz_bs");
    AddBranch(&patElectron_gsfTrack_PCAx_pv  ,"patElectron_gsfTrack_PCAx_pv");
    AddBranch(&patElectron_gsfTrack_PCAy_pv  ,"patElectron_gsfTrack_PCAy_pv");
    AddBranch(&patElectron_gsfTrack_PCAz_pv  ,"patElectron_gsfTrack_PCAz_pv");
    AddBranch(&patElectron_gsfTrackFitErrorMatrix_00     ,"patElectron_gsfTrackFitErrorMatrix_00");
    AddBranch(&patElectron_gsfTrackFitErrorMatrix_01     ,"patElectron_gsfTrackFitErrorMatrix_01");
    AddBranch(&patElectron_gsfTrackFitErrorMatrix_02     ,"patElectron_gsfTrackFitErrorMatrix_02");
    AddBranch(&patElectron_gsfTrackFitErrorMatrix_11     ,"patElectron_gsfTrackFitErrorMatrix_11");
    AddBranch(&patElectron_gsfTrackFitErrorMatrix_12     ,"patElectron_gsfTrackFitErrorMatrix_12");
    AddBranch(&patElectron_gsfTrackFitErrorMatrix_22     ,"patElectron_gsfTrackFitErrorMatrix_22");
  }

  if(debug_) std::cout << "     ElectronPatSelector: Finished setting branches." << std::endl;
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
  passMVAwp1Id_.clear();
  passMVAwp2Id_.clear();
  patElectron_gsfTrack_dxy_pv.clear();
  patElectron_gsfTrack_dxy_bs.clear();
  patElectron_dxyError.clear();
  patElectron_gsfTrack_dz_pv.clear();
  patElectron_gsfTrack_dz_bs.clear();
  patElectron_gsfTrack_normChi2.clear();
  patElectron_gsfTrack_ndof.clear();
  patElectron_gsfTrack_vtx.clear(); 
  patElectron_gsfTrack_vty.clear();
  patElectron_gsfTrack_vtz.clear();
  expectedMissingInnerHits.clear();
  passConversionVeto_.clear();
  isoChargedHadrons_.clear();
  isoNeutralHadrons_.clear();
  isoPhotons_.clear();
  isoPU_.clear();
  patElectron_gsfTrack_PCAx_bs.clear();
  patElectron_gsfTrack_PCAy_bs.clear();
  patElectron_gsfTrack_PCAz_bs.clear();
  patElectron_gsfTrack_PCAx_pv.clear();
  patElectron_gsfTrack_PCAy_pv.clear();
  patElectron_gsfTrack_PCAz_pv.clear();
  patElectron_gsfTrackFitErrorMatrix_00.clear();
  patElectron_gsfTrackFitErrorMatrix_01.clear();
  patElectron_gsfTrackFitErrorMatrix_02.clear();
  patElectron_gsfTrackFitErrorMatrix_11.clear();
  patElectron_gsfTrackFitErrorMatrix_12.clear();
  patElectron_gsfTrackFitErrorMatrix_22.clear();

}

bool ElectronPatSelector::isGoodVertex(const reco::Vertex& vtxx) {
  if (vtxx.isFake()) return false;
  if (vtxx.ndof() < _patElectron_vtx_ndof_min) return false;
  if (vtxx.position().Rho() > _patElectron_vtx_rho_max) return false;
  if (fabs(vtxx.position().Z()) > _patElectron_vtx_position_z_max) return false;
  return true;
}

