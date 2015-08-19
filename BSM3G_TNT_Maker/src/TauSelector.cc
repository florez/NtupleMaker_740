// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#include "NtupleMaker/BSM3G_TNT_Maker/interface/TauSelector.h"

using namespace edm;
using namespace std;
using namespace reco;

TauSelector::TauSelector(std::string name, TTree* tree, bool debug, const pset& iConfig):
 baseTree(name,tree,debug){
  tauToken_       = iConfig.getParameter<edm::InputTag>("taus");
  _vertexInputTag = iConfig.getParameter<edm::InputTag>("vertices");
  _beamSpot                = iConfig.getParameter<edm::InputTag>("beamSpot");
  _Tau_pt_min     = iConfig.getParameter<double>("Tau_pt_min");
  _Tau_eta_max    = iConfig.getParameter<double>("Tau_eta_max");
  _Tau_vtx_ndof_min       = iConfig.getParameter<int>("Tau_vtx_ndof_min");
  _Tau_vtx_rho_max        = iConfig.getParameter<int>("Tau_vtx_rho_max");
  _Tau_vtx_position_z_max = iConfig.getParameter<double>("Tau_vtx_position_z_max");
  _super_TNT      = iConfig.getParameter<bool>("super_TNT");
  SetBranches();
}

TauSelector::~TauSelector(){
  delete tree_;
}

void TauSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();

  int taucoun = 0;

  // grab handle to the tau collection
  edm::Handle<pat::TauCollection> taus;
  iEvent.getByLabel(tauToken_, taus);

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

  // Get BeamSpot Information
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(_beamSpot, beamSpotHandle);
  if ( beamSpotHandle.isValid() ) { beamSpot = *beamSpotHandle; }
  else { edm::LogInfo("MyAnalyzer") << "No beam spot available from EventSetup \n"; }
  math::XYZPoint point(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());
  GlobalPoint thebs(beamSpot.x0(),beamSpot.y0(),beamSpot.z0());

  // get magnetic field and detector geometry information --> used to build transient tracks
//  edm::ESHandle<MagneticField> magfield;
//  iSetup.get<IdealMagneticFieldRecord>().get(magfield);
//  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
//  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
  
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  // loop over taus
  for (const pat::Tau &tau : *taus) {
    if (tau.pt() < _Tau_pt_min) continue;
    if(fabs(tau.eta()) > _Tau_eta_max) continue;
    if(!(tau.leadChargedHadrCand().isNonnull())) continue;
    if(!(tau.leadTrack().isNonnull())) continue;
    
    // fill root tree with "necessary" information:  kinematics, ID, isolation
    Tau_eta.push_back( tau.eta());
    Tau_phi.push_back( tau.phi());
    Tau_pt.push_back( tau.pt());
    Tau_energy.push_back( tau.energy());
    Tau_charge.push_back( tau.charge());
    Tau_nProngs.push_back(tau.signalChargedHadrCands().size());
    Tau_chargedIsoPtSum.push_back(tau.tauID("chargedIsoPtSum"));
    Tau_neutralIsoPtSum.push_back(tau.tauID("neutralIsoPtSum"));
    Tau_puCorrPtSum.push_back(tau.tauID("puCorrPtSum"));

    // tau id. discriminators
    Tau_decayModeFinding.push_back( tau.tauID("decayModeFinding"));
    Tau_decayModeFindingNewDMs.push_back( tau.tauID("decayModeFindingNewDMs"));
    Tau_byLooseIsolationMVA3newDMwLT.push_back(tau.tauID("byLooseIsolationMVA3newDMwLT"));
    Tau_byLooseIsolationMVA3newDMwoLT.push_back(tau.tauID("byLooseIsolationMVA3newDMwoLT"));
    Tau_byLooseIsolationMva3oldDMwLT.push_back(tau.tauID("byLooseIsolationMVA3oldDMwLT"));
    Tau_byLooseIsolationMVA3oldDMwoLT.push_back(tau.tauID("byLooseIsolationMVA3oldDMwoLT"));
    Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
    Tau_byMediumIsolationMVA3newDMwLT.push_back(tau.tauID("byMediumIsolationMVA3newDMwLT"));
    Tau_byMediumIsolationMVA3newDMwoLT.push_back(tau.tauID("byMediumIsolationMVA3newDMwoLT"));
    Tau_byMediumIsolationMva3oldDMwLT.push_back(tau.tauID("byMediumIsolationMVA3oldDMwLT"));
    Tau_byMediumIsolationMVA3oldDMwoLT.push_back(tau.tauID("byMediumIsolationMVA3oldDMwoLT"));
    Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
    Tau_byTightIsolationMVA3newDMwLT.push_back(tau.tauID("byTightIsolationMVA3newDMwLT"));
    Tau_byTightIsolationMVA3newDMwoLT.push_back(tau.tauID("byTightIsolationMVA3newDMwoLT"));
    Tau_byTightIsolationMva3oldDMwLT.push_back(tau.tauID("byTightIsolationMVA3oldDMwLT"));
    Tau_byTightIsolationMVA3oldDMwoLT.push_back(tau.tauID("byTightIsolationMVA3oldDMwoLT"));
    Tau_byTightCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));

    // discriminators against electrons/muons
    Tau_againstMuonLoose2.push_back( tau.tauID("againstMuonLoose2"));
    Tau_againstMuonLoose3.push_back( tau.tauID("againstMuonLoose3"));
    Tau_againstMuonTight2.push_back( tau.tauID("againstMuonTight2"));
    Tau_againstMuonTight3.push_back( tau.tauID("againstMuonTight3"));
    Tau_againstElectronMVALooseMVA5.push_back( tau.tauID("againstElectronLooseMVA5"));
    Tau_againstElectronMVAMediumMVA5.push_back( tau.tauID("againstElectronMediumMVA5"));
    Tau_againstElectronMVATightMVA5.push_back( tau.tauID("againstElectronTightMVA5"));

    // store additional information such as lifetime variables
    if(!_super_TNT){
    
      const reco::CandidatePtr hadTauLeadChargedCand = tau.leadChargedHadrCand();                                                                   
      Tau_leadChargedCandPt.push_back(hadTauLeadChargedCand.isNonnull() ? hadTauLeadChargedCand->pt() : -1.);      
      Tau_leadChargedCandEta.push_back(hadTauLeadChargedCand.isNonnull() ? hadTauLeadChargedCand->eta() : -999);
      Tau_leadChargedCandPhi.push_back(hadTauLeadChargedCand.isNonnull() ? hadTauLeadChargedCand->phi() : -999);
      Tau_leadChargedCandCharge.push_back(hadTauLeadChargedCand.isNonnull() ? hadTauLeadChargedCand->charge() : -2);   
      Tau_leadChargedCandChi2.push_back(tau.leadTrack()->chi2());
      Tau_leadChargedCandValidHits.push_back(tau.leadTrack()->numberOfValidHits());
      Tau_leadChargedCandDxy_pv.push_back(tau.leadTrack()->dxy(firstGoodVertex->position()));
      Tau_leadChargedCandDxy_bs.push_back(-1.*(tau.leadTrack()->dxy(point)));
      Tau_leadChargedCandDz_pv.push_back(tau.leadTrack()->dz(firstGoodVertex->position()));
      Tau_leadChargedCandDz_bs.push_back(tau.leadTrack()->dz(point));
      Tau_leadChargedCandDxyError.push_back(tau.leadTrack()->d0Error());
      Tau_leadChargedCandDzError.push_back(tau.leadTrack()->dzError());
      Tau_leadChargedCandNdof.push_back(tau.leadTrack()->ndof());
      Tau_leadChargedCandVtx.push_back(tau.leadTrack()->vx());
      Tau_leadChargedCandVty.push_back(tau.leadTrack()->vy());
      Tau_leadChargedCandVtz.push_back(tau.leadTrack()->vz());
      Tau_leadChargedCandTrack_pt.push_back(tau.leadTrack()->pt());
      Tau_leadChargedCandTrack_ptError.push_back(tau.leadTrack()->ptError());

      // default tau POG lifetime variables
      Tau_defaultDxy.push_back(tau.dxy());
      Tau_defaultDxyError.push_back(tau.dxy_error());
      Tau_defaultDxySig.push_back(tau.dxy_Sig());
      Tau_defaultFlightLengthX.push_back(tau.flightLength().x());
      Tau_defaultFlightLengthY.push_back(tau.flightLength().y());
      Tau_defaultFlightLengthZ.push_back(tau.flightLength().z());
      Tau_defaultFlightLengthSig.push_back(tau.flightLengthSig());
      Tau_default_PCAx_pv.push_back(tau.dxy_PCA().x());
      Tau_default_PCAy_pv.push_back(tau.dxy_PCA().y());
      Tau_default_PCAz_pv.push_back(tau.dxy_PCA().z());

      // tau lead track point of closest approach (PCA) to the beamspot and primary vertex
      TransientTrack tauTransTkPtr = theB->build(*(tau.leadTrack()));
//      TransientTrack tauTransTkPtr(*(tau.leadTrack()), magfield.product(),theTrackingGeometry);
      GlobalPoint tauLeadTrack_pca_bs = tauTransTkPtr.trajectoryStateClosestToPoint(thebs).position();
      GlobalPoint tauLeadTrack_pca_pv = tauTransTkPtr.trajectoryStateClosestToPoint(thepv).position();
      Tau_leadChargedCandTrack_PCAx_bs.push_back(tauLeadTrack_pca_bs.x());
      Tau_leadChargedCandTrack_PCAy_bs.push_back(tauLeadTrack_pca_bs.y());
      Tau_leadChargedCandTrack_PCAz_bs.push_back(tauLeadTrack_pca_bs.z());
      Tau_leadChargedCandTrack_PCAx_pv.push_back(tauLeadTrack_pca_pv.x());
      Tau_leadChargedCandTrack_PCAy_pv.push_back(tauLeadTrack_pca_pv.y());
      Tau_leadChargedCandTrack_PCAz_pv.push_back(tauLeadTrack_pca_pv.z());

      // extract track fit errors
      const float pionMass = 0.139570;
      float pionSigma = pionMass*1e-6;
      float chi2 = 0.0;
      float ndf = 0.0;
      KinematicParticleFactoryFromTransientTrack pFactory;
      RefCountedKinematicParticle tauParticle = pFactory.particle(tauTransTkPtr, pionMass, chi2, ndf, pionSigma);
      Tau_leadChargedCandTrackFitErrorMatrix_00.push_back(tauParticle->stateAtPoint(tauLeadTrack_pca_bs).kinematicParametersError().matrix()(0,0));
      Tau_leadChargedCandTrackFitErrorMatrix_01.push_back(tauParticle->stateAtPoint(tauLeadTrack_pca_bs).kinematicParametersError().matrix()(0,1));
      Tau_leadChargedCandTrackFitErrorMatrix_02.push_back(tauParticle->stateAtPoint(tauLeadTrack_pca_bs).kinematicParametersError().matrix()(0,2));
      Tau_leadChargedCandTrackFitErrorMatrix_11.push_back(tauParticle->stateAtPoint(tauLeadTrack_pca_bs).kinematicParametersError().matrix()(1,1));
      Tau_leadChargedCandTrackFitErrorMatrix_12.push_back(tauParticle->stateAtPoint(tauLeadTrack_pca_bs).kinematicParametersError().matrix()(1,2));
      Tau_leadChargedCandTrackFitErrorMatrix_22.push_back(tauParticle->stateAtPoint(tauLeadTrack_pca_bs).kinematicParametersError().matrix()(2,2));
    }
    taucoun++;
  }
  
}

void TauSelector::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;

  AddBranch(&Tau_eta,                         "Tau_eta");
  AddBranch(&Tau_phi,                         "Tau_phi");
  AddBranch(&Tau_pt,                          "Tau_pt");
  AddBranch(&Tau_energy,                      "Tau_energy");
  AddBranch(&Tau_charge,                      "Tau_charge");
  AddBranch(&Tau_nProngs,                     "Tau_nProngs");
  AddBranch(&Tau_chargedIsoPtSum,             "Tau_chargedIsoPtSum");
  AddBranch(&Tau_neutralIsoPtSum,             "Tau_neutralIsoPtSum");
  AddBranch(&Tau_puCorrPtSum,                 "Tau_puCorrPtSum");
  AddBranch(&Tau_decayModeFinding,            "Tau_decayModeFinding");
  AddBranch(&Tau_decayModeFindingNewDMs,      "Tau_decayModeFindingNewDMs");
  AddBranch(&Tau_byLooseIsolationMVA3newDMwLT,                  "Tau_byLooseIsolationMVA3newDMwLT");
  AddBranch(&Tau_byLooseIsolationMVA3newDMwoLT,                 "Tau_byLooseIsolationMVA3newDMwoLT");
  AddBranch(&Tau_byLooseIsolationMva3oldDMwLT,                  "Tau_byLooseIsolationMva3oldDMwLT");
  AddBranch(&Tau_byLooseIsolationMVA3oldDMwoLT,                 "Tau_byLooseIsolationMVA3oldDMwoLT");
  AddBranch(&Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits,    "Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits");
  AddBranch(&Tau_byMediumIsolationMVA3newDMwLT,                 "Tau_byMediumIsolationMVA3newDMwLT");
  AddBranch(&Tau_byMediumIsolationMVA3newDMwoLT,                "Tau_byMediumIsolationMVA3newDMwoLT");
  AddBranch(&Tau_byMediumIsolationMva3oldDMwLT,                 "Tau_byMediumIsolationMva3oldDMwLT");
  AddBranch(&Tau_byMediumIsolationMVA3oldDMwoLT,                "Tau_byMediumIsolationMVA3oldDMwoLT");
  AddBranch(&Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits,   "Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits");
  AddBranch(&Tau_byTightIsolationMVA3newDMwLT,                  "Tau_byTightIsolationMVA3newDMwLT");
  AddBranch(&Tau_byTightIsolationMVA3newDMwoLT,                 "Tau_byTightIsolationMVA3newDMwoLT");
  AddBranch(&Tau_byTightIsolationMva3oldDMwLT,                  "Tau_byTightIsolationMva3oldDMwLT");
  AddBranch(&Tau_byTightIsolationMVA3oldDMwoLT,                 "Tau_byTightIsolationMVA3oldDMwoLT");
  AddBranch(&Tau_byTightCombinedIsolationDeltaBetaCorr3Hits,    "Tau_byTightCombinedIsolationDeltaBetaCorr3Hits");
  AddBranch(&Tau_againstMuonLoose2,                             "Tau_againstMuonLoose2");
  AddBranch(&Tau_againstMuonLoose3,                             "Tau_againstMuonLoose3");
  AddBranch(&Tau_againstMuonTight2,                             "Tau_againstMuonTight2");
  AddBranch(&Tau_againstMuonTight3,                             "Tau_againstMuonTight3"); 
  AddBranch(&Tau_againstElectronMVALooseMVA5,                   "Tau_againstElectronMVALooseMVA5");
  AddBranch(&Tau_againstElectronMVAMediumMVA5,                  "Tau_againstElectronMVAMediumMVA5");
  AddBranch(&Tau_againstElectronMVATightMVA5,                   "Tau_againstElectronMVATightMVA5");
  if(!_super_TNT){ 
    AddBranch(&Tau_leadChargedCandPt,                           "Tau_leadChargedCandPt");
    AddBranch(&Tau_leadChargedCandEta,                          "Tau_leadChargedCandEta");
    AddBranch(&Tau_leadChargedCandPhi,                          "Tau_leadChargedCandPhi");
    AddBranch(&Tau_leadChargedCandCharge,                       "Tau_leadChargedCandCharge");
    AddBranch(&Tau_leadChargedCandChi2,                         "Tau_leadChargedCandChi2");
    AddBranch(&Tau_leadChargedCandValidHits,                    "Tau_leadChargedCandValidHits");
    AddBranch(&Tau_leadChargedCandDxy_pv,                       "Tau_leadChargedCandDxy_pv");
    AddBranch(&Tau_leadChargedCandDxy_bs,                       "Tau_leadChargedCandDxy_bs");
    AddBranch(&Tau_leadChargedCandDz_pv,                        "Tau_leadChargedCandDz_pv");
    AddBranch(&Tau_leadChargedCandDz_bs,                        "Tau_leadChargedCandDz_bs");
    AddBranch(&Tau_leadChargedCandDzError,                      "Tau_leadChargedCandDzError");
    AddBranch(&Tau_leadChargedCandDxyError,                     "Tau_leadChargedCandDxyError");
    AddBranch(&Tau_leadChargedCandNdof,                         "Tau_leadChargedCandNdof");
    AddBranch(&Tau_leadChargedCandVtx,                          "Tau_leadChargedCandVtx");
    AddBranch(&Tau_leadChargedCandVty,                          "Tau_leadChargedCandVty");
    AddBranch(&Tau_leadChargedCandVtz,                          "Tau_leadChargedCandVtz");
    AddBranch(&Tau_leadChargedCandTrack_pt,                     "Tau_leadChargedCandTrack_pt");
    AddBranch(&Tau_leadChargedCandTrack_ptError,                "Tau_leadChargedCandTrack_ptError");
    AddBranch(&Tau_leadChargedCandTrack_PCAx_bs              ,"Tau_leadChargedCandTrack_PCAx_bs");
    AddBranch(&Tau_leadChargedCandTrack_PCAy_bs              ,"Tau_leadChargedCandTrack_PCAy_bs");
    AddBranch(&Tau_leadChargedCandTrack_PCAz_bs              ,"Tau_leadChargedCandTrack_PCAz_bs");
    AddBranch(&Tau_leadChargedCandTrack_PCAx_pv              ,"Tau_leadChargedCandTrack_PCAx_pv");
    AddBranch(&Tau_leadChargedCandTrack_PCAy_pv              ,"Tau_leadChargedCandTrack_PCAy_pv");
    AddBranch(&Tau_leadChargedCandTrack_PCAz_pv              ,"Tau_leadChargedCandTrack_PCAz_pv");
    AddBranch(&Tau_leadChargedCandTrackFitErrorMatrix_00     ,"Tau_leadChargedCandTrackFitErrorMatrix_00");
    AddBranch(&Tau_leadChargedCandTrackFitErrorMatrix_01     ,"Tau_leadChargedCandTrackFitErrorMatrix_01");
    AddBranch(&Tau_leadChargedCandTrackFitErrorMatrix_02     ,"Tau_leadChargedCandTrackFitErrorMatrix_02");
    AddBranch(&Tau_leadChargedCandTrackFitErrorMatrix_11     ,"Tau_leadChargedCandTrackFitErrorMatrix_11");
    AddBranch(&Tau_leadChargedCandTrackFitErrorMatrix_12     ,"Tau_leadChargedCandTrackFitErrorMatrix_12");
    AddBranch(&Tau_leadChargedCandTrackFitErrorMatrix_22     ,"Tau_leadChargedCandTrackFitErrorMatrix_22");
    AddBranch(&Tau_defaultDxy                                ,"Tau_defaultDxy");
    AddBranch(&Tau_defaultDxyError                           ,"Tau_defaultDxyError");
    AddBranch(&Tau_defaultDxySig                             ,"Tau_defaultDxySig");
    AddBranch(&Tau_defaultFlightLengthX                      ,"Tau_defaultFlightLengthX");
    AddBranch(&Tau_defaultFlightLengthY                      ,"Tau_defaultFlightLengthY");
    AddBranch(&Tau_defaultFlightLengthZ                      ,"Tau_defaultFlightLengthZ");
    AddBranch(&Tau_defaultFlightLengthSig                    ,"Tau_defaultFlightLengthSig");
    AddBranch(&Tau_default_PCAx_pv                           ,"Tau_default_PCAx_pv");
    AddBranch(&Tau_default_PCAy_pv                           ,"Tau_default_PCAy_pv");
    AddBranch(&Tau_default_PCAz_pv                           ,"Tau_default_PCAz_pv");
  }
}

void TauSelector::Clear(){
  
  Tau_eta.clear(); 
  Tau_phi.clear(); 
  Tau_pt.clear(); 
  Tau_energy.clear(); 
  Tau_charge.clear(); 
  Tau_nProngs.clear();
  Tau_chargedIsoPtSum.clear(); 
  Tau_neutralIsoPtSum.clear(); 
  Tau_puCorrPtSum.clear(); 
  Tau_decayModeFinding.clear(); 
  Tau_decayModeFindingNewDMs.clear(); 
  Tau_byLooseIsolationMVA3newDMwLT.clear();
  Tau_byLooseIsolationMVA3newDMwoLT.clear();
  Tau_byLooseIsolationMva3oldDMwLT.clear();
  Tau_byLooseIsolationMVA3oldDMwoLT.clear();
  Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.clear(); 
  Tau_byMediumIsolationMVA3newDMwLT.clear();
  Tau_byMediumIsolationMVA3newDMwoLT.clear();
  Tau_byMediumIsolationMva3oldDMwLT.clear();
  Tau_byMediumIsolationMVA3oldDMwoLT.clear();
  Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.clear(); 
  Tau_byTightIsolationMVA3newDMwLT.clear();
  Tau_byTightIsolationMVA3newDMwoLT.clear();
  Tau_byTightIsolationMva3oldDMwLT.clear();
  Tau_byTightIsolationMVA3oldDMwoLT.clear();
  Tau_byTightCombinedIsolationDeltaBetaCorr3Hits.clear(); 
  Tau_againstMuonLoose2.clear();
  Tau_againstMuonLoose3.clear();
  Tau_againstMuonTight2.clear();
  Tau_againstMuonTight3.clear();
  Tau_againstElectronMVALooseMVA5.clear();
  Tau_againstElectronMVAMediumMVA5.clear();
  Tau_againstElectronMVATightMVA5.clear();
  Tau_leadChargedCandPt.clear();
  Tau_leadChargedCandCharge.clear();
  Tau_leadChargedCandEta.clear();
  Tau_leadChargedCandPhi.clear();
  Tau_leadChargedCandChi2.clear();
  Tau_leadChargedCandValidHits.clear();
  Tau_leadChargedCandDxy_pv.clear();
  Tau_leadChargedCandDxy_bs.clear();
  Tau_leadChargedCandDz_pv.clear();
  Tau_leadChargedCandDz_bs.clear();
  Tau_leadChargedCandDzError.clear();
  Tau_leadChargedCandDxyError.clear();
  Tau_leadChargedCandNdof.clear();
  Tau_leadChargedCandVtx.clear();
  Tau_leadChargedCandVty.clear();
  Tau_leadChargedCandVtz.clear();
  Tau_leadChargedCandTrack_pt.clear();
  Tau_leadChargedCandTrack_ptError.clear();
  Tau_leadChargedCandTrack_PCAx_bs.clear();
  Tau_leadChargedCandTrack_PCAy_bs.clear();
  Tau_leadChargedCandTrack_PCAz_bs.clear();
  Tau_leadChargedCandTrack_PCAx_pv.clear();
  Tau_leadChargedCandTrack_PCAy_pv.clear();
  Tau_leadChargedCandTrack_PCAz_pv.clear();
  Tau_leadChargedCandTrackFitErrorMatrix_00.clear();
  Tau_leadChargedCandTrackFitErrorMatrix_01.clear();
  Tau_leadChargedCandTrackFitErrorMatrix_02.clear();
  Tau_leadChargedCandTrackFitErrorMatrix_11.clear();
  Tau_leadChargedCandTrackFitErrorMatrix_12.clear();
  Tau_leadChargedCandTrackFitErrorMatrix_22.clear();
  Tau_defaultDxy.clear();
  Tau_defaultDxyError.clear();
  Tau_defaultDxySig.clear();
  Tau_defaultFlightLengthX.clear();
  Tau_defaultFlightLengthY.clear();
  Tau_defaultFlightLengthZ.clear();
  Tau_defaultFlightLengthSig.clear();
  Tau_default_PCAx_pv.clear();
  Tau_default_PCAy_pv.clear();
  Tau_default_PCAz_pv.clear();
}

bool TauSelector::isGoodVertex(const reco::Vertex& vtxxx) {
  if (vtxxx.isFake()) return false;
  if (vtxxx.ndof() < _Tau_vtx_ndof_min) return false;
  if (vtxxx.position().Rho() > _Tau_vtx_rho_max) return false;
  if (fabs(vtxxx.position().Z()) > _Tau_vtx_position_z_max) return false;
  return true;
}

