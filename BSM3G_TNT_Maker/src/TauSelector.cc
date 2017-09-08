// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#include "NtupleMaker/BSM3G_TNT_Maker/interface/TauSelector.h"

using namespace edm;
using namespace std;
using namespace reco;

TauSelector::TauSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector&& iCC):baseTree(name,tree,debug) {
  if(debug) std::cout << "BSM3G TNT Maker: In the TauSelector Constructor --> getting parameters & calling SetBranches()." << std::endl;
  tauToken_                     = iCC.consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("taus"));
  packedPFCandidateToken_       = iCC.consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCandidates"));
  _vertexInputTag               = iCC.consumes<reco::VertexCollection >(iConfig.getParameter<edm::InputTag>("vertices"));
  _beamSpot                     = iCC.consumes<reco::BeamSpot >(iConfig.getParameter<edm::InputTag>("beamSpot"));
  _Tau_pt_min     		= iConfig.getParameter<double>("Tau_pt_min");
  _Tau_eta_max    		= iConfig.getParameter<double>("Tau_eta_max");
  _Tau_vtx_ndof_min       	= iConfig.getParameter<int>("Tau_vtx_ndof_min");
  _Tau_vtx_rho_max        	= iConfig.getParameter<int>("Tau_vtx_rho_max");
  _Tau_vtx_position_z_max 	= iConfig.getParameter<double>("Tau_vtx_position_z_max");
  _super_TNT      		= iConfig.getParameter<bool>("super_TNT");
  SetBranches();
}

TauSelector::~TauSelector(){
  delete tree_;
}

void TauSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();

  if(debug_) std::cout << "     TauSelector: Cleared the vectors, grabbing the tau/PV/BeamSpot collection handles." << std::endl;

  int taucoun = 0;

  // grab handle to the tau collection
  edm::Handle<edm::View<pat::Tau> > taus;
  iEvent.getByToken(tauToken_, taus);

  if(debug_) std::cout << "     TauSelector: Looping over PVs to extract the position of the best PV." << std::endl;

  // grab handle to the packed pf candidate collection
  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(packedPFCandidateToken_, pfs);

  // grab handle to the vertex collection
  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByToken(_vertexInputTag, vtx);
  reco::VertexCollection::const_iterator firstGoodVertex = vtx->end();
  for (reco::VertexCollection::const_iterator it = vtx->begin(); it != vtx->end(); it++) {
    if (isGoodVertex(*it)) {
      firstGoodVertex = it;
      break;
    }
  }
  if (firstGoodVertex == vtx->end()) return; // skip event if there are no good PVs
  GlobalPoint thepv(firstGoodVertex->position().x(),firstGoodVertex->position().y(),firstGoodVertex->position().z());

  if(debug_) std::cout << "     TauSelector: Extracting the beamspot position." << std::endl;

  // Get BeamSpot Information
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(_beamSpot, beamSpotHandle);
  if ( beamSpotHandle.isValid() ) { beamSpot = *beamSpotHandle; }
  else { edm::LogInfo("MyAnalyzer") << "No beam spot available from EventSetup \n"; }
  math::XYZPoint point(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());
  GlobalPoint thebs(beamSpot.x0(),beamSpot.y0(),beamSpot.z0());

  // get magnetic field and detector geometry information --> used to build transient tracks
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  if(debug_) std::cout << "     TauSelector: Looping over taus." << std::endl;

  // loop over taus
  int nTau=0;
  for(edm::View<pat::Tau>::const_iterator tau = taus->begin(); tau != taus->end(); tau++) {
    nTau++;
    if(debug_) std::cout << "          Tau info: Tau #" << nTau << ": pt = " << tau->pt() << ", eta = " << tau->eta() << ", phi = " << tau->phi() << ", isLeadChargedHadrCandNonNull = " << tau->leadChargedHadrCand().isNonnull() <<  std::endl;

    if (tau->pt() < _Tau_pt_min) continue;
    if(fabs(tau->eta()) > _Tau_eta_max) continue;
    if(!(tau->leadChargedHadrCand().isNonnull())) continue;

    if(debug_) std::cout << "          Tau info: Tau #" << nTau << " passed kinematic cuts & lead charged hadron IsNonNull cut" << std::endl;

    // get the embedded packed candidate of the leading PF charged hadron
    const reco::CandidatePtr hadTauLeadChargedCand = tau->leadChargedHadrCand();

    if(debug_) std::cout << "          Tau info: Tau #" << nTau << " --> Tau lead charged hadron: pt = " << hadTauLeadChargedCand->pt() << ", eta = " << hadTauLeadChargedCand->eta() << ", phi = " << hadTauLeadChargedCand->phi() << std::endl;

    // loop over packed PF candidates and find the one which matches the embedded packed candidate within the pat::Tau
    const reco::Track *leadTrack = 0;
    bool isBestTrackNonNull = false;
    bool leadPackedCandidateExists = false;
    for(unsigned int iPF = 0, nPF = pfs->size(); iPF < nPF; ++iPF) {
      const pat::PackedCandidate &pfCandidate = (*pfs)[iPF];
      if( (hadTauLeadChargedCand->pt() == pfCandidate.pt()) && (hadTauLeadChargedCand->eta() == pfCandidate.eta()) &&
          (hadTauLeadChargedCand->phi() == pfCandidate.phi()) ) { // the packed PF candidate and embedded lead candidate within the pat::Tau should be the same
        if(debug_) std::cout << "               PF Candidate info: PF Candidate #" << iPF << ": pt = " << pfCandidate.pt() << ", eta = " << pfCandidate.eta() << ", phi = " << pfCandidate.phi();
        leadPackedCandidateExists = true; // if there is a match between the packed PF candidate and embedded lead candidate within the pat::Tau
        if(pfCandidate.bestTrack() != NULL) {isBestTrackNonNull = true; leadTrack = pfCandidate.bestTrack();} // grab the associated CTF track (if it exists)
      }
    }
    if(debug_) std::cout << "          Tau info: Tau #" << nTau << ", leadPackedCandidateExists = " << leadPackedCandidateExists << ", isBestTrackNonNull = " << isBestTrackNonNull << std::endl;
    if(!(leadPackedCandidateExists)) continue; // throw away the tau if there was no matching packed PF candidate to the embedded lead candidate within the pat::Tau
    if(!(isBestTrackNonNull)) continue; // throw away the tau if it's lead charged hadron has no associated CTF track
    if(debug_) std::cout << "               Track Candidate info: Track Candidate pt = " << leadTrack->pt() << ", eta = " << leadTrack->eta() << ", phi = " << leadTrack->phi() << std::endl;
    if(debug_) std::cout << "          Tau info: Tau #" << nTau << " passed lead track IsNonNull cut" << std::endl;

    // fill root tree with "necessary" information:  kinematics, ID, isolation
    Tau_eta.push_back(tau->eta());
    Tau_phi.push_back(tau->phi());
    Tau_pt.push_back(tau->pt());
    Tau_energy.push_back(tau->energy());
    Tau_charge.push_back(tau->charge());
    Tau_decayMode.push_back(tau->decayMode());
    Tau_nProngs.push_back(tau->signalChargedHadrCands().size());
    //just to keep this for the next change in names
    //cout<<"-----------"<<endl;
    //for(size_t i=0;i< tau->tauIDs().size();i++){
      //cout<<tau->tauIDs()[i].first<<endl;
    //}
    
    Tau_chargedIsoPtSum.push_back(tau->tauID("chargedIsoPtSum"));
    Tau_neutralIsoPtSum.push_back(tau->tauID("neutralIsoPtSum"));
    Tau_puCorrPtSum.push_back(tau->tauID("puCorrPtSum"));
    Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits.push_back(tau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));

    // tau id. discriminators
    Tau_decayModeFinding.push_back(tau->tauID("decayModeFinding"));
    Tau_decayModeFindingNewDMs.push_back(tau->tauID("decayModeFindingNewDMs"));

    //CombinedIsolationDeltaBetaCorr3Hits
    Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.push_back(tau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
    Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(tau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
    Tau_byTightCombinedIsolationDeltaBetaCorr3Hits.push_back(tau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));

    //CombinedIsolationDeltaBetaCorr3Hits DR = 0.3
    //Tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03.push_back(tau->tauID("byLooseCombinedIsolationDeltaBetaCorr3HitsdR03"));
    //Tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03.push_back(tau->tauID("byMediumCombinedIsolationDeltaBetaCorr3HitsdR03"));
    //Tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03.push_back(tau->tauID("byTightCombinedIsolationDeltaBetaCorr3HitsdR03"));

    // MVA ID using iso with old Decay Mode reconstruction:
    Tau_byVLooseIsolationMVArun2v1DBoldDMwLT.push_back(tau->tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"));
    Tau_byLooseIsolationMVArun2v1DBoldDMwLT.push_back(tau->tauID("byLooseIsolationMVArun2v1DBoldDMwLT"));
    Tau_byMediumIsolationMVArun2v1DBoldDMwLT.push_back(tau->tauID("byMediumIsolationMVArun2v1DBoldDMwLT"));
    Tau_byTightIsolationMVArun2v1DBoldDMwLT.push_back(tau->tauID("byTightIsolationMVArun2v1DBoldDMwLT"));
    Tau_byVTightIsolationMVArun2v1DBoldDMwLT.push_back(tau->tauID("byVTightIsolationMVArun2v1DBoldDMwLT"));

    // MVA ID using Iso with new Decay Mode Reconstruction:
    Tau_byVLooseIsolationMVArun2v1DBnewDMwLT.push_back(tau->tauID("byVLooseIsolationMVArun2v1DBnewDMwLT"));
    Tau_byLooseIsolationMVArun2v1DBnewDMwLT.push_back(tau->tauID("byLooseIsolationMVArun2v1DBnewDMwLT"));
    Tau_byTightIsolationMVArun2v1DBnewDMwLT.push_back(tau->tauID("byTightIsolationMVArun2v1DBnewDMwLT"));
    Tau_byVTightIsolationMVArun2v1DBnewDMwLT.push_back(tau->tauID("byVTightIsolationMVArun2v1DBnewDMwLT"));
    Tau_byMediumIsolationMVArun2v1DBnewDMwLT.push_back(tau->tauID("byMediumIsolationMVArun2v1DBnewDMwLT"));

    // Iso MVA with Delta Beta Corrections and DR = 0.3
    Tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau->tauID("byVLooseIsolationMVArun2v1DBdR03oldDMwLT"));
    Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau->tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT"));
    Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau->tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT"));
    Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau->tauID("byTightIsolationMVArun2v1DBdR03oldDMwLT"));
    Tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau->tauID("byVTightIsolationMVArun2v1DBdR03oldDMwLT"));
    
    //get the raw mva output:
    Tau_byIsolationMVArun2v1DBdR03oldDMwLTraw.push_back(tau->tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw"));
    Tau_byIsolationMVArun2v1DBnewDMwLTraw.push_back(tau->tauID("byIsolationMVArun2v1DBnewDMwLTraw"));
    Tau_byIsolationMVArun2v1DBoldDMwLTraw.push_back(tau->tauID("byIsolationMVArun2v1DBoldDMwLTraw"));
    Tau_byIsolationMVArun2v1PWdR03oldDMwLTraw.push_back(tau->tauID("byIsolationMVArun2v1PWdR03oldDMwLTraw"));
    
    

    // discriminators against electrons/muons
    Tau_againstMuonLoose3.push_back(tau->tauID("againstMuonLoose3"));
    Tau_againstMuonTight3.push_back(tau->tauID("againstMuonTight3"));
    Tau_againstElectronMVAVLooseMVA6.push_back(tau->tauID("againstElectronVLooseMVA6"));
    Tau_againstElectronMVALooseMVA6.push_back(tau->tauID("againstElectronLooseMVA6"));
    Tau_againstElectronMVAMediumMVA6.push_back(tau->tauID("againstElectronMediumMVA6"));
    Tau_againstElectronMVATightMVA6.push_back(tau->tauID("againstElectronTightMVA6"));

    // store additional information such as lifetime variables
    if(!_super_TNT){

      Tau_leadChargedCandPt.push_back(hadTauLeadChargedCand.isNonnull() ? hadTauLeadChargedCand->pt() : -1.);
      Tau_leadChargedCandEta.push_back(hadTauLeadChargedCand.isNonnull() ? hadTauLeadChargedCand->eta() : -999);
      Tau_leadChargedCandPhi.push_back(hadTauLeadChargedCand.isNonnull() ? hadTauLeadChargedCand->phi() : -999);
      Tau_leadChargedCandCharge.push_back(hadTauLeadChargedCand.isNonnull() ? hadTauLeadChargedCand->charge() : -2);
      Tau_leadChargedCandChi2.push_back(leadTrack->chi2());
      Tau_leadChargedCandValidHits.push_back(leadTrack->numberOfValidHits());
      Tau_leadChargedCandDxy_pv.push_back(leadTrack->dxy(firstGoodVertex->position()));
      Tau_leadChargedCandDxy_bs.push_back(-1.*(leadTrack->dxy(point)));
      Tau_leadChargedCandDz_pv.push_back(leadTrack->dz(firstGoodVertex->position()));
      Tau_leadChargedCandDz_bs.push_back(leadTrack->dz(point));
      Tau_leadChargedCandDxyError.push_back(leadTrack->d0Error());
      Tau_leadChargedCandDzError.push_back(leadTrack->dzError());
      Tau_leadChargedCandNdof.push_back(leadTrack->ndof());
      Tau_leadChargedCandVtx.push_back(leadTrack->vx());
      Tau_leadChargedCandVty.push_back(leadTrack->vy());
      Tau_leadChargedCandVtz.push_back(leadTrack->vz());
      Tau_leadChargedCandTrack_pt.push_back(leadTrack->pt());
      Tau_leadChargedCandTrack_ptError.push_back(leadTrack->ptError());

      // default tau POG lifetime variables
      Tau_defaultDxy.push_back(tau->dxy());
      Tau_defaultDxyError.push_back(tau->dxy_error());
      Tau_defaultDxySig.push_back(tau->dxy_Sig());
      Tau_defaultFlightLengthX.push_back(tau->flightLength().x());
      Tau_defaultFlightLengthY.push_back(tau->flightLength().y());
      Tau_defaultFlightLengthZ.push_back(tau->flightLength().z());
      Tau_defaultFlightLengthSig.push_back(tau->flightLengthSig());
      Tau_default_PCAx_pv.push_back(tau->dxy_PCA().x());
      Tau_default_PCAy_pv.push_back(tau->dxy_PCA().y());
      Tau_default_PCAz_pv.push_back(tau->dxy_PCA().z());

      // tau lead track point of closest approach (PCA) to the beamspot and primary vertex
      TransientTrack tauTransTkPtr = theB->build(leadTrack);
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

      Tau_signal_cand_chargedHadr_pt.push_back(vector< float >());
      Tau_signal_cand_chargedHadr_eta.push_back(vector< float >());
      Tau_signal_cand_chargedHadr_phi.push_back(vector< float >());
      Tau_signal_cand_chargedHadr_energy.push_back(vector< float >());
      for(unsigned int i=1; i<tau->signalChargedHadrCands().size(); i++){
        Tau_signal_cand_chargedHadr_pt.back().push_back( tau->signalChargedHadrCands()[i]->pt());
        Tau_signal_cand_chargedHadr_eta.back().push_back( tau->signalChargedHadrCands()[i]->eta());
        Tau_signal_cand_chargedHadr_phi.back().push_back( tau->signalChargedHadrCands()[i]->phi());
        Tau_signal_cand_chargedHadr_energy.back().push_back( tau->signalChargedHadrCands()[i]->energy());
      }

      Tau_signal_cand_neutHadr_pt.push_back(vector< float >());
      Tau_signal_cand_neutHadr_eta.push_back(vector< float >());
      Tau_signal_cand_neutHadr_phi.push_back(vector< float >());
      Tau_signal_cand_neutHadr_energy.push_back(vector< float >());
      for(unsigned int i=1; i<tau->signalNeutrHadrCands().size(); i++){
        Tau_signal_cand_neutHadr_pt.back().push_back( tau->signalNeutrHadrCands()[i]->pt());
        Tau_signal_cand_neutHadr_eta.back().push_back( tau->signalNeutrHadrCands()[i]->eta());
        Tau_signal_cand_neutHadr_phi.back().push_back( tau->signalNeutrHadrCands()[i]->phi());
        Tau_signal_cand_neutHadr_energy.back().push_back( tau->signalNeutrHadrCands()[i]->energy());
      }

      Tau_signal_cand_neutGam_pt.push_back(vector< float >());
      Tau_signal_cand_neutGam_eta.push_back(vector< float >());
      Tau_signal_cand_neutGam_phi.push_back(vector< float >());
      Tau_signal_cand_neutGam_energy.push_back(vector< float >());
      for(unsigned int i=1; i<tau->signalGammaCands().size(); i++){
        Tau_signal_cand_neutGam_pt.back().push_back(        tau->signalGammaCands()[i]->pt());
        Tau_signal_cand_neutGam_eta.back().push_back(       tau->signalGammaCands()[i]->eta());
        Tau_signal_cand_neutGam_phi.back().push_back(       tau->signalGammaCands()[i]->phi());
        Tau_signal_cand_neutGam_energy.back().push_back(    tau->signalGammaCands()[i]->energy());
      }
    }
    taucoun++;
  }

}

void TauSelector::SetBranches(){
  if(debug_) std::cout << "     TauSelector: Setting branches by calling AddBranch of baseTree." << std::endl;

  AddBranch(&Tau_eta,                         "Tau_eta");
  AddBranch(&Tau_phi,                         "Tau_phi");
  AddBranch(&Tau_pt,                          "Tau_pt");
  AddBranch(&Tau_energy,                      "Tau_energy");
  AddBranch(&Tau_charge,                      "Tau_charge");
  AddBranch(&Tau_decayMode,                   "Tau_decayMode");
  AddBranch(&Tau_nProngs,                     "Tau_nProngs");
  AddBranch(&Tau_chargedIsoPtSum,             "Tau_chargedIsoPtSum");
  AddBranch(&Tau_neutralIsoPtSum,             "Tau_neutralIsoPtSum");
  AddBranch(&Tau_puCorrPtSum,                 "Tau_puCorrPtSum");
  AddBranch(&Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits,                 "Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits");
  AddBranch(&Tau_decayModeFinding,            "Tau_decayModeFinding");
  AddBranch(&Tau_decayModeFindingNewDMs,      "Tau_decayModeFindingNewDMs");
  AddBranch(&Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits,    "Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits");
  AddBranch(&Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits,   "Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits");
  AddBranch(&Tau_byTightCombinedIsolationDeltaBetaCorr3Hits,    "Tau_byTightCombinedIsolationDeltaBetaCorr3Hits");
//  AddBranch(&Tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03,"Tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03");
//  AddBranch(&Tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03,"Tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03");
//  AddBranch(&Tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03, "Tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03");
  //Iso with Old Decay Mode reconstruction:
  AddBranch(&Tau_byVLooseIsolationMVArun2v1DBoldDMwLT,          "Tau_byVLooseIsolationMVArun2v1DBoldDMwLT");
  AddBranch(&Tau_byLooseIsolationMVArun2v1DBoldDMwLT,           "Tau_byLooseIsolationMVArun2v1DBoldDMwLT");
  AddBranch(&Tau_byMediumIsolationMVArun2v1DBoldDMwLT,          "Tau_byMediumIsolationMVArun2v1DBoldDMwLT");
  AddBranch(&Tau_byTightIsolationMVArun2v1DBoldDMwLT,           "Tau_byTightIsolationMVArun2v1DBoldDMwLT");
  AddBranch(&Tau_byVTightIsolationMVArun2v1DBoldDMwLT,          "Tau_byVTightIsolationMVArun2v1DBoldDMwLT");

  //Iso with new Decay Mode reconstruction:
  AddBranch(&Tau_byVLooseIsolationMVArun2v1DBnewDMwLT,          "Tau_byVLooseIsolationMVArun2v1DBnewDMwLT");
  AddBranch(&Tau_byLooseIsolationMVArun2v1DBnewDMwLT,           "Tau_byLooseIsolationMVArun2v1DBnewDMwLT");
  AddBranch(&Tau_byMediumIsolationMVArun2v1DBnewDMwLT,          "Tau_byMediumIsolationMVArun2v1DBnewDMwLT");
  AddBranch(&Tau_byTightIsolationMVArun2v1DBnewDMwLT,           "Tau_byTightIsolationMVArun2v1DBnewDMwLT");
  AddBranch(&Tau_byVTightIsolationMVArun2v1DBnewDMwLT,          "Tau_byVTightIsolationMVArun2v1DBnewDMwLT");

  // Iso MVA with Delta Beta correction DR = 0.3
  AddBranch(&Tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT,      "Tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT");
  AddBranch(&Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT,       "Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT");
  AddBranch(&Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT,      "Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT");
  AddBranch(&Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT,       "Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT");
  AddBranch(&Tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT,      "Tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT");
  
  
  AddBranch(&Tau_byIsolationMVArun2v1DBdR03oldDMwLTraw,         "Tau_byIsolationMVArun2v1DBdR03oldDMwLTraw");
  AddBranch(&Tau_byIsolationMVArun2v1DBnewDMwLTraw,             "Tau_byIsolationMVArun2v1DBnewDMwLTraw");
  AddBranch(&Tau_byIsolationMVArun2v1DBoldDMwLTraw,             "Tau_byIsolationMVArun2v1DBoldDMwLTraw");
  AddBranch(&Tau_byIsolationMVArun2v1PWdR03oldDMwLTraw,         "Tau_byIsolationMVArun2v1PWdR03oldDMwLTraw");
  

  AddBranch(&Tau_againstMuonLoose3,                             "Tau_againstMuonLoose3");
  AddBranch(&Tau_againstMuonTight3,                             "Tau_againstMuonTight3");
  AddBranch(&Tau_againstElectronMVAVLooseMVA6,                  "Tau_againstElectronMVAVLooseMVA6");
  AddBranch(&Tau_againstElectronMVALooseMVA6,                   "Tau_againstElectronMVALooseMVA6");
  AddBranch(&Tau_againstElectronMVAMediumMVA6,                  "Tau_againstElectronMVAMediumMVA6");
  AddBranch(&Tau_againstElectronMVATightMVA6,                   "Tau_againstElectronMVATightMVA6");
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
    AddBranch(&Tau_signal_cand_chargedHadr_pt,"Tau_signal_cand_chargedHadr_pt");
    AddBranch(&Tau_signal_cand_chargedHadr_eta,"Tau_signal_cand_chargedHadr_eta");
    AddBranch(&Tau_signal_cand_chargedHadr_phi,"Tau_signal_cand_chargedHadr_phi");
    AddBranch(&Tau_signal_cand_chargedHadr_energy,"Tau_signal_cand_chargedHadr_energy");
    AddBranch(&Tau_signal_cand_neutHadr_pt,     "Tau_signal_cand_neutHadr_pt");
    AddBranch(&Tau_signal_cand_neutHadr_eta,    "Tau_signal_cand_neutHadr_eta");
    AddBranch(&Tau_signal_cand_neutHadr_phi,    "Tau_signal_cand_neutHadr_phi");
    AddBranch(&Tau_signal_cand_neutHadr_energy, "Tau_signal_cand_neutHadr_energy");

  }

  if(debug_) std::cout << "     TauSelector: Finished setting branches." << std::endl;
}

void TauSelector::Clear(){

  Tau_eta.clear();
  Tau_phi.clear();
  Tau_pt.clear();
  Tau_energy.clear();
  Tau_charge.clear();
  Tau_nProngs.clear();
  Tau_decayMode.clear();
  Tau_chargedIsoPtSum.clear();
  Tau_neutralIsoPtSum.clear();
  Tau_puCorrPtSum.clear();
  Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits.clear();
  Tau_decayModeFinding.clear();
  Tau_decayModeFindingNewDMs.clear();
  Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.clear();
  Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.clear();
  Tau_byTightCombinedIsolationDeltaBetaCorr3Hits.clear();
//  Tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03.clear();
//  Tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03.clear();
//  Tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03.clear();
  // iso with Old Decay Mode Reco
  Tau_byVLooseIsolationMVArun2v1DBoldDMwLT.clear();
  Tau_byLooseIsolationMVArun2v1DBoldDMwLT.clear();
  Tau_byMediumIsolationMVArun2v1DBoldDMwLT.clear();
  Tau_byTightIsolationMVArun2v1DBoldDMwLT.clear();
  Tau_byVTightIsolationMVArun2v1DBoldDMwLT.clear();

  //iso with new Decay Model Reco
  Tau_byVLooseIsolationMVArun2v1DBnewDMwLT.clear();
  Tau_byLooseIsolationMVArun2v1DBnewDMwLT.clear();
  Tau_byMediumIsolationMVArun2v1DBnewDMwLT.clear();
  Tau_byTightIsolationMVArun2v1DBnewDMwLT.clear();
  Tau_byVTightIsolationMVArun2v1DBnewDMwLT.clear();

  // iso with Delta Beta Corrections and DR = 0.3
  Tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT.clear();
  Tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT.clear();
  
  Tau_byIsolationMVArun2v1DBdR03oldDMwLTraw.clear();   
  Tau_byIsolationMVArun2v1DBnewDMwLTraw.clear();       
  Tau_byIsolationMVArun2v1DBoldDMwLTraw.clear();       
  Tau_byIsolationMVArun2v1PWdR03oldDMwLTraw.clear();   

  Tau_againstMuonLoose3.clear();
  Tau_againstMuonTight3.clear();
  Tau_againstElectronMVAVLooseMVA6.clear();
  Tau_againstElectronMVALooseMVA6.clear();
  Tau_againstElectronMVAMediumMVA6.clear();
  Tau_againstElectronMVATightMVA6.clear();
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
  Tau_signal_cand_chargedHadr_pt.clear();
  Tau_signal_cand_chargedHadr_eta.clear();
  Tau_signal_cand_chargedHadr_phi.clear();
  Tau_signal_cand_chargedHadr_energy.clear();
  Tau_signal_cand_neutHadr_pt.clear();
  Tau_signal_cand_neutHadr_eta.clear();
  Tau_signal_cand_neutHadr_phi.clear();
  Tau_signal_cand_neutHadr_energy.clear();
}

bool TauSelector::isGoodVertex(const reco::Vertex& vtxxx) {
  if (vtxxx.isFake()) return false;
  if (vtxxx.ndof() < _Tau_vtx_ndof_min) return false;
  if (vtxxx.position().Rho() > _Tau_vtx_rho_max) return false;
  if (fabs(vtxxx.position().Z()) > _Tau_vtx_position_z_max) return false;
  return true;
}

