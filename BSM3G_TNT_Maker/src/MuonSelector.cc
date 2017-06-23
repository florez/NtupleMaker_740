// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#include "NtupleMaker/BSM3G_TNT_Maker/interface/MuonSelector.h"

using namespace edm;
using namespace std;
using namespace reco;

MuonSelector::MuonSelector(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector&& iCC):baseTree(name,tree,debug){
  if(debug) std::cout << "BSM3G TNT Maker: In the MuonSelector Constructor --> getting parameters & calling SetBranches()." << std::endl;
  _muonToken               = iCC.consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"));
  _vertexInputTag          = iCC.consumes<reco::VertexCollection >(iConfig.getParameter<edm::InputTag>("vertices"));
  _beamSpot                = iCC.consumes<reco::BeamSpot >(iConfig.getParameter<edm::InputTag>("beamSpot"));
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

void MuonSelector::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();

  if(debug_) std::cout << "     MuonSelector: Cleared the vectors, grabbing the muon/PV/BeamSpot collection handles." << std::endl;

  int mucoun = 0;

  // grab handle to the muon collection
  edm::Handle<edm::View<pat::Muon> > muon_h;
  iEvent.getByToken(_muonToken, muon_h);

  if(debug_) std::cout << "     MuonSelector: Looping over PVs to extract the position of the best PV." << std::endl;

  // grab handle to the vertex collection
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByToken(_vertexInputTag, vtx_h);
  reco::VertexCollection::const_iterator firstGoodVertex = vtx_h->end();
  for (reco::VertexCollection::const_iterator it = vtx_h->begin(); it != vtx_h->end(); it++) {
    if (isGoodVertex(*it)) {
      firstGoodVertex = it;
      break;
    }
  }
  if (firstGoodVertex == vtx_h->end()) return; // skip event if there are no good PVs
  GlobalPoint thepv(firstGoodVertex->position().x(),firstGoodVertex->position().y(),firstGoodVertex->position().z());

  if(debug_) std::cout << "     MuonSelector: Extracting the beamspot position." << std::endl;

  // Get BeamSpot information
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(_beamSpot, beamSpotHandle);
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

  if(debug_) std::cout << "     MuonSelector: Looping over muons." << std::endl;

  // loop over muons
  for(edm::View<pat::Muon>::const_iterator mu = muon_h->begin(); mu != muon_h->end(); mu++) {
    if (mu->pt() < _Muon_pt_min) continue;
    if (fabs(mu->eta()) > _Muon_eta_max) continue;
    if (!(mu->innerTrack().isNonnull())) continue;
    if (!(mu->globalTrack().isNonnull())) continue;

    reco::TrackRef gtk = mu->globalTrack();

    // fill root tree with "necessary" information:  kinematics, ID, isolation
    Muon_pt.push_back(mu->pt());
    Muon_eta.push_back(mu->eta());
    Muon_phi.push_back(mu->phi());
    Muon_energy.push_back(mu->energy());
    Muon_charge.push_back(mu->charge());
    Muon_loose.push_back(mu->isLooseMuon());
    Muon_medium.push_back(mu->isMediumMuon());
    Muon_tight.push_back(mu->isTightMuon(*firstGoodVertex));
    Muon_soft.push_back(mu->isSoftMuon(*firstGoodVertex));
    Muon_isHighPt.push_back(mu->isHighPtMuon(*firstGoodVertex));
    Muon_pf.push_back(mu->isPFMuon());
    Muon_isTrackerMuon.push_back(mu->isTrackerMuon());
    Muon_POGisGood.push_back(muon::isGoodMuon(*mu, muon::TMOneStationTight));
    Muon_isoSum.push_back((mu->trackIso() + mu->ecalIso() + mu->hcalIso()));
    Muon_isoCharParPt.push_back((mu->pfIsolationR04().sumChargedParticlePt));
    Muon_isoCharged.push_back((mu->pfIsolationR04().sumChargedHadronPt));
    Muon_isoNeutralHadron.push_back((mu->pfIsolationR04().sumNeutralHadronEt));
    Muon_isoPhoton.push_back((mu->pfIsolationR04().sumPhotonEt));
    Muon_isoPU.push_back((mu->pfIsolationR04().sumPUPt));
    // PF combined iso with DeltaBeta corrections
    double combined_iso = (mu->pfIsolationR04().sumChargedHadronPt + max(0., mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - 0.5*mu->pfIsolationR04().sumPUPt))/mu->pt();
    Muon_combinedIso.push_back(combined_iso);
    Muon_trackRe_iso.push_back(mu->isolationR03().sumPt/mu->pt());

    // store additional information such as lifetime variables (needed for tau analyses)
    if (!_super_TNT){
      Muon_chi2.push_back(gtk->normalizedChi2());
      Muon_validHits.push_back(gtk->hitPattern().numberOfValidMuonHits());
      Muon_matchedStat.push_back(mu->numberOfMatchedStations());
      Muon_isGlobal.push_back(mu->isGlobalMuon());
      Muon_validHitsInner.push_back(mu->innerTrack()->hitPattern().numberOfValidPixelHits());
      Muon_TLayers.push_back(mu->innerTrack()->hitPattern().trackerLayersWithMeasurement());
      Muon_dxy_pv.push_back(mu->innerTrack()->dxy(firstGoodVertex->position()));
      Muon_dxy_bs.push_back(-1.*(mu->innerTrack()->dxy(point)));
      Muon_dz_pv.push_back(mu->innerTrack()->dz(firstGoodVertex->position()));
      Muon_dz_bs.push_back(mu->innerTrack()->dz(point));
      Muon_dzError.push_back(mu->innerTrack()->dzError());
      Muon_dxyError.push_back(mu->innerTrack()->d0Error());
      Muon_ndof.push_back(mu->innerTrack()->ndof());
      Muon_vtx.push_back(mu->innerTrack()->vx());
      Muon_vty.push_back(mu->innerTrack()->vy());
      Muon_vtz.push_back(mu->innerTrack()->vz());
      Muon_track_pt.push_back(mu->innerTrack()->pt());
      Muon_track_ptError.push_back(mu->innerTrack()->ptError());
      Muon_dB.push_back(mu->dB());
      Muon_besttrack_pt.push_back(mu->muonBestTrack()->pt());
      Muon_besttrack_ptError.push_back(mu->muonBestTrack()->ptError());
      reco::TrackRef tunePBestTrack = mu->tunePMuonBestTrack();
      Muon_tunePBestTrack_pt.push_back(tunePBestTrack.isNonnull() ? tunePBestTrack->pt() : -1.);
      reco::Muon::MuonTrackType tunePBestTrackType = mu->tunePMuonBestTrackType();
      Muon_tunePBestTrackType.push_back(tunePBestTrackType);
      Muon_chi2LocalPosition.push_back(mu->combinedQuality().chi2LocalPosition);
      Muon_trkKink.push_back(mu->combinedQuality().trkKink);
      Muon_segmentCompatibility.push_back(mu->segmentCompatibility());
      Muon_validFraction.push_back(mu->innerTrack()->validFraction());
      Muon_pixelLayersWithMeasurement.push_back(mu->innerTrack()->hitPattern().pixelLayersWithMeasurement());
      Muon_qualityhighPurity.push_back(mu->innerTrack()->quality(reco::TrackBase::highPurity));

      // point of closest approach (PCA) to the beamspot and primary vertex
      TransientTrack muonTransTkPtr = theB->build(*(mu->innerTrack()));
      GlobalPoint mu_pca_bs = muonTransTkPtr.trajectoryStateClosestToPoint(thebs).position();
      GlobalPoint mu_pca_pv = muonTransTkPtr.trajectoryStateClosestToPoint(thepv).position();
      Muon_track_PCAx_bs.push_back(mu_pca_bs.x());
      Muon_track_PCAy_bs.push_back(mu_pca_bs.y());
      Muon_track_PCAz_bs.push_back(mu_pca_bs.z());
      Muon_track_PCAx_pv.push_back(mu_pca_pv.x());
      Muon_track_PCAy_pv.push_back(mu_pca_pv.y());
      Muon_track_PCAz_pv.push_back(mu_pca_pv.z());

      // extract track fit errors
      const float muonMass = 0.1056583715;
      float muonSigma = muonMass*1e-6;
      float chi2 = 0.0;
      float ndf = 0.0;
      KinematicParticleFactoryFromTransientTrack pFactory;
      RefCountedKinematicParticle muonParticle = pFactory.particle(muonTransTkPtr, muonMass, chi2, ndf, muonSigma);
      Muon_trackFitErrorMatrix_00.push_back(muonParticle->stateAtPoint(mu_pca_bs).kinematicParametersError().matrix()(0,0));
      Muon_trackFitErrorMatrix_01.push_back(muonParticle->stateAtPoint(mu_pca_bs).kinematicParametersError().matrix()(0,1));
      Muon_trackFitErrorMatrix_02.push_back(muonParticle->stateAtPoint(mu_pca_bs).kinematicParametersError().matrix()(0,2));
      Muon_trackFitErrorMatrix_11.push_back(muonParticle->stateAtPoint(mu_pca_bs).kinematicParametersError().matrix()(1,1));
      Muon_trackFitErrorMatrix_12.push_back(muonParticle->stateAtPoint(mu_pca_bs).kinematicParametersError().matrix()(1,2));
      Muon_trackFitErrorMatrix_22.push_back(muonParticle->stateAtPoint(mu_pca_bs).kinematicParametersError().matrix()(2,2));
    }

    mucoun++;
  }

}

void MuonSelector::SetBranches(){
  if(debug_) std::cout << "     MuonSelector: Setting branches by calling AddBranch of baseTree." << std::endl;

  AddBranch(&Muon_pt                ,"Muon_pt");
  AddBranch(&Muon_eta               ,"Muon_eta");
  AddBranch(&Muon_phi               ,"Muon_phi");
  AddBranch(&Muon_energy            ,"Muon_energy");
  AddBranch(&Muon_charge            ,"Muon_charge");
  AddBranch(&Muon_tight             ,"Muon_tight");
  AddBranch(&Muon_loose             ,"Muon_loose");
  AddBranch(&Muon_medium            ,"Muon_medium");
  AddBranch(&Muon_soft              ,"Muon_soft");
  AddBranch(&Muon_isHighPt          ,"Muon_isHighPt");
  AddBranch(&Muon_pf                ,"Muon_pf");
  AddBranch(&Muon_isoCharged        ,"Muon_isoCharged");
  AddBranch(&Muon_isoSum            ,"Muon_isoSum");
  AddBranch(&Muon_isoCharParPt      ,"Muon_isoCharParPt");
  AddBranch(&Muon_isTrackerMuon    ,"Muon_isTrackerMuon");
  AddBranch(&Muon_POGisGood        ,"Muon_POGisGood");

  if (!_super_TNT){
    AddBranch(&Muon_chi2              ,"Muon_chi2");
    AddBranch(&Muon_validHits         ,"Muon_validHits");
    AddBranch(&Muon_validHitsInner    ,"Muon_validHitsInner");
    AddBranch(&Muon_matchedStat       ,"Muon_matchedStat");
    AddBranch(&Muon_dxy_pv            ,"Muon_dxy_pv");
    AddBranch(&Muon_dxy_bs            ,"Muon_dxy_bs");
    AddBranch(&Muon_dz_pv             ,"Muon_dz_pv");
    AddBranch(&Muon_dz_bs             ,"Muon_dz_bs");
    AddBranch(&Muon_dxyError           ,"Muon_dxyError");
    AddBranch(&Muon_dzError           ,"Muon_dzError");
    AddBranch(&Muon_ndof              ,"Muon_ndof");
    AddBranch(&Muon_vtx               ,"Muon_vtx");
    AddBranch(&Muon_vty               ,"Muon_vty");
    AddBranch(&Muon_vtz               ,"Muon_vtz");
    AddBranch(&Muon_track_pt          ,"Muon_track_pt");
    AddBranch(&Muon_track_ptError     ,"Muon_track_ptError");
    AddBranch(&Muon_isGlobal          ,"Muon_isGlobal");
    AddBranch(&Muon_TLayers           ,"Muon_TLayers");
    AddBranch(&Muon_isoNeutralHadron  ,"Muon_isoNeutralHadron");
    AddBranch(&Muon_isoPhoton         ,"Muon_isoPhoton");
    AddBranch(&Muon_isoPU             ,"Muon_isoPU");
    AddBranch(&Muon_combinedIso       ,"Muon_combinedIso");
    AddBranch(&Muon_trackRe_iso       ,"Muon_trackRe_iso");
    AddBranch(&Muon_dB               ,"Muon_dB");
    AddBranch(&Muon_besttrack_pt     ,"Muon_besttrack_pt");
    AddBranch(&Muon_besttrack_ptError ,"Muon_besttrack_ptError");
    AddBranch(&Muon_tunePBestTrack_pt      ,"Muon_tunePBestTrack_pt");
    AddBranch(&Muon_tunePBestTrackType         ,"Muon_tunePBestTrackType");
    AddBranch(&Muon_chi2LocalPosition          ,"Muon_chi2LocalPosition");
    AddBranch(&Muon_trkKink                    ,"Muon_trkKink");
    AddBranch(&Muon_segmentCompatibility       ,"Muon_segmentCompatibility");
    AddBranch(&Muon_validFraction              ,"Muon_validFraction");
    AddBranch(&Muon_pixelLayersWithMeasurement ,"Muon_pixelLayersWithMeasurement");
    AddBranch(&Muon_qualityhighPurity          ,"Muon_qualityhighPurity");
    AddBranch(&Muon_track_PCAx_bs              ,"Muon_track_PCAx_bs");
    AddBranch(&Muon_track_PCAy_bs              ,"Muon_track_PCAy_bs");
    AddBranch(&Muon_track_PCAz_bs              ,"Muon_track_PCAz_bs");
    AddBranch(&Muon_track_PCAx_pv              ,"Muon_track_PCAx_pv");
    AddBranch(&Muon_track_PCAy_pv              ,"Muon_track_PCAy_pv");
    AddBranch(&Muon_track_PCAz_pv              ,"Muon_track_PCAz_pv");
    AddBranch(&Muon_trackFitErrorMatrix_00     ,"Muon_trackFitErrorMatrix_00");
    AddBranch(&Muon_trackFitErrorMatrix_01     ,"Muon_trackFitErrorMatrix_01");
    AddBranch(&Muon_trackFitErrorMatrix_02     ,"Muon_trackFitErrorMatrix_02");
    AddBranch(&Muon_trackFitErrorMatrix_11     ,"Muon_trackFitErrorMatrix_11");
    AddBranch(&Muon_trackFitErrorMatrix_12     ,"Muon_trackFitErrorMatrix_12");
    AddBranch(&Muon_trackFitErrorMatrix_22     ,"Muon_trackFitErrorMatrix_22");
  }

  if(debug_) std::cout << "     MuonSelector: Finished setting branches." << std::endl;
}

void MuonSelector::Clear(){

  Muon_pt.clear();
  Muon_eta.clear();
  Muon_phi.clear();
  Muon_energy.clear();
  Muon_tight.clear();
  Muon_loose.clear();
  Muon_medium.clear();
  Muon_soft.clear();
  Muon_pf.clear();
  Muon_isHighPt.clear();
  Muon_isoCharged.clear();
  Muon_isoSum.clear();
  Muon_isoCharParPt.clear();
  Muon_isoNeutralHadron.clear();
  Muon_isoPhoton.clear();
  Muon_isoPU.clear();
  Muon_combinedIso.clear();
  Muon_trackRe_iso.clear();
  Muon_charge.clear();
  Muon_chi2.clear();
  Muon_validHits.clear();
  Muon_validHitsInner.clear();
  Muon_matchedStat.clear();
  Muon_dxy_pv.clear();
  Muon_dxy_bs.clear();
  Muon_dz_pv.clear();
  Muon_dz_bs.clear();
  Muon_dxyError.clear();
  Muon_dzError.clear();
  Muon_ndof.clear();
  Muon_vtx.clear();
  Muon_vty.clear();
  Muon_vtz.clear();
  Muon_isGlobal.clear();
  Muon_track_pt.clear();
  Muon_track_ptError.clear();
  Muon_TLayers.clear();
  Muon_dB.clear();
  Muon_besttrack_pt.clear();
  Muon_besttrack_ptError.clear();
  Muon_tunePBestTrack_pt.clear();
  Muon_isTrackerMuon.clear();
  Muon_POGisGood.clear();
  Muon_chi2LocalPosition.clear();
  Muon_trkKink.clear();
  Muon_segmentCompatibility.clear();
  Muon_validFraction.clear();
  Muon_pixelLayersWithMeasurement.clear();
  Muon_qualityhighPurity.clear();
  Muon_tunePBestTrackType.clear();
  Muon_track_PCAx_bs.clear();
  Muon_track_PCAy_bs.clear();
  Muon_track_PCAz_bs.clear();
  Muon_track_PCAx_pv.clear();
  Muon_track_PCAy_pv.clear();
  Muon_track_PCAz_pv.clear();
  Muon_trackFitErrorMatrix_00.clear();
  Muon_trackFitErrorMatrix_01.clear();
  Muon_trackFitErrorMatrix_02.clear();
  Muon_trackFitErrorMatrix_11.clear();
  Muon_trackFitErrorMatrix_12.clear();
  Muon_trackFitErrorMatrix_22.clear();
}

bool MuonSelector::isGoodVertex(const reco::Vertex& vtx) {
  if (vtx.isFake()) return false;
  if (vtx.ndof() < _Muon_vtx_ndof_min) return false;
  if (vtx.position().Rho() > _Muon_vtx_rho_max) return false;
  if (fabs(vtx.position().Z()) > _Muon_vtx_position_z_max) return false;
  return true;
}
