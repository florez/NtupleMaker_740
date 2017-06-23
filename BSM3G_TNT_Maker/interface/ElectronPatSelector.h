// Authors:  Alfredo Gurrola (Vanderbilt University)
// Andres Florez: Universidad de los Andes, Colombia.
// kaur amandeepkalsi: Panjab University, India.

#ifndef __ELECTRON_PAT_H_
#define __ELECTRON_PAT_H_


#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <string>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TBranch.h>
#include <TClonesArray.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "baseTree.h"

using namespace std;
using namespace pat;
using namespace edm;
//
// class declaration
//

class ElectronPatSelector : public  baseTree{

public:

  ElectronPatSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector && iC);
  ~ElectronPatSelector();
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void SetBranches();
  void Clear();
  bool isGoodVertex(const reco::Vertex& vtxx);
 private:
  ElectronPatSelector(){};

  edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronMVAwp1Token_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronMVAwp2Token_;
  // ----------member data ---------------------------

  vector<double> patElectron_pt , patElectron_eta, patElectron_phi, patElectron_energy, patElectron_charge;
  vector<double> patElectron_gsfTrack_ndof, patElectron_gsfTrack_dxy_pv, patElectron_dxyError, patElectron_gsfTrack_normChi2, patElectron_gsfTrack_dz_pv, patElectron_gsfTrack_dz_bs;
  vector<double> patElectron_gsfTrack_vtx, patElectron_gsfTrack_vty, patElectron_gsfTrack_vtz;
  vector<double> patElectron_gsfTrack_dxy_bs, isoChargedHadrons_, isoNeutralHadrons_, isoPhotons_, isoPU_;
  vector<double> patElectron_gsfTrack_PCAx_bs, patElectron_gsfTrack_PCAy_bs, patElectron_gsfTrack_PCAz_bs;
  vector<double> patElectron_gsfTrack_PCAx_pv, patElectron_gsfTrack_PCAy_pv, patElectron_gsfTrack_PCAz_pv;
  vector<int>  passVetoId_, passLooseId_, passMediumId_, passTightId_, passHEEPId_, passMVAwp1Id_, passMVAwp2Id_, passConversionVeto_, expectedMissingInnerHits;
  vector<double> patElectron_gsfTrackFitErrorMatrix_00, patElectron_gsfTrackFitErrorMatrix_01, patElectron_gsfTrackFitErrorMatrix_02, patElectron_gsfTrackFitErrorMatrix_11, patElectron_gsfTrackFitErrorMatrix_12, patElectron_gsfTrackFitErrorMatrix_22;

  // super tiny ntuple?
  bool _super_TNT;

  // confit variables
//  edm::InputTag _patElectronToken;
  edm::EDGetTokenT<edm::View<pat::Electron> > _patElectronToken;
//  edm::InputTag _vertexInputTag;
//  edm::InputTag _beamSpot;
  edm::EDGetTokenT<reco::VertexCollection> _vertexInputTag;
  edm::EDGetTokenT<reco::BeamSpot> _beamSpot;
  double _patElectron_pt_min;
  double _patElectron_eta_max;
  int _patElectron_vtx_ndof_min;
  int _patElectron_vtx_rho_max;
  double _patElectron_vtx_position_z_max;

};

#endif
