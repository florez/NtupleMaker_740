//
//  Authors:  Andres Florez: Universidad de los Andes, Colombia.
//  kaur amandeepkalsi: Panjab University, India.
//

#ifndef __TAU_MU_H_
#define __TAU_MU_H_

//#include <memory>
//#include <iostream>
//#include <cmath>
//#include <vector>
//#include <TBranch.h>
//#include <TTree.h>
//#include <TFile.h>
//#include <TH1.h>
//#include <TH2.h>
//#include <string>
//#include <map>
//#include <sstream>
//#include <stdio.h>
//#include <stdlib.h>
//#include <TRandom3.h>
//#include <TBranch.h>
//#include <TClonesArray.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
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

class TauSelector : public baseTree{

public:

  TauSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector&& icc);
  ~TauSelector();
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void SetBranches();
  void Clear();
  bool isGoodVertex(const reco::Vertex& vtxxx);

private:
  TauSelector(){};

  vector <double> Tau_eta, Tau_phi, Tau_pt, Tau_energy, Tau_charge, Tau_chargedIsoPtSum, Tau_neutralIsoPtSum, Tau_puCorrPtSum ;
  vector<double>  Tau_leadChargedCandPt, Tau_leadChargedCandCharge, Tau_leadChargedCandEta, Tau_leadChargedCandPhi,Tau_nProngs;
  vector<double> Tau_leadChargedCandChi2,  Tau_leadChargedCandValidHits,  Tau_leadChargedCandDxy_pv,  Tau_leadChargedCandDxy_bs,  Tau_leadChargedCandDz_bs;
  vector<double> Tau_leadChargedCandDz_pv,  Tau_leadChargedCandDzError,  Tau_leadChargedCandDxyError,  Tau_leadChargedCandNdof,  Tau_leadChargedCandVtx;
  vector<double> Tau_leadChargedCandVty,  Tau_leadChargedCandVtz,  Tau_leadChargedCandTrack_pt,  Tau_leadChargedCandTrack_ptError;
  vector<double> Tau_leadChargedCandTrack_PCAx_bs, Tau_leadChargedCandTrack_PCAy_bs, Tau_leadChargedCandTrack_PCAz_bs;
  vector<double> Tau_leadChargedCandTrack_PCAx_pv, Tau_leadChargedCandTrack_PCAy_pv, Tau_leadChargedCandTrack_PCAz_pv;
  vector<double> Tau_leadChargedCandTrackFitErrorMatrix_00, Tau_leadChargedCandTrackFitErrorMatrix_01, Tau_leadChargedCandTrackFitErrorMatrix_02, Tau_leadChargedCandTrackFitErrorMatrix_11, Tau_leadChargedCandTrackFitErrorMatrix_12, Tau_leadChargedCandTrackFitErrorMatrix_22;
  vector<double> Tau_defaultDxy, Tau_defaultDxyError, Tau_defaultDxySig, Tau_defaultFlightLengthSig, Tau_default_PCAx_pv, Tau_default_PCAy_pv, Tau_default_PCAz_pv;
  vector<double> Tau_defaultFlightLengthX, Tau_defaultFlightLengthY, Tau_defaultFlightLengthZ;
  vector <int> Tau_decayModeFinding, Tau_decayModeFindingOldDMs, Tau_decayModeFindingNewDMs;
  vector <int> Tau_byVLooseIsolationMVArun2v1DBoldDMwLT, Tau_byLooseIsolationMVArun2v1DBoldDMwLT, Tau_byMediumIsolationMVArun2v1DBoldDMwLT;
  vector <int> Tau_byTightIsolationMVArun2v1DBoldDMwLT, Tau_byVTightIsolationMVArun2v1DBoldDMwLT;
  vector <int> Tau_byVLooseIsolationMVArun2v1DBnewDMwLT, Tau_byLooseIsolationMVArun2v1DBnewDMwLT, Tau_byMediumIsolationMVArun2v1DBnewDMwLT;
  vector <int> Tau_byTightIsolationMVArun2v1DBnewDMwLT, Tau_byVTightIsolationMVArun2v1DBnewDMwLT;

  vector <int> Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits, Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits, Tau_byTightCombinedIsolationDeltaBetaCorr3Hits;
  vector <int> Tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03, Tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03, Tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03;

  vector <int> Tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT, Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT, Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
  vector <int> Tau_byTightIsolationMVArun2v1DBdR03oldDMwLT, Tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT;
  vector <int> Tau_againstMuonLoose2, Tau_againstMuonLoose3, Tau_againstMuonTight2, Tau_againstMuonTight3;
  vector <int> Tau_againstElectronMVALooseMVA6, Tau_againstElectronMVAMediumMVA6, Tau_againstElectronMVATightMVA6, Tau_againstElectronMVAVLooseMVA6;
  vector <int> Tau_byVLooseCombinedIsolationDeltaBetaCorr, Tau_byVLooseIsolationMVA3newDMwLT, Tau_byVLooseIsolationMva3oldDMwLT;
  vector <int> Tau_decayMode;
  vector<double> Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits;
  vector<double> Tau_byIsolationMVArun2v1DBdR03oldDMwLTraw;
  vector<double> Tau_byIsolationMVArun2v1DBnewDMwLTraw;
  vector<double> Tau_byIsolationMVArun2v1DBoldDMwLTraw;
  vector<double> Tau_byIsolationMVArun2v1PWdR03oldDMwLTraw;
  
  vector< vector< float > > Tau_signal_cand_chargedHadr_pt,  Tau_signal_cand_chargedHadr_eta,  Tau_signal_cand_chargedHadr_phi,  Tau_signal_cand_chargedHadr_energy;
  vector< vector< float > > Tau_signal_cand_neutGam_pt,  Tau_signal_cand_neutGam_eta,  Tau_signal_cand_neutGam_phi,  Tau_signal_cand_neutGam_energy;
  vector< vector< float > > Tau_signal_cand_neutHadr_pt,  Tau_signal_cand_neutHadr_eta,  Tau_signal_cand_neutHadr_phi,  Tau_signal_cand_neutHadr_energy;


  // config variables
//  edm::InputTag tauToken_;
//  edm::InputTag packedPFCandidateToken_;
//  edm::InputTag _vertexInputTag;
//  edm::InputTag _beamSpot;
  edm::EDGetTokenT<edm::View<pat::Tau> > tauToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> packedPFCandidateToken_;
  edm::EDGetTokenT<reco::VertexCollection> _vertexInputTag;
  edm::EDGetTokenT<reco::BeamSpot> _beamSpot;

  double _Tau_pt_min;
  double _Tau_eta_max;
  int _Tau_vtx_ndof_min;
  int _Tau_vtx_rho_max;
  double _Tau_vtx_position_z_max;
  bool   _super_TNT;
};

#endif

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
