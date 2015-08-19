// 
//  Authors:  Andres Florez: Universidad de los Andes, Colombia. 
//  kaur amandeepkalsi: Panjab University, India. 
//

#ifndef __TAU_MU_H_                                                                                                                                  
#define __TAU_MU_H_

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
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
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

  TauSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
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
  vector <int> Tau_byLooseCombinedIsolationDeltaBetaCorr,Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;
  vector <int> Tau_byMediumCombinedIsolationDeltaBetaCorr,Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits; 
  vector <int> Tau_byTightCombinedIsolationDeltaBetaCorr,Tau_byTightCombinedIsolationDeltaBetaCorr3Hits; 
  vector <int> Tau_byLooseIsolationMVA3newDMwLT,Tau_byLooseIsolationMVA3newDMwoLT,Tau_byMediumIsolationMVA3newDMwLT; 
  vector <int> Tau_byMediumIsolationMVA3newDMwoLT, Tau_byLooseIsolationMva3oldDMwLT, Tau_byLooseIsolationMVA3oldDMwoLT; 
  vector <int> Tau_byMediumIsolationMva3oldDMwLT, Tau_byMediumIsolationMVA3oldDMwoLT, Tau_byTightIsolationMVA3newDMwLT;  
  vector <int> Tau_byTightIsolationMVA3newDMwoLT, Tau_byTightIsolationMva3oldDMwLT, Tau_byTightIsolationMVA3oldDMwoLT; 
  vector <int> Tau_againstMuonLoose2, Tau_againstMuonLoose3, Tau_againstMuonTight2, Tau_againstMuonTight3; 
  vector <int> Tau_againstElectronMVALooseMVA5, Tau_againstElectronMVAMediumMVA5, Tau_againstElectronMVATightMVA5; 
  vector <int> Tau_byVLooseCombinedIsolationDeltaBetaCorr; 

  // config variables
  edm::InputTag tauToken_;
  edm::InputTag _vertexInputTag;
  edm::InputTag _beamSpot;

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
