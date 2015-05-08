// 
//  Authors:  Andres Florez: Universidad de los Andes, Colombia. 
//  kaur amandeepkalsi: Panjab University, India. 
// 

#ifndef __MUON_MU_H_                                                                                                                                  
#define __MUON_MU_H_

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
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "baseTree.h"

using namespace std;
using namespace pat;
using namespace edm;
//
// class declaration
//

class MuonSelector : public  baseTree{

public:

  MuonSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg);
  ~MuonSelector();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();
  bool isGoodVertex(const reco::Vertex& vtx);
private:
  MuonSelector(){};
  
  // ----------member data ---------------------------

  vector <double> Muon_isoSum , Muon_isoCharParPt ;
  vector <double> Muon_pt ,Muon_eta,Muon_phi, Muon_dz, Muon_energy, Muon_iso;
  vector <double> Muon_isoCharged, Muon_isoNeutralHadron , Muon_isoPhoton, Muon_isoPU;
  vector <double> Muon_charge, Muon_chi2, Muon_p, Muon_matchedStat, Muon_dxy, Muon_validHits, Muon_validHitsInner, Muon_TLayers; 
  
  vector <bool> Muon_tight, Muon_soft, Muon_pf ;   

  // super tiny ntuple?
  bool _super_TNT;

  // confit variables
  edm::InputTag _muonToken;
  edm::InputTag _vertexInputTag; 
  double _Muon_pt_min;
  double _Muon_eta_max;
  int _Muon_vtx_ndof_min;
  int _Muon_vtx_rho_max;
  double _Muon_vtx_position_z_max;
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
