// 
// Authors:  Andres Florez: Universidad de los Andes, Colombia. 
// kaur amandeepkalsi: Panjab University, India. 
//
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
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "Math/VectorUtil.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/Math/interface/deltaR.h"
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

class ElectronPatSelector : public  baseTree{

public:

  ElectronPatSelector(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector && iC);
  ~ElectronPatSelector();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();
 private:
  ElectronPatSelector(){};
  
  edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;
  // ----------member data ---------------------------
  
  vector<double> patElectron_pt , patElectron_eta, patElectron_phi, patElectron_energy, patElectron_charge;
  vector<double> patElectron_d0, patElectron_dz, isoChargedHadrons_, isoNeutralHadrons_, isoPhotons_;
  vector<int> passVetoId_, passLooseId_, passMediumId_, passTightId_, passConversionVeto_, expectedMissingInnerHits;
  
  // confit variables
  edm::InputTag _patElectronToken;
  edm::InputTag _vertexInputTag;
  double _patElectron_pt_min;
  double _patElectron_eta_max;
};

#endif
